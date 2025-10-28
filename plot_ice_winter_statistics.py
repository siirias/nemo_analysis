#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
from typing import Optional, Dict, Any

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from mpl_toolkits.axes_grid1.inset_locator import inset_axes


# ---------------- helpers for lon/lat healing ----------------

def _fill_1d_nearest(a: np.ndarray) -> np.ndarray:
    """Fill NaNs in a 1D array by nearest-value (edge-hold) via index interpolation."""
    a = a.astype(float)
    good = np.isfinite(a)
    if not good.any():
        # Fallback to a simple ramp to avoid NaNs in edges
        return np.arange(a.size, dtype=float)
    idx = np.arange(a.size, dtype=float)
    return np.interp(idx, idx[good], a[good])


def _centers_to_edges(c: np.ndarray) -> np.ndarray:
    """Convert 1D center coordinates to strictly monotonic edges for pcolormesh."""
    c = np.asarray(c, dtype=float)
    d = np.diff(c)
    if not np.isfinite(d).any():
        d = np.ones_like(c[:-1])
    else:
        med = np.nanmedian(d[np.isfinite(d)])
        d = np.where(np.isfinite(d), d, med)
    mids = 0.5 * (c[:-1] + c[1:])
    left  = c[0]  - 0.5 * d[0]
    right = c[-1] + 0.5 * d[-1]
    edges = np.concatenate(([left], mids, [right]))
    # enforce strictly increasing (guard against tiny numerical wobbles)
    eps = np.finfo(float).eps * max(1.0, np.nanmax(np.abs(edges)))
    for k in range(1, edges.size):
        if edges[k] <= edges[k - 1]:
            edges[k] = edges[k - 1] + eps
    return edges


# ---------------- data preparation ----------------

def prepare_ice_field(nc_file: str,
                      bathy_file: str,
                      var_name: str = "ice_season_length") -> Dict[str, Any]:
    """
    Load and prepare a NEMO-style ice field for plotting.

    Steps:
      - open data and bathy
      - build sea mask from bathy (bdy_msk != 0)
      - mask variable on land
      - heal nav_lon/nav_lat (zeros -> NaN, derive lon(x)/lat(y), fill gaps)
      - compute 1D lon/lat edges for pcolormesh

    Returns a dict with:
      Z       : 2D np.ndarray (masked with NaN on land)
      lon_e   : 1D np.ndarray (X+1) pcolormesh edges
      lat_e   : 1D np.ndarray (Y+1) pcolormesh edges
      lon_1d  : 1D centers along x
      lat_1d  : 1D centers along y
      extent  : (lon_min, lon_max, lat_min, lat_max)
      meta    : dict with central_lon, central_lat, var_name, units
    """
    ds = xr.open_dataset(nc_file)
    bathy = xr.open_dataset(bathy_file)

    if var_name not in ds:
        raise KeyError(f"Variable '{var_name}' not in {nc_file}. "
                       f"Available: {list(ds.data_vars)}")

    Z_da = ds[var_name]  # (y,x)

    # Align bathy mask to data grid
    msk = xr.DataArray(np.array(bathy["bdy_msk"]), dims=Z_da.dims).reindex_like(Z_da)
    sea = (msk != 0)
    Z_da = Z_da.where(sea)  # NaN on land

    # Pull lon/lat and mark zeros as NaN
    lon2d = ds["nav_lon"].values.copy()
    lat2d = ds["nav_lat"].values.copy()
    lon2d[lon2d == 0] = np.nan
    lat2d[lat2d == 0] = np.nan

    Y, X = Z_da.shape

    # Recover 1D lon(x) & lat(y) robustly by column/row median
    lon_1d = np.full(X, np.nan)
    for i in range(X):
        col = lon2d[:, i]
        if np.isfinite(col).any():
            lon_1d[i] = np.nanmedian(col)

    lat_1d = np.full(Y, np.nan)
    for j in range(Y):
        row = lat2d[j, :]
        if np.isfinite(row).any():
            lat_1d[j] = np.nanmedian(row)

    # Fill gaps by nearest across index
    lon_1d = _fill_1d_nearest(lon_1d)
    lat_1d = _fill_1d_nearest(lat_1d)

    # Edges for pcolormesh
    lon_e = _centers_to_edges(lon_1d)
    lat_e = _centers_to_edges(lat_1d)

    Z = Z_da.values
    extent = (float(lon_e.min()), float(lon_e.max()),
              float(lat_e.min()), float(lat_e.max()))

    units = Z_da.attrs.get("units", "")
    meta = dict(
        central_lon=float(np.nanmean(lon_1d)),
        central_lat=float(np.nanmean(lat_1d)),
        var_name=var_name,
        units=units
    )

    return dict(Z=Z, lon_e=lon_e, lat_e=lat_e,
                lon_1d=lon_1d, lat_1d=lat_1d,
                extent=extent, meta=meta)


# ---------------- plotting ----------------

def plot_ice_map(prep: Dict[str, Any],
                 output_path: Optional[str] = None,
                 projection: str = "albers",   # "albers" | "lambert" | "merc"
                 coastlines: bool = False,
                 vmin: float = 0.0, vmax: float = 200.0,
                 title: Optional[str] = None,
                 figsize=(8, 9), dpi: int = 300,
                 show: bool = True):
    """
    Plot prepared ice map with Cartopy. Saves to output_path if given.
    """
    Z = prep["Z"]
    lon_e = prep["lon_e"]
    lat_e = prep["lat_e"]
    extent = prep["extent"]
    meta = prep["meta"]

    # Projection
    if projection == "albers":
        proj_map = ccrs.AlbersEqualArea(
            central_longitude=meta["central_lon"],
            central_latitude=meta["central_lat"],
            standard_parallels=(60, 66),
        )
    elif projection == "lambert":
        proj_map = ccrs.LambertConformal(
            central_longitude=meta["central_lon"],
            central_latitude=meta["central_lat"],
        )
    elif projection == "merc":
        proj_map = ccrs.Mercator()
    else:
        raise ValueError("projection must be one of: 'albers', 'lambert', 'merc'")
    proj_data = ccrs.PlateCarree()

    # Colormap
    cmap = plt.cm.viridis.copy()
    cmap.set_bad("#6a3e25")  # land/background as brown

    fig, ax = plt.subplots(
        figsize=figsize,
        subplot_kw={"projection": proj_map},
        constrained_layout=True,
    )
    ax.set_facecolor("#6a3e25")
    ax.set_anchor("W")

    pc = ax.pcolormesh(
        lon_e, lat_e, Z,
        transform=proj_data, shading="flat",
        cmap=cmap, vmin=vmin, vmax=vmax
    )

    if coastlines:
        ax.add_feature(cfeature.COASTLINE.with_scale("10m"), linewidth=0.6)

    ax.set_extent(list(extent), crs=proj_data)

    # Gridlines: bottom longitudes, left latitudes; 1° steps; no inline labels
    gl = ax.gridlines(
        crs=proj_data, draw_labels=True,
        linewidth=0.5, alpha=0.4, linestyle=":",
        x_inline=False, y_inline=False
    )
    gl.right_labels   = False
    gl.top_labels     = False
    gl.left_labels    = True
    gl.bottom_labels  = True
    gl.rotate_labels  = False
    gl.xformatter     = LONGITUDE_FORMATTER
    gl.yformatter     = LATITUDE_FORMATTER
    gl.xlocator       = mticker.MultipleLocator(1)
    lat_min, lat_max  = float(lat_e.min()), float(lat_e.max())
    gl.ylocator       = mticker.FixedLocator(np.arange(np.floor(lat_min),
                                                       np.ceil(lat_max) + 1, 1))
    gl.xlabel_style   = {"size": 9}
    gl.ylabel_style   = {"size": 9}
    gl.xpadding       = 3
    gl.ypadding       = 3

    # Title
    if title is None:
        label_units = f" ({meta['units']})" if meta["units"] else ""
        title = meta["var_name"].replace("_", " ").title() + label_units
    ax.set_title(title, pad=12)

    # Inset colorbar (doesn't consume outer margins)
    cax = inset_axes(ax, width="3%", height="85%", loc="center right", borderpad=-1.0)
    cb = fig.colorbar(pc, cax=cax, orientation="vertical")
    cb.set_label(title)

    if output_path:
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        plt.savefig(output_path, bbox_inches="tight", dpi=dpi)

    if show:
        plt.show()
    else:
        plt.close(fig)

    return fig, ax


# ---------------- main ----------------

def main():
    # --- user-editables (paths & options) ---
    bathy_file = r"c:/Data/NemoTest/bathy_meter.nc"
    input_dir  = r"C:\Data\VanhataloEtAl\ice_seasons_2006to2100/"
    nc_file    = os.path.join(input_dir, "ice_season_2007-2008_set_B002.nc")
    var_name   = "ice_season_length"
    out_dir    = r"C:\Data\VanhataloEtAl\plots/"
    out_name   = f"{os.path.splitext(os.path.basename(nc_file))[0]}_{var_name}.png"
    out_path   = os.path.join(out_dir, out_name)

    prep = prepare_ice_field(nc_file, bathy_file, var_name=var_name)
    plot_ice_map(
        prep,
        output_path=out_path,
        projection="albers",     # 'albers'|'lambert'|'merc'
        coastlines=False,
        vmin=0.0, vmax=200.0,
        title="Ice season length",
        figsize=(5,5), dpi=300,
        show=True
    )


if __name__ == "__main__":
    main()
