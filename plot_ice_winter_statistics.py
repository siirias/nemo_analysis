#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import re
import glob
from matplotlib.animation import FuncAnimation, FFMpegWriter, PillowWriter

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
    cax = inset_axes(ax, width="3%", height="85%", loc="center right", borderpad=-2.0)
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

def _parse_set_years_from_name(fname: str):
    """
    Return (y1, y2, set_code) parsed from a filename like:
    ice_season_2007-2008_set_A002.nc
    """
    m = re.search(r"ice_season_(\d{4})-(\d{4})_set_([A-Z]\d{3})\.nc$", os.path.basename(fname))
    if not m:
        return None
    y1, y2, code = int(m.group(1)), int(m.group(2)), m.group(3)
    return y1, y2, code


def list_set_files(input_dir: str, set_code: str):
    """
    Find and sort all files of a given set (e.g., 'A002') by start year.
    """
    pattern = os.path.join(input_dir, f"ice_season_*_set_{set_code}.nc")
    files = [f for f in glob.glob(pattern) if _parse_set_years_from_name(f)]
    files.sort(key=lambda f: _parse_set_years_from_name(f)[0])  # sort by y1
    return files


def make_ice_animation(input_dir: str,
                       bathy_file: str,
                       set_code: str = "A002",
                       var_name: str = "ice_season_length",
                       projection: str = "albers",
                       coastlines: bool = False,
                       vmin: float = 0.0, vmax: float = 200.0,
                       figsize=(5, 5),
                       fps: int = 4,
                       out_dir: Optional[str] = None,
                       out_basename: Optional[str] = None,
                       writer: Optional[str] = None):
    """
    Build an animation over all 'ice_season_YYYY-YYYY_set_<set_code>.nc' files.
    Saves MP4 if ffmpeg writer is available, else GIF.
    """
    files = list_set_files(input_dir, set_code)
    if not files:
        raise FileNotFoundError(f"No files found for set {set_code} in {input_dir}")

    # Prepare first frame (defines grid/projection/extent)
    prep0 = prepare_ice_field(files[0], bathy_file, var_name=var_name)

    # Choose projection & create persistent figure/axes
    meta = prep0["meta"]
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

    cmap = plt.cm.viridis.copy()
    cmap.set_bad("#6a3e25")

    fig, ax = plt.subplots(figsize=figsize, subplot_kw={"projection": proj_map}, constrained_layout=True)
    ax.set_facecolor("#6a3e25")
    ax.set_anchor("W")
    ax.set_extent(list(prep0["extent"]), crs=proj_data)

    if coastlines:
        ax.add_feature(cfeature.COASTLINE.with_scale("10m"), linewidth=0.6)

    # Gridlines (same as your static map)
    gl = ax.gridlines(crs=proj_data, draw_labels=True, linewidth=0.5, alpha=0.4, linestyle=":",
                      x_inline=False, y_inline=False)
    gl.right_labels = False
    gl.top_labels = False
    gl.left_labels = True
    gl.bottom_labels = True
    gl.rotate_labels = False
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.xlocator = mticker.MultipleLocator(1)
    lat_min, lat_max = float(prep0["lat_e"].min()), float(prep0["lat_e"].max())
    gl.ylocator = mticker.FixedLocator(np.arange(np.floor(lat_min), np.ceil(lat_max) + 1, 1))
    gl.xlabel_style = {"size": 9}
    gl.ylabel_style = {"size": 9}
    gl.xpadding = 3
    gl.ypadding = 3

    # Persistent pcolormesh; we’ll only update its color array
    pc = ax.pcolormesh(prep0["lon_e"], prep0["lat_e"], prep0["Z"],
                       transform=proj_data, shading="flat",
                       cmap=cmap, vmin=vmin, vmax=vmax)

    # Inset colorbar
    from mpl_toolkits.axes_grid1.inset_locator import inset_axes
    cax = inset_axes(ax, width="3%", height="85%", loc="center right", borderpad=-1.0)
    cb = fig.colorbar(pc, cax=cax, orientation="vertical")
    title_base = f"{var_name.replace('_',' ')} — set {set_code}"
    cb.set_label("Ice season length (days)")
    ax.set_title(f"{title_base}", pad=12)

    # Preload Z fields and labels (faster & stable)
    Z_list = []
    labels = []
    for f in files:
        prep = prepare_ice_field(f, bathy_file, var_name=var_name)
        # sanity: grid shape must match initial
        if prep["Z"].shape != prep0["Z"].shape:
            raise ValueError(f"Grid shape differs in {f}: {prep['Z'].shape} vs {prep0['Z'].shape}")
        Z_list.append(prep["Z"])
        y1, y2, _code = _parse_set_years_from_name(f)
        labels.append(f"{y1}–{y2}")

    # Update function: set new color array and title
    def update(i):
        Zi = Z_list[i]
        # QuadMesh expects 1D array; mask NaNs:
        pc.set_array(np.ma.masked_invalid(Zi).ravel())
        ax.set_title(f"{title_base} — {labels[i]}", pad=12)
        return (pc,)

    anim = FuncAnimation(fig, update, frames=len(Z_list), interval=1000//fps, blit=False, repeat=True)

    # Choose writer
    if writer is None:
        try:
            _ = FFMpegWriter(fps=fps)  # probe
            writer = "ffmpeg"
        except Exception:
            writer = "pillow"

    if out_dir is None:
        out_dir = os.path.join(input_dir, "animations")
    os.makedirs(out_dir, exist_ok=True)

    if out_basename is None:
        out_basename = f"{set_code}_{var_name}_animation"

    if writer == "ffmpeg":
        out_path = os.path.join(out_dir, out_basename + ".mp4")
        anim.save(out_path, writer=FFMpegWriter(fps=fps, bitrate=3000))
    else:
        out_path = os.path.join(out_dir, out_basename + ".gif")
        anim.save(out_path, writer=PillowWriter(fps=fps))

    print(f"Saved animation to: {out_path}")
    return out_path


# ---------------- main ----------------

# def main():
#     # --- user-editables (paths & options) ---
#     bathy_file = r"c:/Data/NemoTest/bathy_meter.nc"
#     input_dir  = r"C:\Data\VanhataloEtAl\ice_seasons_2006to2100/"
#     nc_file    = os.path.join(input_dir, "ice_season_2007-2008_set_B002.nc")
#     var_name   = "ice_season_length"
#     out_dir    = r"C:\Data\VanhataloEtAl\plots/"
#     out_name   = f"{os.path.splitext(os.path.basename(nc_file))[0]}_{var_name}.png"
#     out_path   = os.path.join(out_dir, out_name)

#     prep = prepare_ice_field(nc_file, bathy_file, var_name=var_name)
#     plot_ice_map(
#         prep,
#         output_path=out_path,
#         projection="albers",     # 'albers'|'lambert'|'merc'
#         coastlines=False,
#         vmin=0.0, vmax=200.0,
#         title="Ice season length",
#         figsize=(5,5), dpi=300,
#         show=True
#     )

def main():
    bathy_file = r"c:/Data/NemoTest/bathy_meter.nc"
    input_dir  = r"C:\Data\VanhataloEtAl\ice_seasons_2006to2100/"
    var_name   = "ice_season_length"

    # one-off static example (kept, if you want)
    # nc_file = os.path.join(input_dir, "ice_season_2007-2008_set_A002.nc")
    # prep = prepare_ice_field(nc_file, bathy_file, var_name=var_name)
    # plot_ice_map(prep, output_path=None, projection="albers", coastlines=False,
    #              vmin=0.0, vmax=200.0, title="Ice season length", figsize=(5,5), show=True)

    # animation for a chosen set:
    make_ice_animation(
        input_dir=input_dir,
        bathy_file=bathy_file,
        set_code="B005",            # ← change to A005, B002, etc.
        var_name=var_name,
        projection="albers",
        coastlines=False,
        vmin=0.0, vmax=200.0,
        figsize=(5, 5),
        fps=4,                      # 4 frames per second
        out_dir=os.path.join(input_dir, "animations"),
        out_basename=None,          # auto name
        writer=None                 # auto choose ffmpeg→mp4, else Pillow→gif
    )




if __name__ == "__main__":
    main()
