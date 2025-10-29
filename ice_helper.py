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
import cmocean as cmo
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.colors import Colormap


    
# ---------------- helpers ----------------

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


def _var_label(meta: Dict[str, Any]) -> str:
    """Human-readable label from var name + units."""
    units = f" ({meta.get('units','')})" if meta.get('units') else ""
    return meta["var_name"].replace("_", " ").title() + units

def _resolve_cmap(cmap: Optional[Colormap or str] = None,
                  nan_color: str = "#6a3e25") -> Colormap:
    """
    Accepts a Matplotlib Colormap object or a string name (e.g. 'viridis', 'plasma', 'cividis', 'magma_r').
    Returns a copy with NaN color set to `nan_color`.
    """
    import matplotlib.pyplot as plt
    if cmap is None:
        cm = plt.get_cmap("viridis").copy()
    elif isinstance(cmap, str):
        cm = plt.get_cmap(cmap).copy()
    else:
        cm = cmap.copy() if hasattr(cmap, "copy") else cmap
    try:
        cm.set_bad(nan_color)
    except Exception:
        pass
    return cm

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

# ============== shared plotting factory ==============
def create_ice_plot(prep: Dict[str, Any],
                    projection: str = "albers",
                    coastlines: bool = False,
                    vmin: float = 0.0, vmax: float = 200.0,
                    figsize=(5, 5),
                    cbar_borderpad: float = -1.0,
                    cbar_label: Optional[str] = None,
                    cmap: Optional[Colormap or str] = None): 
    """
    Build the Cartopy figure/axes + pcolormesh + inset colorbar.
    Returns (fig, ax, pc, set_title) where set_title(text) updates the title.
    """


    Z      = prep["Z"]
    lon_e  = prep["lon_e"]
    lat_e  = prep["lat_e"]
    extent = prep["extent"]
    meta   = prep["meta"]

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
    cm = _resolve_cmap(cmap)
    fig, ax = plt.subplots(figsize=figsize,
                           subplot_kw={"projection": proj_map},
                           constrained_layout=True)
    ax.set_facecolor("#6a3e25")
    ax.set_anchor("W")
    ax.set_extent(list(extent), crs=proj_data)

    if coastlines:
        ax.add_feature(cfeature.COASTLINE.with_scale("10m"), linewidth=0.6)

    # Gridlines
    gl = ax.gridlines(crs=proj_data, draw_labels=True,
                      linewidth=0.5, alpha=0.4, linestyle=":",
                      x_inline=False, y_inline=False)
    gl.right_labels  = False
    gl.top_labels    = False
    gl.left_labels   = True
    gl.bottom_labels = True
    gl.rotate_labels = False
    gl.xformatter    = LONGITUDE_FORMATTER
    gl.yformatter    = LATITUDE_FORMATTER
    gl.xlocator      = mticker.MultipleLocator(1)

    lat_min, lat_max = float(lat_e.min()), float(lat_e.max())
    gl.ylocator      = mticker.FixedLocator(
        np.arange(np.floor(lat_min), np.ceil(lat_max) + 1, 1)
    )
    gl.xlabel_style  = {"size": 9}
    gl.ylabel_style  = {"size": 9}
    gl.xpadding      = 3
    gl.ypadding      = 3

    # pcolormesh (edges) + colorbar inset
    pc = ax.pcolormesh(lon_e, lat_e, Z,
                       transform=proj_data, shading="flat",
                       cmap=cm, vmin=vmin, vmax=vmax)

    cax = inset_axes(ax, width="3%", height="85%",
                     loc="center right", borderpad=cbar_borderpad)
    cb  = fig.colorbar(pc, cax=cax, orientation="vertical")
    cb.set_label(cbar_label or _var_label(meta))
    # Title setter
    def set_title(text: str):
        ax.set_title(text, pad=12)

    return fig, ax, pc, set_title


# ============== still image wrapper (uses shared factory) ==============
def plot_ice_map(prep: Dict[str, Any],
                 output_path: Optional[str] = None,
                 projection: str = "albers",
                 coastlines: bool = False,
                 vmin: float = 0.0, vmax: float = 200.0,
                 title: Optional[str] = None,
                 figsize=(5, 5),
                 cbar_borderpad: float = -1.0,
                 cmap: Optional[Colormap or str] = None,
                 dpi: Optional[int] = 300,
                 show: bool = True):
    """
    Plot prepared ice map. Saves to output_path if given.
    """
    meta = prep["meta"]
    if title is None:
        units = f" ({meta['units']})" if meta["units"] else ""
        title = meta["var_name"].replace("_", " ").title() + units

    fig, ax, pc, set_title = create_ice_plot(
        prep,
        projection=projection,
        coastlines=coastlines,
        vmin=vmin, vmax=vmax,
        figsize=figsize,
        cbar_borderpad=cbar_borderpad,
        cmap = cmap
    )
    set_title(title)

    if output_path:
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        plt.savefig(output_path, bbox_inches="tight", dpi=dpi or 300)

    if show:
        plt.show()
    else:
        plt.close(fig)

    return fig, ax


# ============== animation helpers ==============
def _parse_set_years_from_name(fname: str):
    """Return (y1, y2, set_code) for names like ice_season_2007-2008_set_A002.nc"""
    m = re.search(
        r"ice_season_(\d{4})-(\d{4})_set_([A-Z]\d{3})\.nc$",
        os.path.basename(fname)
    )
    if not m:
        return None
    return int(m.group(1)), int(m.group(2)), m.group(3)


def list_set_files(input_dir: str, set_code: str):
    """Find and sort all files for a given set code by start year."""
    pattern = os.path.join(input_dir, f"ice_season_*_set_{set_code}.nc")
    files = [f for f in glob.glob(pattern) if _parse_set_years_from_name(f)]
    files.sort(key=lambda f: _parse_set_years_from_name(f)[0])
    return files


# ============== animation wrapper (uses shared factory) ==============
def make_ice_animation(input_dir: str,
                       bathy_file: str,
                       set_code: str = "A002",
                       var_name: str = "ice_season_length",
                       projection: str = "albers",
                       coastlines: bool = False,
                       vmin: float = 0.0, vmax: float = 200.0,
                       figsize=(5, 5),
                       cbar_borderpad: float = -1.0,
                       fps: int = 4,
                       out_dir: Optional[str] = None,
                       out_basename: Optional[str] = None,
                       writer: Optional[str] = None,
                       cmap: Optional[Colormap or str] = None):
    """
    Build an animation over all 'ice_season_YYYY-YYYY_set_<set_code>.nc' files.
    Saves MP4 if ffmpeg available, else GIF.
    """
    files = list_set_files(input_dir, set_code)
    if not files:
        raise FileNotFoundError(f"No files found for set {set_code} in {input_dir}")

    # Prepare first for geometry/projection
    prep0 = prepare_ice_field(files[0], bathy_file, var_name=var_name)

    # Create shared figure/axes/pc via the SAME factory
    fig, ax, pc, set_title = create_ice_plot(
        prep0,
        projection=projection,
        coastlines=coastlines,
        vmin=vmin, vmax=vmax,
        figsize=figsize,
        cbar_borderpad=cbar_borderpad,
        cmap = cmap
    )

    title_base = f"{var_name.replace('_',' ')} — set {set_code}"
    set_title(title_base)

    # Preload fields + labels
    Z_list, labels = [], []
    for f in files:
        prep = prepare_ice_field(f, bathy_file, var_name=var_name)
        if prep["Z"].shape != prep0["Z"].shape:
            raise ValueError(f"Grid shape differs in {f}: {prep['Z'].shape} vs {prep0['Z'].shape}")
        Z_list.append(prep["Z"])
        y1, y2, _ = _parse_set_years_from_name(f)
        labels.append(f"{y1}–{y2}")

    # Update callback
    def update(i):
        Zi = Z_list[i]
        pc.set_array(np.ma.masked_invalid(Zi).ravel())
        set_title(f"{title_base} — {labels[i]}")
        return (pc,)

    anim = FuncAnimation(fig, update, frames=len(Z_list), interval=1000 // fps,
                         blit=False, repeat=True)

    # Pick writer
    if writer is None:
        try:
            _ = FFMpegWriter(fps=fps)
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


# ======== MASK & REDUCTION UTILITIES ========

def load_bathy_depth(bathy_file: str, like_da: xr.DataArray) -> xr.DataArray:
    """
    Load a depth field from bathy_file and align it to like_da (same dims/order).
    Tries common variable names if you don't know the exact one.
    Returns an xarray.DataArray aligned like like_da (no copy to np unless needed).
    """
    dsb = xr.open_dataset(bathy_file)

    # Try common candidates
    candidates = [
        "bathy_meter", "Bathymetry", "bathymetry", "depth", "DEPTH", "BATHY",
        "Bathy", "Bathy_Depth"
    ]
    depth_var = None
    for name in candidates:
        if name in dsb.data_vars:
            depth_var = name
            break
    if depth_var is None:
        # Fall back: pick the first non-mask variable
        non_mask = [k for k in dsb.data_vars if "msk" not in k.lower()]
        if not non_mask:
            raise KeyError("No bathymetry depth variable found in bathy file.")
        depth_var = non_mask[0]

    depth = xr.DataArray(np.array(dsb[depth_var]), dims=like_da.dims).reindex_like(like_da)
    return depth


def make_mask_bbox(prep: Dict[str, Any],
                   lon_min: float, lon_max: float,
                   lat_min: float, lat_max: float) -> np.ndarray:
    """
    Boolean mask (True=selected) for a lon/lat bounding box using center coords.
    """
    lon_1d = prep["lon_1d"]; lat_1d = prep["lat_1d"]
    # Broadcast to 2D grid
    LON, LAT = np.meshgrid(lon_1d, lat_1d)
    m = (LON >= lon_min) & (LON <= lon_max) & (LAT >= lat_min) & (LAT <= lat_max)
    # Also exclude NaNs in data
    return m & np.isfinite(prep["Z"])


def make_mask_bathy(prep: Dict[str, Any],
                    bathy_file: str,
                    threshold_m: float,
                    mode: str = "shallower") -> np.ndarray:
    """
    Boolean mask from bathymetry: 'shallower' selects depth <= threshold,
    'deeper' selects depth >= threshold.
    """
    # Build a dummy like_da for alignment
    like_da = xr.DataArray(prep["Z"], dims=("y", "x"))
    depth = load_bathy_depth(bathy_file, like_da)  # xarray.DataArray aligned to Z
    D = np.array(depth)

    if mode.lower() in ("shallower", "shallow", "<=", "le"):
        m = (D <= threshold_m)
    elif mode.lower() in ("deeper", "deep", ">=", "ge"):
        m = (D >= threshold_m)
    else:
        raise ValueError("mode must be 'shallower' or 'deeper'")

    return m & np.isfinite(prep["Z"])


def combine_masks(*masks: np.ndarray, how: str = "and") -> np.ndarray:
    """
    Combine multiple boolean masks (same shape) via 'and' or 'or'.
    """
    if not masks:
        raise ValueError("No masks provided")
    out = masks[0].astype(bool).copy()
    for m in masks[1:]:
        if m.shape != out.shape:
            raise ValueError("All masks must have the same shape")
        if how.lower() == "and":
            out &= m
        elif how.lower() == "or":
            out |= m
        else:
            raise ValueError("how must be 'and' or 'or'")
    return out


def reduce_over_mask(Z: np.ndarray,
                     mask: np.ndarray,
                     stat: str = "mean",
                     q: float = 0.5,
                     min_valid: int = 5) -> float:
    """
    Reduce Z over the True cells in mask, ignoring NaNs.
      stat: 'mean' | 'median' | 'percentile'
      q   : percentile in [0,1] if stat='percentile'
      min_valid: require at least this many finite samples, else return np.nan
    """
    values = Z[mask]
    values = values[np.isfinite(values)]
    if values.size < min_valid:
        return np.nan
    if stat == "mean":
        return float(np.nanmean(values))
    elif stat == "median":
        return float(np.nanmedian(values))
    elif stat == "percentile":
        return float(np.nanpercentile(values, 100.0 * q))
    else:
        raise ValueError("Unknown stat. Use 'mean', 'median', or 'percentile'.")
