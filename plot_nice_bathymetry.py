#!/usr/bin/env python3
"""
Simple, clean bathymetry plot for /mnt/data/bathy_meter.nc using Cartopy.

- Uses known names: lon, lat, Bathymetry (m), with land == 0
- Equal-area projection (LAEA) centered on the dataset to minimize distortion
- Defaults to 0–155 m (98th percentile) for better contrast; --full-range shows full 0–281 m
- Adds labeled grid, coastline, readable colorbar

Usage:
    python plot_bathy_simple.py /mnt/data/bathy_meter.nc -o bathy.png
    # or to show full dynamic range:
    python plot_bathy_simple.py /mnt/data/bathy_meter.nc -o bathy_full.png --full-range
"""

import argparse
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature

# Fixed domain/projection from the file (inspected once so we don't need to guess)
LAEA_LON0 = 21.11085625
LAEA_LAT0 = 62.9164688
EXTENT = (16.40257, 25.8191425, 59.92485, 65.9080876)  # lon_min, lon_max, lat_min, lat_max
V_MAX_98PCT = 155.0  # good-contrast upper limit (rounded from ~155.39)
V_MAX_FULL = 281.0   # near the true max (~280.8)

def main(path, out_png, use_full, dpi, figsize):
    ds = xr.open_dataset(path)
    lon = ds["lon"].values
    lat = ds["lat"].values
    Z = ds["Bathymetry"].values  # meters, sea > 0, land == 0

    # Mask land
    Zm = np.where(Z > 0, Z, np.nan)

    # Build coordinates
    if lon.ndim == 1 and lat.ndim == 1:
        LON, LAT = np.meshgrid(lon, lat)
    else:
        LON, LAT = lon, lat

    # Colormap (cmocean if available)
    try:
        import cmocean
        cmap = cmocean.cm.deep
    except Exception:
        cmap = plt.get_cmap("viridis")

    # Depth range
    vmin = 0.0
    vmax = V_MAX_FULL if use_full else V_MAX_98PCT

    # Plot
    fig = plt.figure(figsize=figsize, dpi=dpi)
    proj = ccrs.LambertAzimuthalEqualArea(central_longitude=LAEA_LON0, central_latitude=LAEA_LAT0)
    ax = plt.axes(projection=proj)

    # Tight extent with a tiny pad
    pad_lon, pad_lat = 0.2, 0.2
    lon_min, lon_max, lat_min, lat_max = EXTENT
    ax.set_extent([lon_min - pad_lon, lon_max + pad_lon, lat_min - pad_lat, lat_max + pad_lat],
                  crs=ccrs.PlateCarree())

    # Context
    ax.add_feature(cfeature.LAND.with_scale("10m"), facecolor="#f3f3f3", edgecolor="none", zorder=1)
    ax.add_feature(cfeature.COASTLINE.with_scale("10m"), linewidth=0.5, zorder=3)

    # Raster
    im = ax.pcolormesh(LON, LAT, Zm,
                       transform=ccrs.PlateCarree(),
                       cmap=cmap, vmin=vmin, vmax=vmax,
                       shading="auto", zorder=2)

    # Optional crisp isobaths for orientation
    for iso in ( 15, 50, 100, 200):
        cs = ax.contour(LON, LAT, Zm, levels=[iso], transform=ccrs.PlateCarree(),
                        linewidths=0.4, linestyles="-", colors="k", alpha=0.35, zorder=4)
        ax.clabel(cs, fmt=lambda v: f"{int(v)} m", fontsize=7, inline=1)

    # Colorbar
    cb = plt.colorbar(im, ax=ax, orientation="vertical", fraction=0.030, pad=0.04)
    cb.set_label("Depth (m)")

    # Gridlines with labels
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                      linewidth=0.6, linestyle=":", alpha=0.75)
    gl.right_labels = False
    gl.top_labels = False

    # Title
    rng = f"0–{int(vmax)} m"
#    ax.set_title(f"Bathymetry ({rng}) • LAEA @ {LAEA_LAT0:.2f}°, {LAEA_LON0:.2f}°", fontsize=11)
    ax.set_title(f"Gulf of Bothnia Bathymetry, 1 NM model setup", fontsize=11)

    plt.tight_layout()
    plt.savefig(out_png, bbox_inches="tight")
    print(f"Saved: {out_png}")

if __name__ == "__main__":
    p = argparse.ArgumentParser(description="Plot bathymetry (simplified, equal-area).")
    p.add_argument("-n","--netcdf", default = "c:/Data/NemoTest/bathy_meter.nc", help="Path to bathy_meter.nc")
    p.add_argument("-o", "--out", dest="out_png", default="bathymetry_map.png", help="Output PNG filename")
    p.add_argument("--full-range", action="store_true", help="Use full 0–281 m range (less contrast)")
    p.add_argument("--dpi", type=int, default=200, help="Figure DPI")
    p.add_argument("--figsize", type=float, nargs=2, default=(10, 8), help="Figure size, e.g. --figsize 10 8")
    args = p.parse_args()
    main(args.netcdf, args.out_png, args.full_range, args.dpi, tuple(args.figsize))
