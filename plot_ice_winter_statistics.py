import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

# ---------------- Paths ----------------
bathy_file = r"c:/Data/NemoTest/bathy_meter.nc"
input_dir  = r"C:\Data\VanhataloEtAl\ice_seasons_2006to2100/"
nc_file    = input_dir + "ice_season_2007-2008_set_B002.nc"
# ---------------- Switches ----------------
draw_coastlines = False

# ---------------- Open data ----------------
ds    = xr.open_dataset(nc_file)
bathy = xr.open_dataset(bathy_file)

# Data & sea mask aligned to grid
Z = ds["ice_season_length"]                               # (y,x)
msk = xr.DataArray(np.array(bathy["bdy_msk"]), dims=Z.dims).reindex_like(Z)
sea = (msk != 0)
Zp = Z.where(sea)                                         # NaN on land

# Pull lon/lat and mark zeros as NaN
lon2d = ds["nav_lon"].values.copy()                       # (y,x)
lat2d = ds["nav_lat"].values.copy()
lon2d[lon2d == 0] = np.nan
lat2d[lat2d == 0] = np.nan

Y, X = Z.shape

# -------- Helpers --------
def fill_1d_nearest(a):
    a = a.astype(float)
    good = np.isfinite(a)
    if not good.any():
        return np.arange(a.size, dtype=float)
    idx = np.arange(a.size, dtype=float)
    return np.interp(idx, idx[good], a[good])

def centers_to_edges(c):
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
    eps = np.finfo(float).eps * max(1.0, np.nanmax(np.abs(edges)))
    for k in range(1, edges.size):
        if edges[k] <= edges[k-1]:
            edges[k] = edges[k-1] + eps
    return edges

# -------- Recover 1D lon(x) and lat(y) robustly --------
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

lon_1d = fill_1d_nearest(lon_1d)
lat_1d = fill_1d_nearest(lat_1d)

# Build edges for pcolormesh (length X+1, Y+1)
lon_e = centers_to_edges(lon_1d)
lat_e = centers_to_edges(lat_1d)

# ---------------- Plot (Cartopy) ----------------
central_lon = float(np.nanmean(lon_1d))
central_lat = float(np.nanmean(lat_1d))
proj_map  = ccrs.LambertConformal(central_longitude=central_lon,
                                  central_latitude=central_lat)
proj_data = ccrs.PlateCarree()  # data are lon/lat degrees

cmap = plt.cm.viridis.copy()
cmap.set_bad("#6a3e25")  # land/background as brown

# Use constrained layout; DO NOT call tight_layout()
fig, ax = plt.subplots(
    figsize=(8, 9),
    subplot_kw={"projection": proj_map},
    constrained_layout=True
)
ax.set_facecolor("#6a3e25")
ax.set_anchor("W")  # keep the map frame snug to the left

# pcolormesh with explicit edges (stable & fast)
pc = ax.pcolormesh(
    lon_e, lat_e, Zp.values,
    transform=proj_data,
    shading="flat",
    cmap=cmap, vmin=0, vmax=200
)

# Coastlines
if draw_coastlines:
    ax.add_feature(cfeature.COASTLINE.with_scale("10m"), linewidth=0.6)

# Extent to data bounds
ax.set_extent([float(lon_e.min()), float(lon_e.max()),
               float(lat_e.min()), float(lat_e.max())],
              crs=proj_data)

# Gridlines: ONLY bottom longitudes and left latitudes; no inline labels
gl = ax.gridlines(
    crs=proj_data, draw_labels=True,
    linewidth=0.5, alpha=0.4, linestyle=":",
    x_inline=False, y_inline=False
)
gl.right_labels   = False
gl.top_labels     = False
gl.left_labels    = True
gl.bottom_labels  = True
gl.xformatter     = LONGITUDE_FORMATTER
gl.yformatter     = LATITUDE_FORMATTER
gl.xlocator       = mticker.MultipleLocator(1)      # integer 1° steps
# gl.ylocator     = mticker.MultipleLocator(1)      # optional
gl.xlabel_style   = {"size": 9}
gl.ylabel_style   = {"size": 9}
gl.xpadding       = 3
gl.ypadding       = 3

ax.set_title(
    "Ice season length",
    pad=12   # push it slightly up
)

# Colorbar as an inset (doesn't consume figure margin)
cax = inset_axes(ax, width="3%", height="85%", loc="center right", borderpad=-2.0)
cb = fig.colorbar(pc, cax=cax, orientation="vertical")
cb.set_label("Ice season length (days)")

plt.show()
