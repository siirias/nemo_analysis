#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import re
import glob
from typing import List, Tuple, Optional

import numpy as np
import matplotlib.pyplot as plt

# Import from your existing module (rename file if needed)
# from ice_viz import prepare_ice_field
import ice_helper as ih


# ---------------- file utilities (same naming logic you already use) ----------------

def parse_set_years(fname: str):
    """Return (y1, y2, set_code) from 'ice_season_2007-2008_set_A002.nc'."""
    m = re.search(r"ice_season_(\d{4})-(\d{4})_set_([A-Z]\d{3})\.nc$", os.path.basename(fname))
    if not m:
        return None
    return int(m.group(1)), int(m.group(2)), m.group(3)

def list_set_files(input_dir: str, set_code: str):
    """Find and sort all files for a given set by start year."""
    pattern = os.path.join(input_dir, f"ice_season_*_set_{set_code}.nc")
    files = [f for f in glob.glob(pattern) if parse_set_years(f)]
    files.sort(key=lambda f: parse_set_years(f)[0])
    return files


# ---------------- sampling helpers ----------------

def _nearest_index(vec: np.ndarray, value: float) -> int:
    """Index of nearest element in 1D array (assumes finite values)."""
    return int(np.nanargmin(np.abs(vec - value)))

def sample_point_from_prep(prep: dict, lon: float, lat: float) -> float:
    """
    Sample the prepared field at (lon, lat) using nearest cell.
    Returns np.nan if the nearest cell is NaN (e.g., land).
    """
    lon_1d = prep["lon_1d"]
    lat_1d = prep["lat_1d"]
    Z      = prep["Z"]

    ii = _nearest_index(lon_1d, lon)
    jj = _nearest_index(lat_1d, lat)
    val = Z[jj, ii]
    if not np.isfinite(val):
        val = np.nan
    return float(val)


# ---------------- core: extract series ----------------

def extract_point_series(input_dir: str,
                         bathy_file: str,
                         set_code: str,
                         var_name: str,
                         lon: float, lat: float):
    """
    For a given (lon, lat), return (years, values) over all seasons in set_code.
    'years' are start years (e.g., 2006 for 2006–2007).
    """
    files = list_set_files(input_dir, set_code)
    if not files:
        raise FileNotFoundError(f"No files for set {set_code} in {input_dir}")

    years = []
    values = []

    for f in files:
        y1, y2, _ = parse_set_years(f)
        prep = ih.prepare_ice_field(f, bathy_file, var_name=var_name)
        v = sample_point_from_prep(prep, lon, lat)
        years.append(y1)
        values.append(v)

    return np.array(years, dtype=int), np.array(values, dtype=float)


# ---------------- plotting ----------------

def plot_point_series(years: np.ndarray,
                      values: np.ndarray,
                      lon: float, lat: float,
                      set_code: str,
                      var_label: str = "Ice season length (days)",
                      out_path: Optional[str] = None,
                      ylim: Optional[Tuple[float, float]] = None,
                      figsize=(7, 3.2),
                      grid: bool = True):
    """
    Simple time-series plot (year vs value). Saves PNG if out_path provided.
    """
    fig, ax = plt.subplots(figsize=figsize, constrained_layout=True)
    ax.plot(years, values, marker="o", linewidth=1.2)

    ax.set_xlabel("Season start year")
    ax.set_ylabel(var_label)
    title = f"{var_label} at ({lon:.3f}°E, {lat:.3f}°N) — set {set_code}"
    ax.set_title(title, pad=8)

    if ylim is not None:
        ax.set_ylim(*ylim)

    if grid:
        ax.grid(True, linewidth=0.5, alpha=0.5)

    if out_path:
        os.makedirs(os.path.dirname(out_path), exist_ok=True)
        fig.savefig(out_path, dpi=200, bbox_inches="tight")

    return fig, ax


# ---------------- optional: CSV writer ----------------

def save_series_csv(years: np.ndarray, values: np.ndarray, out_csv: str):
    import csv
    os.makedirs(os.path.dirname(out_csv), exist_ok=True)
    with open(out_csv, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["season_start_year", "value"])
        for y, v in zip(years, values):
            w.writerow([int(y), "" if np.isnan(v) else float(v)])


# ---------------- main example ----------------

def main():
    bathy_file = r"c:/Data/NemoTest/bathy_meter.nc"
    input_dir  = r"C:\Data\VanhataloEtAl\ice_seasons_2006to2100/"
    set_code   = "A002"
    var_name   = "ice_season_length"

    # Choose a point (example: near 22.0E, 63.5N)
    lon, lat   = 22.0, 63.5

    years, values = extract_point_series(input_dir, bathy_file, set_code, var_name, lon, lat)

    # Build a nice y-label from metadata name + units (lightweight, manual here)
    var_label = "Ice season length (days)"  # if you want, you can import _var_label via a tiny helper

    out_dir  = os.path.join(input_dir, "timeseries")
    png_name = f"timeseries_{set_code}_{var_name}_{lon:.3f}E_{lat:.3f}N.png"
    csv_name = f"timeseries_{set_code}_{var_name}_{lon:.3f}E_{lat:.3f}N.csv"

    plot_point_series(years, values, lon, lat, set_code,
                      var_label=var_label,
                      out_path=os.path.join(out_dir, png_name),
                      ylim=(0, 250), figsize=(8, 3.5))

    save_series_csv(years, values, os.path.join(out_dir, csv_name))


if __name__ == "__main__":
    main()
