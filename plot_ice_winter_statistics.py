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

import ice_helper as ih

def main():
    bathy_file = r"c:/Data/NemoTest/bathy_meter.nc"
    input_dir  = r"C:\Data\VanhataloEtAl\ice_seasons_2006to2100/"
    var_name   = "first_continuous_ice_day"
#    var_name   = "ice_season_length"

    # one-off static example (kept, if you want)
    # nc_file = os.path.join(input_dir, "ice_season_2007-2008_set_A002.nc")
    # prep = prepare_ice_field(nc_file, bathy_file, var_name=var_name)
    # plot_ice_map(prep, output_path=None, projection="albers", coastlines=False,
    #              vmin=0.0, vmax=200.0, title="Ice season length", figsize=(5,5), show=True)
    for plot_set in ['A002','A005', 'B002', 'B005', 'D002', 'D005']:
        # animation for a chosen set:
        ih.make_ice_animation(
            input_dir=input_dir,
            bathy_file=bathy_file,
            set_code= plot_set, 
            var_name=var_name,
            projection="albers",
            coastlines=False,
#            cmap = cmo.cm.curl,
            cmap = 'hsv',
            vmin=150.0, vmax=360.0,
            figsize=(5, 5),
            fps=4,                      # 4 frames per second
            out_dir=os.path.join(input_dir, "animations"),
            out_basename=None,          # auto name
            writer=None                 # auto choose ffmpeg→mp4, else Pillow→gif
        )




if __name__ == "__main__":
    main()
