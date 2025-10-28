# !/usr/bin/env ipython
# -*- coding: utf-8 -*-
"""
@author: siirias
"""

import datetime
import matplotlib as mp
import matplotlib.pyplot as plt
import numpy as np
from scipy.io import netcdf
from smartseahelper import smh
import os
import cmocean
import pandas as pd
import iris
import iris.plot as iplt
import iris.quickplot as iqplt
import siri_omen
import siri_omen.utility as sou
import siri_omen.nemo_reader as nrd
import cf_units
import iris.util
import gsw  # TEOS-10
import re

in_file = "/arch/smartsea/analysis/experimental_run/long_test/SS-GOB_1d_19850101_19851231_grid_T.nc"
#in_file = "/lustre/tmp/siirias/c/nemo_gob_new_forcing/OUTPUT/SS-GOB_1h_19800101_19801231_grid_T.nc"
out_dir = "/arch/smartsea/analysis/derived_data/test/tmp_images/SSS/"
variable = "SSS"
#variable = "sea_surface_height"
vmin =  None
vmax =  None
#vmin = -0.2  # None
#vmax = 0.2   # None
print('1')
data = iris.load(in_file,variable)[0]
print('2')
data = nrd.remove_null_indices(data,fill_value=0.0)
print('3')
nrd.fix_cube_coordinates(data)
print('hep!')
time_coord = data.coord_dims('time')[0]

for step in range(data.shape[time_coord]):
    current_slice = [slice(None)]*len(data.shape)
    current_slice[time_coord] = step
    current_slice = tuple(current_slice)
    current = data[current_slice]
    this_time = re.search("\[.*\]",str(data.coord('time')[step])).group()
    plt.close()
    iqplt.contourf(current, vmin=vmin, vmax=vmax)
    plt.title("{} {}".format(variable,this_time))
    plt.draw()
    plt.savefig(out_dir+"{}_{}.png".format(variable,step))
    print(step)
plt.show()
