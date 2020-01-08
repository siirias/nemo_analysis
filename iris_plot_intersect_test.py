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


data_dir = "/arch/smartsea/analysis/derived_data/"
out_dir = "/arch/smartsea/analysis/derived_data/series_salt/"
files = [
         "average_salinity_new_REANALYSIS_depthlat.nc",
         "average_salinity_REANALYSIS_SMHI_depthlat.nc",
         "average_salinity_REANALYSIS_depthlat.nc",
         "average_salinity_REANALYSIS_ERA40_new_ice_depthlat.nc"
        ]
ifile = files[3]
#set_name = "heat_content"
set_name = "salinity"
data = iris.load(data_dir+ifile)[0]
#data_now = data.extract(iris.Constraint(time = \
#             lambda t: datemin < t < datemax))
contours = np.arange(0,23e17,1e17)
contours = np.arange(0.5,9.5,0.25)
for time_frame in range(data.shape[0]):
    data_now = data[time_frame,:,:]

    plt.clf()
    iqplt.contourf(data_now, contours, coords =['latitude','depth'])
    
    contour_set = iqplt.contour(data_now, [1.,3.,5.,7.,8.],\
                                 coords =['latitude','depth'],\
                                 colors = 'k',\
                                 linewidths=1.0,\
                                 alpha=0.3)
    plt.clabel(contour_set, fmt='%1.0f')
    #iqplt.contourf(data_now, coords =['latitude','depth'])
    #plt.gca().xaxis.axis_date()
    plt.title("{} ({})".format(set_name,data.coord('time').cell(time_frame)))
    plt.gca().invert_yaxis()
    plt.draw()
    plt.savefig(out_dir+"{}_{:05d}.png".format(set_name,time_frame))
    print("{} of {}".format(time_frame+1,data.shape[0]))
print("All done!")
