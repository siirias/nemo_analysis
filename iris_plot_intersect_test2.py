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
import datetime

comparison = True
#plot_type = "longitude_crop"  # Plot boundary area 
plot_type = "latitude_mean"   # plot dep,lat over ageraged area
data_dir = "/arch/smartsea/analysis/derived_data/"
out_dir = "/arch/smartsea/analysis/derived_data/series_salt/"
datemin = datetime.datetime(1980,1,1)
datemax = datetime.datetime(2008,12,31)
if plot_type in ['latitude_mean']:
    axis_to_plot = ['latitude','depth']
    set_name = "average_diff_T"
    files = [
             "average_salinity_new_REANALYSIS_depthlat.nc",
             "average_salinity_REANALYSIS_SMHI_depthlat.nc",
             "average_salinity_REANALYSIS_depthlat.nc",
             "average_salinity_REANALYSIS_ERA40_new_ice_depthlat.nc",
             "average_salinity_SSR_5_depthlat.nc",
             "average_salinity_SSR_6_depthlat.nc",
             "border_slice_potential_temperature_SSR_5_depthlat.nc",
             "border_slice_potential_temperature_SSR_6_depthlat.nc"
            ]
if plot_type in ['longitude_crop']:
    axis_to_plot = ['longitude','depth']
    set_name = "boundary_diff_T"
    files = [
             "border_slice_potential_temperature_new_REANALYSIS_depthlon.nc",
             "border_slice_potential_temperature_REANAyyLYSIS_SMHI_depthlon.nc",
             "border_slice_potential_temperature_REANALYSIS_depthlon.nc",
             "average_salinity_SSR_5_depthlat.nc",
             "average_salinity_SSR_6_depthlat.nc"
            ]
colormap = None
ifile = files[6]
comp_file = files[7]
#set_name = "heat_content"
data = iris.load(data_dir+ifile)[0]
if plot_type in ['longitude_crop']:
    data = data.extract(iris.Constraint(coord_values = \
             {'longitude':lambda l: 18.5 < l < 23.}))
data = data.extract(iris.Constraint(time = \
             lambda t: datemin < t < datemax))
contour_lines = [1.,3.,5.,7.,8.]
#contours = np.arange(0,23e17,1e17)
contours = np.arange(0.5,9.5,0.25)
variable_name = data.name()
if comparison:
    data2 = iris.load(data_dir+comp_file)[0]
    if plot_type in ['longitude_crop']:
        data2 = data2.extract(iris.Constraint(coord_values = \
                     {'longitude':lambda l: 18.5 < l < 23.}))
    data2 = data2.extract(iris.Constraint(time = \
             lambda t: datemin < t < datemax))
    data.data = data.data - data2.data
    contour_lines = [-3.,-2., -1.5,-1.0,-0.6,-0.4,0.,0.4, 0.6,1.0,1.5, 2.0,3.]
    contour_lines = list(np.arange(-3,3,0.1)) 
#    contour_lines = list(np.array(contour_lines)*1.0)  # quick scaling.
    #contours = np.arange(0,23e17,1e17)
    #contours = np.arange(-3.5,3.5,0.1)
    contours = np.arange(-2.,2.,0.1)
#    contours = np.arange(-0.2,0.2,0.005)
    #contours = np.arange(-1.0,1.0,0.05)
    colormap =  'PiYG' #  'RdBu_r'  # 'PiYG'
for time_frame in range(data.shape[0]):
    data_now = data[time_frame,:,:]

    plt.clf()
    iqplt.contourf(data_now, contours,\
                     coords =axis_to_plot,\
                     cmap = colormap)
    
    contour_set = iqplt.contour(data_now, contour_lines,\
                                 coords =axis_to_plot,\
                                 colors = 'k',\
                                 linewidths=1.0,\
                                 alpha=0.3)
    plt.clabel(contour_set, fmt='%1.1f')
    #iqplt.contourf(data_now, coords =['latitude','depth'])
    #plt.gca().xaxis.axis_date()
    plt.title("{} {} ({})".format(\
                    variable_name,\
                    set_name,\
                    data.coord('time').cell(time_frame)))
    plt.gca().invert_yaxis()
    plt.draw()
    plt.savefig(out_dir+"{}_{:05d}.png".format(set_name,time_frame))
    print("{} of {}".format(time_frame+1,data.shape[0]))
print("All done!")
