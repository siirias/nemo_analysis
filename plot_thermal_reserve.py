#!/usr/bin/env ipython
# -*- coding: utf-8 -*-
import datetime
import matplotlib as mp
import matplotlib.pyplot as plt
import numpy as np
from scipy.io import netcdf
from scipy.signal import savgol_filter
from mpl_toolkits.basemap import Basemap
from mpl_toolkits.basemap import Basemap, shiftgrid, cm
import netCDF4 as nc4
from netCDF4 import Dataset
from smartseahelper import smh
import os
import cmocean
import iris
import iris.plot as iplt
import iris.quickplot as qplt
import cf_units
input_dir="~/smartsea/derived_data/"
output_dir="/home/siirias/smartsea/derived_data/"
input_file="reserve_votemper_REANALYSIS.nc"
other_input_dir="~/smartsea/derived_data/"
other_input_file="reserve_votemper_new_REANALYSIS.nc"
cube_T=iris.load_cube(input_dir+input_file)
other_cube_T=iris.load_cube(other_input_dir+other_input_file)
for i in cube_T.coords('depth')[0].cells():
  print(i)

thermal_sum=cube_T.collapsed('depth',iris.analysis.SUM)
other_thermal_sum=other_cube_T.collapsed('depth',iris.analysis.SUM)
t1=iris.time.PartialDateTime(year=1980, month=1)
t2=iris.time.PartialDateTime(year=2010,month=12,day=31)
thermal_sum=thermal_sum.extract(iris.Constraint(time=lambda x: t1<x<t2))
other_thermal_sum=other_thermal_sum.extract(iris.Constraint(time=lambda x: t1<x<t2))
fig=plt.figure(figsize=(10,5))
plt.subplot(211)
plt.title('Heat content (non-unit) Old-New')
iplt.plot(thermal_sum)
iplt.plot(other_thermal_sum)
plt.subplot(212)
difference=thermal_sum-other_thermal_sum
smooth_diff=difference.copy()
smooth_diff.data=savgol_filter(difference.data,49,1)
iplt.plot(difference)
iplt.plot(smooth_diff)
iplt.plot(smooth_diff*0.0)
plt.savefig(output_dir+'heat_content_comparison.png')
plt.show()

#iplt.plot(cube_T[:,1]);plt.show()
