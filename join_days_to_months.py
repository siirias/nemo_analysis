#!/usr/bin/env ipython
# -*- coding: utf-8 -*-
"""

@author: siirias
"""

import datetime
import numpy as np
from scipy.io import netcdf
import netCDF4 as nc4
from netCDF4 import Dataset
from smartseahelper import smh
import os

input_files =[  "climatology_A001_d_icecon.nc",\
                "climatology_A002_d_icecon.nc",\
                "climatology_A005_d_icecon.nc",\
                 "climatology_A001_d_icevolume.nc",\
                "climatology_A002_d_icevolume.nc",\
                "climatology_A005_d_icevolume.nc"
            ]
input_dir = "/arch/smartsea/analysis/derived_data/test/"
output_dir = "/arch/smartsea/analysis/derived_data/test/"
for input_file in input_files:
    output_file = "monthly_"+input_file
    original_data = Dataset(input_dir + input_file)
    variables = original_data.variables.keys()
    # get the first variable with "mean" in it's name
    variable = [i for i in variables if 'mean' in i][0] 
    data = original_data.variables[variable][:,:,:]
    size_lat = data.shape[1]
    size_lon = data.shape[2]
    new_data = np.ma.masked_array(np.zeros((12,size_lat,size_lon)))
    month_lengths = [31,29,31,30,31,30,31,31,30,31,30,31]
    start_num = 0
    end_num = 0
    month_num = 1
    for m in month_lengths:
        end_num += m
        print(start_num,end_num)
        new_data[month_num-1,:,:] = np.mean(data[start_num:end_num,:,:],0)
        start_num = end_num+1
        month_num +=1
    new_data.mask[new_data==0.0]= True  # make sure, areas with no ice are masked.
    #create new netcdf
    new_file=nc4.Dataset(output_dir+output_file,'w',format='NETCDF4')
    new_file.createDimension('x',original_data.dimensions['x'].size)
    new_file.createDimension('y',original_data.dimensions['y'].size)
    new_file.createDimension('time',12)
    variable_axis=('time','y','x')
    clim_lat=new_file.createVariable('nav_lat','f4',('y','x'))
    clim_lon=new_file.createVariable('nav_lon','f4',('y','x'))
    clim_samples=new_file.createVariable('sample_number','i4','time')
    clim_day=new_file.createVariable('month_of_year','i4','time')
    clim_day[:]=np.arange(1,13)
    clim_samples[:]=np.zeros((12,))
    clim_variable_mean=new_file.createVariable(variable,'f4',variable_axis)
    clim_lat[:]=original_data.variables['nav_lat'][:]
    clim_lon[:]=original_data.variables['nav_lon'][:]
    clim_variable_mean[:,:,:] = new_data
    clim_lon.units = 'degrees east'
    clim_lat.units = 'degrees north'
    clim_day.units = 'month of year'
    clim_variable_mean.units = original_data.variables[variable].units
    clim_variable_mean.coordinates = "time_centered nav_lon nav_lat"                
    new_file.close()
    original_data.close()
