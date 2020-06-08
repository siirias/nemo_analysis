import datetime
import numpy as np
from scipy.io import netcdf
from smartseahelper import smh
import os
import warnings
import multiprocessing as mpr
import time
import random 
import xarray as xr
import re
ss = smh()
ss.grid_type = 'T'
ss.interval = 'd'
ss.save_interval = 'year'
#ss.root_data_in = "/lustre/tmp/siirias/o/tmp/"  # gludge as the main disk is not sure enough.
#ss.root_data_in = "/arch/smartsea/analysis/"
ss.root_data_in = "/scratch/project_2001635/siiriasi/smartsea_data/"
ss.root_data_out = "/scratch/project_2001635/siiriasi/smartsea_data/"
startdate = datetime.datetime(1975, 1, 1)
enddate = datetime.datetime(2060, 12, 31)
folder_start = ''
# folder_start = 'OUTPUT'
#name_markers = ['new_REANALYSIS']
#name_markers = ['D001','C001','D002', 'D005', 'C002']
#name_markers = ['D001']
name_markers = ['A001']
variable_cover = 'soicecov'
ice_limit = 0.2
def empty_field_like(sample_data, new_name = 'noname'):
    empty = sample_data[0,:,:] # assumes time, x, y
    empty[:,:] = np.zeros(empty.shape)
    return empty.rename(new_name)

for name_marker in name_markers:
    ss.file_name_format = 'NORDIC-GOB_1{}_{}_{}_grid_{}.nc'
    if 'D' in name_marker or 'C' in name_marker:
        ss.file_name_format = 'SS-GOB_1{}_{}_{}_grid_{}.nc'
    if '1' in name_marker: # the 001 series are hindcasts, all other scenarios
        startdate = datetime.datetime(1975, 1, 1)
        enddate = datetime.datetime(2005, 12, 31)
    elif 'REANALYSIS' in name_marker:
        startdate = datetime.datetime(1980, 1, 1)
        enddate = datetime.datetime(2012, 12, 31)
        ss.save_interval = 'year'
        ss.file_name_format = 'NORDIC-GoB_1{}_{}_{}_grid_{}.nc'
    ss.main_data_folder= ss.root_data_in+"/{}{}/".format(folder_start, name_marker)
    filenames = ss.filenames_between(startdate, enddate)
    ok_files = 0
    files_working = []
    for f in filenames:
        if(os.path.isfile(ss.main_data_folder+f)):
            ok_files += 1
            files_working.append(f)
        else:
            print(f)
    print()
    print("ok {} out of {}".format(ok_files, len(filenames)))
    print(ss.main_data_folder, f)

    first_year = True
    for f in files_working:
        print(f)
        this_years_data = xr.open_dataset(ss.root_data_in+name_marker+'/'+f)[variable_cover]
        first_water = empty_field_like(this_years_data) 
        last_water = empty_field_like(this_years_data) 

        if(first_year):
            first_year = False
            out_name = re.search('(.*)\.nc',f).groups()[0]+'_icetest.nc'
            prev_data = this_years_data[:,:,:]
        else: 
            data = xr.concat([prev_data, this_years_data], 'time_counter')
            #first find real first and last ice (which is used to find the mid-winter)
            first_real_ice = empty_field_like(data) 
            last_real_ice = empty_field_like(data) 
            for time in range(data.shape[0]):
                free_water = empty_field_like(data,'first_real_ice_day')
                free_water = free_water.where(data[time,:,:] < ice_limit, time)
                first_real_ice = free_water.where(first_water < 0.001 , first_real_ice)
            for time in range(data.shape[0]-1,-1,-1):
                free_water = empty_field_like(data,'last_real_ice_day')
                free_water = free_water.where(data[time,:,:] > ice_limit, time) # add last year's days here
                last_real_ice = free_water.where(last_water < 0.001 , last_real_ice)
            
            for time in range(data.shape[0]-1,-1,-1):
                free_water = empty_field_like(data,'first_ice_day')
                free_water = free_water.where(data[time,:,:] > ice_limit, time)
                first_water = free_water.where(first_water < 0.001 , first_water)
            for time in range(data.shape[0]):
                free_water = empty_field_like(data,'last_ice_day')
                free_water = free_water.where(data[time,:,:] > ice_limit, time+prev_data_shape[0]) # add last year's days here
                last_water = free_water.where(last_water < 0.001 , last_water)
            first_ice = first_water + 1
            first_ice = data[0,:,:].where(np.isnan(data[0,:,:]),first_ice)
            last_ice = last_water - 1
            last_ice = data[0,:,:].where(np.isnan(data[0,:,:]),last_ice)
            out_dir = ss.root_data_out + '/tmp/'
            tmp = xr.merge([prev_first_water, last_water], compat='override')
            tmp.to_netcdf(out_dir + out_name)
            out_name = re.search('(.*)\.nc',f).groups()[0]+'_icetest.nc'
        prev_first_water = first_water[:,:]
        prev_data_shape = this_years_data.shape
