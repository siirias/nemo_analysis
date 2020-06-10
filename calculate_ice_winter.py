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
#name_markers = ['A001', 'A002','A005','B001','B002','B005','D001','D002','D005']
name_markers = ['A002','A005','B002','B005','D002','D005']
variable_cover = 'soicecov'
continuous_ice_limit = 0.2
ice_limit = 0.3
crop_days = int(365/2)  # just to crop this much from start and end of two year set, to get the winter time.
def empty_field_like(sample_data, new_name = 'noname'):
    empty = sample_data[0,:,:] # assumes time, x, y
    empty[:,:] = np.zeros(empty.shape)
    return empty.rename(new_name)

for name_marker in name_markers:
    startdate = datetime.datetime(2006, 1, 1)
    enddate = datetime.datetime(2060, 12, 31)
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
            year = re.search('_1d_(\d\d\d\d)',f).groups()[0]
            out_name = 'ice_season_{}-{}_set_{}.nc'.format(year, str(int(year)+1), name_marker)
            prev_data = this_years_data[:,:,:]
        else: 
            data = xr.concat([prev_data, this_years_data], 'time_counter')
            # First let's just find the first and last ice day for the winter,
            # without thinking the continuous ice-time.
            # This is then used to guess the mid-winter,
            # and that is used to calculate the contnuous ice-season.
            first_ice = empty_field_like(data) 
            last_ice = empty_field_like(data) 
            for time in range(crop_days, data.shape[0]-crop_days):
                tmp = empty_field_like(data)
                tmp = tmp.where(data[time,:,:] < ice_limit, time)
                first_ice = tmp.where(first_ice < 0.001 , first_ice)
            first_ice = first_ice.rename('first_ice_day')+1 #+1 to fix the index 0

            for time in range(data.shape[0]-crop_days, crop_days, -1):
                tmp = empty_field_like(data)
                tmp = tmp.where(data[time,:,:] < ice_limit, time) 
                last_ice = tmp.where(last_ice < 0.001 , last_ice)
            last_ice = last_ice.rename('last_ice_day')+1 #+1 to fix the index 0   
            
            mid_winter = (first_ice + last_ice)/2.0
            mid_winter = mid_winter.rename('middle_day_of_the_ice_season') 
            ice_season_length = last_ice - first_ice
            ice_season_length = ice_season_length.rename('ice_season_length') 
            # Now we have mid winter, and first and last ice-occurences.
            # Next we need to figure out the continuous ice-season, around the mid-winter
            first_continuous_ice = empty_field_like(data) 
            last_continuous_ice = empty_field_like(data) 
            for time in range(crop_days, data.shape[0]-crop_days):
                tmp = empty_field_like(data)
                tmp = tmp.where(data[time,:,:] > continuous_ice_limit, time)
                # tmp has time in the ones with water
                tmp = tmp.where(mid_winter<time,0)
                # tmp has time in ones with water, where time is past midwinter
                last_continuous_ice = tmp.where(last_continuous_ice < 0.001 , last_continuous_ice)
            last_continuous_ice = last_continuous_ice.rename('last_continuous_ice_day')+1 #+1 to fix the index 0

            for time in range(data.shape[0]-crop_days, crop_days, -1):
                tmp = empty_field_like(data)
                tmp = tmp.where(data[time,:,:] > continuous_ice_limit, time)
                # tmp has time in the ones with water
                tmp = tmp.where(mid_winter>time,0)
                # tmp has time in ones with water, where time before  midwinter
                first_continuous_ice = tmp.where(first_continuous_ice < 0.001 , first_continuous_ice)
            first_continuous_ice = first_continuous_ice.rename('first_continuous_ice_day')+1 #+1 to fix the index 0
            continuous_mid_winter = (first_continuous_ice + last_continuous_ice)/2.0
            continuous_mid_winter = continuous_mid_winter.rename('middle_day_of_the_continuous_ice_season') 
            continuous_ice_season_length = last_continuous_ice - first_continuous_ice
            continuous_ice_season_length = continuous_ice_season_length.rename('continuous_ice_season_length') 

            length_diff = (ice_season_length - continuous_ice_season_length).rename('length_diff')
            mid_diff = (mid_winter - continuous_mid_winter).rename('mid_diff')
            start_diff = (first_ice - first_continuous_ice).rename('start_diff')
            end_diff = (last_ice - last_continuous_ice).rename('end_diff')

            out_dir = ss.root_data_out + '/tmp/'
            tmp = xr.merge([first_ice, last_ice, ice_season_length, mid_winter,\
                            first_continuous_ice, last_continuous_ice, \
                            continuous_ice_season_length, continuous_mid_winter,
                            length_diff, mid_diff, start_diff, end_diff], compat='override')
            tmp.to_netcdf(out_dir + out_name)
            year = re.search('_1d_(\d\d\d\d)',f).groups()[0]
#            out_name = re.search('(.*)\.nc',f).groups()[0]+'_icetest.nc'
            out_name = 'ice_season_{}-{}_set_{}.nc'.format(year, str(int(year)+1), name_marker)
        prev_first_water = first_water[:,:]
        prev_data_shape = this_years_data.shape
