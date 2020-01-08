#!/usr/bin/env ipython
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  3 14:11:45 2018

@author: siirias
"""

import datetime
import numpy as np
from scipy.io import netcdf
import netCDF4 as nc4
from netCDF4 import Dataset
from smartseahelper import smh
import os
import re
import pandas as pd

ss = smh()  # create an object, to use for calling

data_sets = [\

{'name':'A001',\
'dir':'/lustre/tmp/siirias/c/smartsea_forcing_smhi/E-Hype_CSC_REMO2009_MPI-ESM-LR_rcp45/'},
{'name':'A002',\
'dir':'/lustre/tmp/siirias/c/smartsea_forcing_smhi/E-Hype_CSC_REMO2009_MPI-ESM-LR_rcp45/'},
{'name':'A005',\
'dir':'/lustre/tmp/siirias/c/smartsea_forcing_smhi/E-Hype_CSC_REMO2009_MPI-ESM-LR_rcp85/'},
{'name':'B001',\
'dir':'/lustre/tmp/siirias/c/smartsea_forcing_smhi/E-Hype_SMHI_RCA4_EC-EARTH_rcp45/'},
{'name':'B002',\
'dir':'/lustre/tmp/siirias/c/smartsea_forcing_smhi/E-Hype_SMHI_RCA4_EC-EARTH_rcp45/'},
{'name':'B005',\
'dir':'/lustre/tmp/siirias/c/smartsea_forcing_smhi/E-Hype_SMHI_RCA4_EC-EARTH_rcp85/'},
{'name':'C001',\
'dir':'/lustre/tmp/siirias/c/smartsea_forcing_smhi/E-Hype_IPSL-IPSL-CM5A-MR_rcp45/'},
{'name':'C002',\
'dir':'/lustre/tmp/siirias/c/smartsea_forcing_smhi/E-Hype_IPSL-IPSL-CM5A-MR_rcp45/'},
{'name':'D001',\
'dir':'/lustre/tmp/siirias/c/smartsea_forcing_smhi/E-Hype_SMHI_RCA4_HadGEM2-ES_rcp45/'},
{'name':'D002',\
'dir':'/lustre/tmp/siirias/c/smartsea_forcing_smhi/E-Hype_SMHI_RCA4_HadGEM2-ES_rcp45/'},
{'name':'D005',\
'dir':'/lustre/tmp/siirias/c/smartsea_forcing_smhi/E-Hype_SMHI_RCA4_HadGEM2-ES_rcp85/'},
{'name':'hindcast',\
'dir':'/lustre/tmp/siirias/c/smartsea_forcing_smhi/rco_riverdat_1850-2009_observed.I/'}

]
output_dir = \
"/arch/smartsea/analysis/derived_data/inflow/"
for data_set in data_sets:
    total_water_in = np.zeros(0)
    time_axis = []
    name_marker = data_set['name']
    input_dir = data_set['dir']
    files = os.listdir(input_dir)
    files = [x for x in files if re.search('runoff_.*\.nc', x)] 
    files.sort()
#    files = files[:20] # just a few for testing
    for f in files:   
        year = int(re.search('y(\d\d\d\d)',f).groups()[0])
        print("file: {}, set {}".format(f,name_marker))
        data = Dataset(input_dir+f)
        lats = data.variables['nav_lat'][:,:]
        lons = data.variables['nav_lon'][:,:]
        time = list(data.variables['time'][:])
        time_multiplier = 1  # gludge needed, as hindcast is monthly, not daily.
        if name_marker == 'hindcast':
            time_multiplier = 366.0/12.0
        runoff = data.variables['sorunoff'][:,:,:]
        data.close()
        areas = ss.give_areas(lats,lons)    
        water_in = runoff * areas

        time = [datetime.datetime(year,1,1)+\
                datetime.timedelta(int(float(x)*time_multiplier)) for x in time]
        time_axis+=time

        if total_water_in.size == 0:
            total_water_in = np.sum(water_in,(1,2))
        else:
           total_water_in =  np.concatenate((total_water_in,np.sum(water_in,(1,2))))
    full_data=pd.DataFrame({'time':time_axis,'inflow':total_water_in})
    full_data.to_csv(output_dir+'river_water_in_{}.csv'.format(name_marker),index=False)
    print(total_water_in.shape, len(time_axis))
