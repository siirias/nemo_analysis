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

do_average =False 
do_points = True
ss = smh()  # create an object, to use for calling
variable_name = "SOLRAD"  # 'TAIR', 'LWRAD_DN', 'SOLRAD'
file_start = "solrad"  # 'tair', 'lwdrad', 'solrad'
#variable_name = "TAIR"  # 'TAIR', 'LWRAD_DN', 'SOLRAD'
#file_start = "tair"  # 'tair', 'lwdrad', 'solrad'
specific_points =[\
{'lat':65.0,'lon':23.0,'name':'BayofBotnia', 'data':[]},
{'lat':61.5,'lon':19.5,'name':'BothnianSea', 'data':[]},
{'lat':65.67,'lon':24.52,'name':'KemiAjos', 'data':[]},
{'lat':65.39,'lon':24.1,'name':'KemiLighthouse', 'data':[]},
{'lat':65.67,'lon':25.7,'name':'KemiLAND', 'data':[]},
{'lat':65.6,'lon':24.3,'name':'KemiOFHSHORE', 'data':[]},
{'lat':63.44,'lon':21.07,'name':'Mustasaari', 'data':[]},
{'lat':61.14,'lon':21.52,'name':'Rauma', 'data':[]},
{'lat':60.30,'lon':19.13,'name':'Hammarland', 'data':[]},
]

lat_min = 60.
lat_max = 66.
lon_min = 16.
lon_max = 26.

#specific_points =[\
#{'lat':65.67,'lon':25.7,'name':'KemiLAND', 'data':[]},
#{'lat':65.6,'lon':24.3,'name':'KemiOFHSHORE', 'data':[]}
#]

#base_directory = '/arch/smartsea/analysis/forcings/'
base_directory = '/lustre/tmp/siirias/c/smartsea_forcing_smhi/'

data_sets = [\

{'name':'A001',\
'dir':base_directory+'RCA4_NEMO_run_471-windamp/'},
{'name':'A002',\
'dir':base_directory+'RCA4_NEMO_run_475-windamp/'},
{'name':'A005',\
'dir':base_directory+'RCA4_NEMO_run_473-windamp/'},
{'name':'B001',\
'dir':base_directory+'RCA4_NEMO_run_472-windamp/'},
{'name':'B002',\
'dir':base_directory+'RCA4_NEMO_run_476-windamp/'},
{'name':'B005',\
'dir':base_directory+'RCA4_NEMO_run_474-windamp/'},
{'name':'C001',\
'dir':base_directory+'RCA4_NEMO_run_502-gregorian-windamp/'},
{'name':'C002',\
'dir':base_directory+'RCA4_NEMO_run_506-gregorian-windamp/'},
{'name':'D001',\
'dir':base_directory+'RCA4_NEMO_run_501-gregorian-windamp/'},
{'name':'D002',\
'dir':base_directory+'RCA4_NEMO_run_505-gregorian-windamp/'},
{'name':'D005',\
'dir':base_directory+'RCA4_NEMO_run_503-gregorian-windamp/'},
{'name':'hindcast',\
'dir':base_directory+'RCA4_NEMO_run_477-windamp/'}

]
#data_sets=[data_sets[-3]] ##remove!
output_dir = \
"/arch/smartsea/analysis/derived_data/atm_forcing/new/"
for data_set in data_sets:
    total_mean_var = np.zeros(0)
    time_axis = []
    name_marker = data_set['name']
    input_dir = data_set['dir']
    files = os.listdir(input_dir)
    files = [x for x in files if re.search(file_start+'.*\.nc$', x)] 
    files.sort()
    if do_points:
        for point in specific_points:
            point['data']=[]  # clear old data
    #files = files[:5] # just a few for testing
    for f in files:   
        year = int(re.search('y(\d\d\d\d)',f).groups()[0])
        print("file: {}, set {}".format(f,name_marker))
        data = Dataset(input_dir+f)
        lats = data.variables['lat'][:,:]
        lons = data.variables['lon'][:,:]
        time = list(data.variables['time'][:])
        the_var = data.variables[variable_name][:,:,:]
        data.close()
        if  np.min(the_var.shape)>10:  # otherwise the file is broken, ignore.
            time = [datetime.datetime(year,1,1)+datetime.timedelta(hours=int(x)) for x in time]
            time_axis+=time
            #mean temperature
            if do_average:
                areas = ss.give_areas(lats,lons)    
                areas[lats<lat_min] = 0.0
                areas[lats>lat_max] = 0.0
                areas[lons<lon_min] = 0.0
                areas[lons>lon_max] = 0.0
                total_area = np.sum(areas)
                mean_var = np.sum(the_var * areas,(1,2))/total_area
                if total_mean_var.size == 0:
                    total_mean_var = mean_var
                else:
                   total_mean_var =  \
                        np.concatenate((total_mean_var,mean_var))
            # specific points:
            if do_points:
                for point in specific_points:
                    point_index = (np.abs(lons[:,:]-point['lon'])+\
                                   np.abs(lats[:,:]-point['lat'])).argmin()
                    (point_j, point_i) = np.unravel_index(point_index,lons.shape)
                    point['data']=np.concatenate((point['data'],the_var[:,point_j,point_i]))
    if do_average:
        full_data=pd.DataFrame({'time':time_axis,'avg_temp':total_mean_var})
        full_data.to_csv(output_dir+'forcing_mean_{}_{}.csv'.format(variable_name, name_marker),\
                index=False)
    if do_points:
        for point in specific_points:
            full_data=pd.DataFrame({'time':time_axis,'avg_temp':point['data']})
            full_data.to_csv(output_dir+'{}_{}_{}.csv'\
                .format(point['name'],variable_name, name_marker),index=False)
    print(total_mean_var.shape, len(time_axis))
