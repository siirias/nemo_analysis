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

do_average = True 
do_depths = True
ss = smh()  # create an object, to use for calling
specific_depths =[\
{'depth':5.0,'name':'5meter', 'data':[]},
{'depth':20.0,'name':'20meter', 'data':[]},
{'depth':80.0,'name':'80meter', 'data':[]},
{'depth':120.0,'name':'120meter', 'data':[]}
]
depths = \
np.array([1.501665, 4.516333, 7.549625, 10.60677, 13.69444, 16.82112, \
19.99763, 23.23766, 26.55851, 29.98192, 33.53508, 37.25171, 41.17324, \
45.34998, 49.84202, 54.71985, 60.06403, 65.96383, 72.51433, 79.81163, \
87.94651, 96.99683, 107.0199, 118.0461, 130.0755, 143.0776, 156.9952, \
171.7504, 187.2523, 203.405, 220.1133, 237.2879, 254.8483, 272.724, \
290.8549, 309.1907])
layer_thicknessess = depths.copy()
layer_thicknessess[1:]=np.diff(depths)
layer_thicknessess[0] = layer_thicknessess[1]
data_sets = [\

{'name':'A001',\
'dir':'/arch/smartsea/analysis/forcings/boundary/boundary_GoB_1nm_BalticApp_A001/'},
{'name':'A002',\
'dir':'/arch/smartsea/analysis/forcings/boundary/boundary_GoB_1nm_BalticApp_A002/'},
{'name':'A005',\
'dir':'/arch/smartsea/analysis/forcings/boundary/boundary_GoB_1nm_BalticApp_A005/'},
{'name':'B001',\
'dir':'/arch/smartsea/analysis/forcings/boundary/boundary_GoB_1nm_BalticApp_B001/'},
{'name':'B002',\
'dir':'/arch/smartsea/analysis/forcings/boundary/boundary_GoB_1nm_BalticApp_B002/'},
{'name':'B005',\
'dir':'/arch/smartsea/analysis/forcings/boundary/boundary_GoB_1nm_BalticApp_B005/'},
{'name':'C001',\
'dir':'/arch/smartsea/analysis/forcings/boundary/boundary_GoB_1nm_BalticApp_C001/'},
{'name':'C002',\
'dir':'/arch/smartsea/analysis/forcings/boundary/boundary_GoB_1nm_BalticApp_C002/'},
{'name':'D001',\
'dir':'/arch/smartsea/analysis/forcings/boundary/boundary_GoB_1nm_BalticApp_D001/'},
{'name':'D002',\
'dir':'/arch/smartsea/analysis/forcings/boundary/boundary_GoB_1nm_BalticApp_D002/'},
{'name':'D005',\
'dir':'/arch/smartsea/analysis/forcings/boundary/boundary_GoB_1nm_BalticApp_D005/'},
{'name':'hindcast',\
'dir':'/arch/smartsea/analysis/forcings/boundary/boundary_GoB_1nm_BalticApp_hindcast0011/'}

]
#data_sets=[data_sets[-3]] ##remove!
output_dir = \
"/arch/smartsea/analysis/derived_data/boundary/"

variable = 'vosaline'

for data_set in data_sets:
    if do_depths:
        for point in specific_depths:
            point['data']=[]
    total_mean_variable = np.zeros(0)
    time_axis = []
    name_marker = data_set['name']
    input_dir = data_set['dir']
    files = os.listdir(input_dir)
    files = [x for x in files if re.search('bdy_phys.*\.nc', x)] 
    files.sort()
    #files = files[:5] # just a few for testing
    for f in files:   
        year = int(re.search('y(\d\d\d\d)',f).groups()[0])
        print("file: {}, set {}".format(f,name_marker))
        data = Dataset(input_dir+f)
        var_data = data.variables[variable][:,:,:,:] # time,depth,lon(1),lat
        var_data.mask = True  # create mask for this
        var_data.mask[var_data.data!=0.0]=False
        time = list(range(var_data.shape[0]))
        data.close()
        if  var_data.shape[0]>10:  # otherwise the file is broken, ignore.
            time = [datetime.datetime(year,1,1)+datetime.timedelta(days=int(x)) for x in time]
            time_axis+=time
            #mean temperature
            if do_average:
                var_data_depth = np.mean(var_data,(2,3))  # now time and depth remaining
                                                          # has averages of each depth  
                cells_per_depth = np.sum(~var_data[0,:,:,:].mask,(1,2))
                l_t = layer_thicknessess*cells_per_depth # multiplier for each layer weight
                l_t = np.repeat(l_t[np.newaxis,:], var_data_depth.shape[0],0)
                var_data_depth = var_data_depth*l_t   # each depth average is multiplied by
                                                      # the layer thickness and cell amount
                mean_variable = np.sum(var_data_depth,(1))
                mean_variable = mean_variable/np.sum(l_t[0])
                if total_mean_variable.size == 0:
                    total_mean_variable = mean_variable
                else:
                   total_mean_variable =  \
                        np.concatenate((total_mean_variable,mean_variable))
            # specific points:
            if do_depths:
                for point in specific_depths:
                    point_depth=np.abs(depths-point['depth']).argmin()
                    mean_point = np.mean(var_data[:,point_depth,:,:],(1,2))
#                    mean_point = np.sum(var_data[:,point_depth,:,:],(1,2))/\
#                           (np.sum(~var_data.mask[:,point_depth,:,:],(1,2))).astype(float)
                    point['data']=np.concatenate((point['data'],mean_point))
    if do_average:
        full_data=pd.DataFrame({'time':time_axis,variable:total_mean_variable})
        full_data.to_csv(output_dir+'boundary_mean_{}_{}.csv'.format(variable,name_marker),\
                index=False)
    if do_depths:
        for point in specific_depths:
            full_data=pd.DataFrame({'time':time_axis,variable:point['data']})
            full_data.to_csv(output_dir+'{}_{}_{}.csv'\
                .format(point['name'],variable,name_marker),index=False)
    print(total_mean_variable.shape, len(time_axis))
