# -*- coding: utf-8 -*-
"""
Created on Sat Apr 17 18:52:58 2021

Quick test to compare single measurement point data between two setups
@author: siirias
"""


import sys
import os
import re
import matplotlib as mp
import matplotlib.pyplot as plt
import numpy as np
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import datetime as dt
import time
from netCDF4 import Dataset
import xarray as xr
import pandas as pd
import cmocean as cmo
from smartseahelper import smh
from collections import Counter # Needed for the histogram stuff


out_dir = "D:\\Data\\Figures\\SmartSeaNew\\compare\\"
compare_data_dir = "E:\\SmartSea\\all_shark_files\\"  
data_dir = "E:\\SmartSea\\new_dataset\\"

recalculate = True

#datasets = [{'n':'A','ref':"A001", '4.5':"A002", '8.5':"A005"},\
#            {'n':'B','ref':"B001", '4.5':"B002", '8.5':"B005"},\
#            {'n':'D','ref':"D001", '4.5':"D002", '8.5':"D005"},\
#            ]

datasets = [{'n':'A','ref':"A001", '8.5':"A005"},\
            {'n':'B','ref':"B001", '8.5':"B005"},\
            {'n':'D','ref':"D001", '8.5':"D005"},\
            ]

study_depths = [0.0, 100.0]
parameter = 'vosaline'  # 'vosaline', 'votemper'
group_p = '10AS'
points = ['F64', 'SR5', 'US5B', 'BO3']
colors = {'ref':'k', '4.5':'b', '8.5':'r'}

fig_size = (15,10)

def gather_dataset(dataset, part, parameter, comparison = False):
    if(comparison):    
        tmp_data_dir = compare_data_dir
    else:
        tmp_data_dir = data_dir
    filenames = os.listdir(os.path.join(tmp_data_dir,dataset[part]))
    filenames = [f for f in filenames if f.endswith("shark"+point + '.nc')]
    data = []
    
    for f in filenames:
        # open, and convert to pandas dataframe
        with xr.open_dataset(os.path.join(tmp_data_dir,dataset[part],f)) as d:
            the_depth_arg = np.argmin(np.abs(d.deptht.values-study_depth))
            tmp = pd.DataFrame({
                    'time':d.time_counter,
                    parameter:d[parameter][:,the_depth_arg,0,0]})
        data.append(tmp)
    data = pd.concat(data)
    data.set_index('time',inplace = True)
    return data


for dataset in datasets:
    for point in points:
        for study_depth in study_depths:
            if(recalculate):
                data = {}
                for d in dataset:
                    if(d != 'n'):
                        data[d] = gather_dataset(dataset, d, parameter)
                compare_data = {}
                for d in dataset:
                    if(d != 'n'):
                        compare_data[d] = gather_dataset(dataset, d, parameter, comparison = True)
                print("Calculated: {},{},{} m".format(dataset, point, study_depth))
                
            
            fig, axes_list = plt.subplots(2,1,figsize = fig_size)
            plt.axes(axes_list[0])
            plt.title("{}, {}  {} m\n {}".format(dataset['n'],parameter, study_depth, point))
            for d in data:
                plt.plot(data[d].index, data[d][parameter],colors[d])
            for d in compare_data:
                plt.plot(compare_data[d].index, compare_data[d][parameter],colors[d], alpha = 0.3)
            plt.grid()
        
        # DIFF FIGURE        
            plt.axes(axes_list[1])
            plt.title("Diff, {}, {} {} m\n {}".format(dataset['n'],parameter, study_depth, point))
            for d in data:
                tmp_dat = data[d]- compare_data[d]
                plt.plot(tmp_dat.index, tmp_dat[parameter],colors[d])
                plt.fill_between(tmp_dat.index, tmp_dat[parameter],color = colors[d],alpha = 0.2)
            plt.grid()
            out_filename = "{}_{}_{}_{}m_comparison.png".format(\
                    dataset['n'], parameter, point, int(study_depth))
            plt.savefig(out_dir+out_filename, dpi = 300, bbox_inches='tight')
            print("saved: {} {}".format(out_dir, out_filename))
                            