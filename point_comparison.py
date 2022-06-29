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
import scipy.stats as sp

out_dir = "D:\\Data\\Figures\\SmartSeaNew\\compare\\"
compare_data_dir = "E:\\SmartSea\\all_shark_files\\"  
data_dir = "E:\\SmartSea\\new_dataset\\"

recalculate = True

datasets = [{'n':'A','ref':"A001", '4.5':"A002", '8.5':"A005"},\
            {'n':'B','ref':"B001", '4.5':"B002", '8.5':"B005"},\
            {'n':'D','ref':"D001", '4.5':"D002", '8.5':"D005"},\
            ]

#datasets = [{'n':'A','ref':"A001", '8.5':"A005"},\
#            {'n':'B','ref':"B001", '8.5':"B005"},\
#            {'n':'D','ref':"D001", '8.5':"D005"},\
#            ]

study_depths = [0.0, 100.0]
parameter = 'vosaline'  # 'vosaline', 'votemper'
group_p = '10AS'
points = ['F64', 'SR5', 'US5B', 'BO3']
colors = {'ref':'k', '4.5':'b', '8.5':'r'}

#fig_size = (15,10)
fig_size = (10,15)

def gather_dataset(dataset, part, parameter, study_depth, comparison = False):
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
            tmp.set_index('time',inplace = True)
            tmp = tmp[tmp.index < dt.datetime(2060,12,31)]
        data.append(tmp)
    data = pd.concat(data)
    return data

def calculate_trend(data, parameter):
    y=np.array(data[parameter], dtype=float)
    x=np.array(pd.to_datetime(data.index.values), dtype=float)
    slope, intercept, r_value, p_value, std_err =sp.linregress(x,y)
    xf = np.linspace(min(x),max(x),100)
    xf1 = xf.copy()
    xf1 = pd.to_datetime(xf1)
    yf = (slope*xf)+intercept
                
#    f, ax = plt.subplots(1, 1)
#    ax.plot(xf1, yf,label='Linear fit', lw=3)
#    plt.plot(data.index,data[parameter])
    return slope*24.0*60.0*60.0*1.0e9*3650.0  # trend in decade


trends = {}
compare_trends = {}
for point in points:
    trends[point] = {}
    compare_trends[point] = {}
    for study_depth in study_depths:
        trends[point][study_depth] = {}
        compare_trends[point][study_depth] = {}
        for dataset in datasets:
            for d in dataset:
                if(d != 'n'):
                    trends[point][study_depth][d] = {}
                    compare_trends[point][study_depth][d] = {}

for dataset in datasets:
    for point in points:
        for study_depth in study_depths:
            if(recalculate):
                data = {}
                for d in dataset:
                    if(d != 'n'):
                        data[d] = gather_dataset(dataset, d, parameter, study_depth)
                compare_data = {}
                for d in dataset:
                    if(d != 'n'):
                        compare_data[d] = gather_dataset(dataset, d, parameter, study_depth, comparison = True)
                        print("{} Trend: {}".format(parameter, calculate_trend(data[d],parameter)))
                        trends[point][study_depth][d][dataset['n']] = calculate_trend(data[d],parameter)
                        compare_trends[point][study_depth][d][dataset['n']] = calculate_trend(compare_data[d],parameter)
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
            
            
# Print out the trends
trend_file = open(out_dir + "trends_{}.txt".format(parameter),'w')
datas_options = [ x for x in dataset.keys() if x != 'n']
trend_file.write("NEW RUN, PARAM: {}\n".format(parameter))
for study_depth in study_depths:
    trend_file.write(f"## DEPTH: {study_depth:0.1f}\n")
    for point in points:
        trend_file.write(f"{point}\n")
        for d in datas_options:
            if(d != 'n' and len(trends[point][study_depth][d])>0):
                trend_file.write("{}\n".format(d))
                vals = []
                for dataset in datasets:
                    trend_file.write("\t\t{}: {} //10 years\n".format(dataset['n'],
                                  trends[point][study_depth][d][dataset['n']]))
                    vals.append(trends[point][study_depth][d][dataset['n']])
                mean_val = np.mean(vals)
                trend_file.write("\tEnsemble: {} //10 years\n".format(mean_val))

trend_file.write("*"*20+'\n')
trend_file.write("OLD RUN, PARAM: {}\n".format(parameter))
for study_depth in study_depths:
    trend_file.write(f"## DEPTH: {study_depth:0.1f}\n")
    for point in points:
        trend_file.write(f"{point}\n")
        for d in datas_options:
            if(d != 'n' and len(trends[point][study_depth][d])>0):
                trend_file.write("{}\n".format(d))
                vals = []
                for dataset in datasets:
                    trend_file.write("\t\t{}: {} //10 years\n".format(dataset['n'],
                                  compare_trends[point][study_depth][d][dataset['n']]))
                    vals.append(compare_trends[point][study_depth][d][dataset['n']])
            mean_val = np.mean(vals)
            trend_file.writelines("\tEnsemble: {} //10 years\n".format(mean_val))
trend_file.close()

for i in open(out_dir + "trends_{}.txt".format(parameter)).readlines():
    print(i)
#y=np.array(df['OW2 As(mg/L)'].dropna().values, dtype=float)
#x=np.array(pd.to_datetime(df['OW2 As(mg/L)'].dropna()).index.values, dtype=float)
#slope, intercept, r_value, p_value, std_err =sp.linregress(x,y)
#xf = np.linspace(min(x),max(x),100)
#xf1 = xf.copy()
#xf1 = pd.to_datetime(xf1)
#yf = (slope*xf)+intercept
#print('r = ', r_value, '\n', 'p = ', p_value, '\n', 's = ', std_err)
            
#f, ax = plt.subplots(1, 1)
#ax.plot(xf1, yf,label='Linear fit', lw=3)
#plt.plot(data[d].index,data[d][parameter])
#plt.ylabel('Arsenic concentration')
#ax.legend();
            
##https://mohammadimranhasan.com/linear-regression-of-time-series-data-with-pandas-library-in-python/
