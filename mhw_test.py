# -*- coding: utf-8 -*-
"""
Created on Sat Apr 17 18:52:58 2021

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

from netCDF4 import Dataset
import xarray as xr
import pandas as pd
import cmocean as cmo
from smartseahelper import smh


output_dir = "D:\\Data\\Figures\\SmartSea\\mhw\\"
data_dir = "E:\\SmartSea\\all_shark_files\\"  
plot_figures = False

datasets = [{'n':'A','ref':"B001", '4.5':"B002", '8.5':"B005"}]
points = ["F64", "SR5", "MS4", "US5B", "F16", "BO3", "F3" ]
study_depths = [0.0]
parameter = 'votemper'


def gather_dataset(dataset, part, parameter):
            filenames = os.listdir(os.path.join(data_dir,dataset[part]))
            filenames = [f for f in filenames if f.endswith("shark"+point + '.nc')]
            data = []
            
            for f in filenames:
                # open, and convert to pandas dataframe
                with xr.open_dataset(os.path.join(data_dir,dataset[part],f)) as d:
                    the_depth_arg = np.argmin(np.abs(d.deptht.values-study_depth))
                    tmp = pd.DataFrame({
                            'time':d.time_counter,
                            parameter:d[parameter][:,the_depth_arg,0,0]})
                data.append(tmp)
            data = pd.concat(data)
            data.set_index('time',inplace = True)
            return data

def clim_data_to_axis(reference, clim_data_one):
        clim_data = pd.DataFrame(reference.index)
        clim_data.set_index('time',inplace = True)
        clim_data[parameter] = None
        for i in range(len(clim_data.index)):
            clim_data.loc[clim_data.index[i]][parameter] = \
                     clim_data_one[clim_data_one.index == \
                                   clim_data.index[i].dayofyear]\
                     [parameter].values[0]
        return clim_data

def gather_extremes(treshold, data, parameter):
    max_break = 2
    min_heat_wave = 5
    comparison = np.array(data[parameter]).astype('float') - treshold
    in_heat_wave = False
    break_length = 0
    heat_wave_length = 0
    heat_wave_start = None
    # init the heatwaves. Add one dummy to get statistics start right,
    # even if there is no heatwaves in first year.
    heat_waves = [{'time':data.index[0], 'length':0}] 
    for i  in range(len(data.index)):
        if comparison[i] > 0.0:
            if(not in_heat_wave):
                heat_wave_start = data.index[i]
            in_heat_wave = True
            if(break_length == 0):
                heat_wave_length +=1
            else:
                heat_wave_length +=1+ break_length
            break_length = 0
        else:
            if in_heat_wave:
                if(heat_wave_length<min_heat_wave):
                    in_heat_wave = False
                else:
                    break_length += 1
                    if(break_length > max_break):
                        in_heat_wave = False
                        if(heat_wave_length>=min_heat_wave):
                            heat_waves.append({'time':heat_wave_start, \
                                               'length':heat_wave_length})
            if( not in_heat_wave):
                heat_wave_start = None
                heat_wave_length = 0
                break_length = 0
    return heat_waves
        
heatwaves_ref = {}
heatwaves_45 = {}
heatwaves_85 = {}
for dataset in datasets:
    heatwaves_ref[dataset['n']] = {}
    heatwaves_45[dataset['n']] = {}
    heatwaves_85[dataset['n']] = {}
    for point in points:
        for study_depth in study_depths:
            data_ref = gather_dataset(dataset, 'ref', parameter)
            data_45 = gather_dataset(dataset, '4.5', parameter)
            data_85 = gather_dataset(dataset, '8.5', parameter)
            #calculate the climatology values:
            clim_data_one = data_ref.groupby(data_ref.index.dayofyear).mean() # one year
            clim_data_90_perc = data_ref.groupby(data_ref.index.dayofyear).quantile(0.95) # one year
            clim_data_10_perc = data_ref.groupby(data_ref.index.dayofyear).quantile(0.05) # one year
            # then copy this to cover the whole series:
            clim_data_long =  clim_data_to_axis(data_85,clim_data_one)
            clim_data_trigg_ref = clim_data_to_axis(data_ref, clim_data_90_perc)
            clim_data_trigg = clim_data_to_axis(data_85, clim_data_90_perc)
            clim_data_trigg_low = clim_data_to_axis(data_85, clim_data_10_perc)

            koe_ref = np.array(clim_data_trigg_ref[parameter]).astype('float')
            koe = np.array(clim_data_trigg[parameter]).astype('float')
            koe2 = np.array(clim_data_trigg_low[parameter]).astype('float')
            cd_long = np.array(clim_data_long[parameter]).astype('float')

            heatwaves_ref[dataset['n']][point] = gather_extremes(koe_ref, data_ref, parameter)
            heatwaves_45[dataset['n']][point] = gather_extremes(koe, data_45, parameter)
            heatwaves_85[dataset['n']][point] = gather_extremes(koe, data_85, parameter)
            if(plot_figures):
                fig = plt.figure()
                plt.title("{}:DIFF  depth {}".format(point, study_depth))
                plt.plot(data_45.index, data_45[parameter]-cd_long,'c')
                plt.plot(data_85.index, data_85[parameter]-cd_long,'m')
                plt.plot(data_85.index,np.zeros(len(data_85.index)) ,'k')
                plt.fill_between(data_85.index,\
                                 data_85[parameter]-cd_long, \
                                 koe-cd_long, \
                                 where = (data_85[parameter] > koe),\
                                 color = 'r', alpha = 0.3)
    
                plt.fill_between(data_45.index,\
                                 data_45[parameter]-cd_long, \
                                 koe-cd_long, \
                                 where = (data_45[parameter] > koe),\
                                 color = 'b', alpha = 0.3)
    
                
                fig = plt.figure()
                
                plt.title("{}: depth {}".format(point, study_depth))
                plt.plot(data_45.index, data_45[parameter],'b')
                plt.plot(data_85.index, data_85[parameter],'r')
                plt.fill_between(data_85.index, koe,koe2, alpha = 0.2)
                plt.fill_between(data_85.index,\
                                 data_85[parameter], \
                                 koe, \
                                 where = (data_85[parameter] > koe),\
                                 color = 'r', alpha = 0.3)
    
                plt.fill_between(data_45.index,\
                                 data_45[parameter], \
                                 koe, \
                                 where = (data_45[parameter] > koe),\
                                 color = 'b', alpha = 0.3)
    
                plt.plot(data_85.index, clim_data_long[parameter],'k', alpha = 0.3)


# Gather all data in one nice form:
for dataset in datasets:
    for point in points:
        heatwaves_ref[dataset['n']][point] =\
            pd.DataFrame(heatwaves_ref[dataset['n']][point])        
        heatwaves_45[dataset['n']][point] =\
            pd.DataFrame(heatwaves_45[dataset['n']][point])        
        heatwaves_85[dataset['n']][point] =\
            pd.DataFrame(heatwaves_85[dataset['n']][point])
            
        heatwaves_ref[dataset['n']][point].set_index('time',inplace = True)
        heatwaves_45[dataset['n']][point].set_index('time',inplace = True)
        heatwaves_85[dataset['n']][point].set_index('time',inplace = True)

# calculate decadal trends:
decadals_ref = {}
decadals_45 = {}
decadals_85 = {}
for dataset in datasets:
    decadals_ref[dataset['n']] = []
    decadals_45[dataset['n']] = []
    decadals_85[dataset['n']] = []
    for point in points:
        decadals_ref[dataset['n']].append(\
                   heatwaves_ref[dataset['n']][point].groupby(\
                    pd.Grouper(freq='10AS')).sum())
        decadals_45[dataset['n']].append(\
                   heatwaves_45[dataset['n']][point].groupby(\
                    pd.Grouper(freq='10AS')).sum())
        decadals_85[dataset['n']].append(\
                   heatwaves_85[dataset['n']][point].groupby(\
                    pd.Grouper(freq='10AS')).sum())
    # calculate average over all points of interest
    decadals_ref[dataset['n']] = \
        pd.concat(decadals_ref[dataset['n']])\
        .groupby(pd.Grouper(freq='10AS')).mean()

    decadals_45[dataset['n']] = \
        pd.concat(decadals_45[dataset['n']])\
        .groupby(pd.Grouper(freq='10AS')).mean()

    decadals_85[dataset['n']] = \
        pd.concat(decadals_85[dataset['n']])\
        .groupby(pd.Grouper(freq='10AS')).mean()
        
