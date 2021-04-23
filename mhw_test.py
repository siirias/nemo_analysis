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
import datetime as dt
import time
from netCDF4 import Dataset
import xarray as xr
import pandas as pd
import cmocean as cmo
from smartseahelper import smh


output_dir = "D:\\Data\\Figures\\SmartSea\\mhw\\"
data_dir = "E:\\SmartSea\\all_shark_files\\"  
plot_figures = False
recalculate = True

datasets = [{'n':'A','ref':"A001", '4.5':"A002", '8.5':"A005"},\
            {'n':'B','ref':"B001", '4.5':"B002", '8.5':"B005"},\
            {'n':'D','ref':"D001", '4.5':"D002", '8.5':"D005"},\
            ]
#datasets = [{'n':'A','ref':"A001", '4.5':"A002", '8.5':"A005"},\
#            {'n':'B','ref':"B001", '4.5':"B002", '8.5':"B005"},\
#            ]
#points = ["F64", "SR5", "MS4", "US5B", "F16", "BO3", "F3" ]  # all
#points = ["F64", "SR5", "MS4", "C3", "US5B" ]  # BS
#setup_name = "Bothnian Sea"

#points = ["F16","BO5", "BO3", "F9", "F3" ]  # BoB
#setup_name = "Bay of Bothnia"

cases = {'Bothnian Sea':["F64", "SR5", "MS4", "C3", "US5B" ],
         'Bay of Bothnia':["F16","BO5", "BO3", "F9", "F3" ],
         'Gulf of Bothnia':["F64", "SR5", "MS4", "C3", "US5B", "F16","BO5", "BO3", "F9", "F3"]}

#cases = {'Bothnian Sea':["F64", "SR5" ],
#         'Bay of Bothnia':["F16","BO5" ]}

study_depths = [0.0]
parameter = 'votemper'
group_p = '10AS'


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
    heat_wave_peak_diff = 0.0
    # init the heatwaves. Add one dummy to get statistics start right,
    # even if there is no heatwaves in first year.
    heat_waves = [] 
    for i  in range(len(data.index)):
        if comparison[i] > 0.0:
            if(not in_heat_wave):
                heat_wave_start = data.index[i]
            in_heat_wave = True
            if(comparison[i]>heat_wave_peak_diff):
                heat_wave_peak_diff = comparison[i]
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
                        if(heat_wave_start != data.index[0] and\
                               len(heat_waves) == 0):
                            heat_waves.append({'time':data.index[0], 'length':0})
                            #THis is done to ensure data set starts at day 1.
                        if(heat_wave_length>=min_heat_wave):
                            heat_waves.append({'time':heat_wave_start, \
                                               'length':heat_wave_length,\
                                               'peak_diff':heat_wave_peak_diff})
            if( not in_heat_wave):
                heat_wave_start = None
                heat_wave_length = 0
                break_length = 0
                heat_wave_peak_diff = 0.0
    
    heat_waves.append({'time':data.index[-1], 'length':0})
    return heat_waves


if(recalculate):
    start_time = time.time()
    last_time = start_time
    case_data = {}
    for case in cases:
        points = cases[case]
        setup_name = case
        heatwaves = {'ref':{}, '45':{}, '85':{}}        
        for dataset in datasets:
            heatwaves['ref'][dataset['n']] = {}
            heatwaves['45'][dataset['n']] = {}
            heatwaves['85'][dataset['n']] = {}
            for point in points:
                for study_depth in study_depths:
                    print("set: {}, point: {}, depth: {} m".format(dataset,point,study_depth))
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
    
                    heatwaves['ref'][dataset['n']][point] = gather_extremes(koe_ref, data_ref, parameter)
                    heatwaves['45'][dataset['n']][point] = gather_extremes(koe, data_45, parameter)
                    heatwaves['85'][dataset['n']][point] = gather_extremes(koe, data_85, parameter)
                    print("{},{},{} done, {} sec, total {}, sec".format(\
                         case,dataset['n'],point, \
                         time.time() - last_time, time.time() - start_time))
                    last_time = time.time()
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
                for i in ['ref','45','85']:
                    heatwaves[i][dataset['n']][point] =\
                        pd.DataFrame(heatwaves[i][dataset['n']][point])        
                    heatwaves[i][dataset['n']][point].set_index('time',inplace = True)
        
        # calculate decadal trends:
        decadals = {'ref':{}, '45':{}, '85':{}}        
        decadals['ref'] = {}
        decadals['45'] = {}
        decadals['85'] = {}
        for dataset in datasets:
            decadals['ref'][dataset['n']] = {}
            decadals['45'][dataset['n']] = {}
            decadals['85'][dataset['n']] = {}
            for point in points:
                for i in ['ref','45','85']:
                    tmp = \
                               heatwaves[i][dataset['n']][point].groupby(\
                                pd.Grouper(freq=group_p)).mean()
                    tmp['peak_diff_mean'] = \
                               heatwaves[i][dataset['n']][point].groupby(\
                                pd.Grouper(freq=group_p)).mean()['peak_diff']
                    tmp['peak_diff_min'] = \
                               heatwaves[i][dataset['n']][point].groupby(\
                                pd.Grouper(freq=group_p)).quantile(0.1)['peak_diff']
                    tmp['peak_diff_max'] = \
                               heatwaves[i][dataset['n']][point].groupby(\
                                pd.Grouper(freq=group_p)).quantile(0.9)['peak_diff']
                    
                    tmp['len_std'] = \
                               heatwaves[i][dataset['n']][point].groupby(\
                                pd.Grouper(freq=group_p)).std()['length']
                    tmp['len_75'] = \
                               heatwaves[i][dataset['n']][point].groupby(\
                                pd.Grouper(freq=group_p)).quantile(0.75)['length']
                    tmp['len_25'] = \
                               heatwaves[i][dataset['n']][point].groupby(\
                                pd.Grouper(freq=group_p)).quantile(0.25)['length']


                    tmp['len_high'] = \
                               heatwaves[i][dataset['n']][point].groupby(\
                                pd.Grouper(freq=group_p)).quantile(0.95)['length']
                    tmp['len_low'] = \
                               heatwaves[i][dataset['n']][point].groupby(\
                                pd.Grouper(freq=group_p)).quantile(0.05)['length']
                    decadals[i][dataset['n']][point] = tmp.copy()
                
        # calculate average over all points of interest
        for i in ['ref', '45', '85']:
            decadals[i][dataset['n']] = \
                pd.concat([decadals[i][dataset['n']][x] \
                           for x in decadals[i][dataset['n']]]).\
                           groupby(pd.Grouper(freq=group_p)).mean()        
    #            pd.concat(decadals[i][dataset['n']])\
    #            .groupby(pd.Grouper(freq=group_p)).mean()
        # calculate synthesis:
        for i in ['ref', '45', '85']:
            tmp_all = []
            for the_set in decadals[i]:
                tmp = \
                    pd.concat([decadals[i][the_set][a] for a in decadals[i][the_set]])\
                    .groupby(pd.Grouper(freq=group_p)).mean()
                tmp_all.append(tmp.copy())
            decadals[i]['all'] = pd.concat(tmp_all).groupby(pd.Grouper(freq=group_p)).mean()
    #    decadals['45']['all'] = \
    #        pd.concat([decadals['45'][a] for a in decadals['45']])\
    #        .groupby(pd.Grouper(freq=group_p)).mean()
    #    decadals['85']['all'] = \
    #        pd.concat([decadals['85'][a] for a in decadals['85']])\
    #        .groupby(pd.Grouper(freq=group_p)).mean()
        
        case_data[case] = {'heatwaves':heatwaves, 'decadals':decadals}
    
    #Calculate how many days in year in heatwave:
    for case in case_data:
        case_data[case]['in_hw'] = {}
        for part in ['ref', '45', '85']:
            case_data[case]['in_hw'][part] = {}
            for dataset in datasets:
                case_data[case]['in_hw'][part][dataset['n']] = {}
                for point in case_data[case]['heatwaves'][part][dataset['n']]:
                    dat = case_data[case]['heatwaves'][part][dataset['n']][point]
                    the_days = pd.date_range(dat.index.min(),\
                                  dat.index.max()+\
                                  dt.timedelta(days=\
                                  dat.loc[dat.index.max()]['length']))
                    tmp = pd.DataFrame(zip(the_days,np.zeros(len(the_days))),\
                                 columns = ['date','isinHW'])
                    tmp.set_index('date',inplace=True)
                    
                    #next, lets fill the days actually in heat wave:
                    for hw_ind in case_data[case]['heatwaves'][part][dataset['n']][point].index:
                        hw = case_data[case]['heatwaves'][part][dataset['n']][point].loc[hw_ind]
                        if(hw['length']>0):
                            dates = pd.date_range(hw.name,hw.name+\
                                              dt.timedelta(days=hw['length']))
                            for i in dates:
                                tmp.loc[i] = 1
                    case_data[case]['in_hw'][part][dataset['n']][point] = tmp.copy()
    #sum up the heatwave days:
    for case in case_data:
        case_data[case]['hw_days'] = {}
        for part in ['ref', '45', '85']:
            case_data[case]['hw_days'][part] = {}
            for dataset in datasets:
                
                tmp_all = []
                for point in case_data[case]['in_hw'][part][dataset['n']]:
                    tmp = case_data[case]['in_hw'][part][dataset['n']][point].\
                            groupby(pd.Grouper(freq=group_p)).sum()
    
    
    
                    
                    tmp_all.append(tmp.copy())
    
                tmp = \
                    pd.concat(tmp_all).\
                    groupby(pd.Grouper(freq=group_p)).median()
    
                tmp['days_std'] =  pd.concat(tmp_all).\
                    groupby(pd.Grouper(freq=group_p)).std()['isinHW']
                tmp['days_75'] =  pd.concat(tmp_all).\
                    groupby(pd.Grouper(freq=group_p)).quantile(0.75)['isinHW']
                tmp['days_25'] =  pd.concat(tmp_all).\
                    groupby(pd.Grouper(freq=group_p)).quantile(0.25)['isinHW']
                tmp['days_high'] =  pd.concat(tmp_all).\
                    groupby(pd.Grouper(freq=group_p)).quantile(0.95)['isinHW']
                tmp['days_low'] =  pd.concat(tmp_all).\
                    groupby(pd.Grouper(freq=group_p)).quantile(0.05)['isinHW']
                case_data[case]['hw_days'][part][dataset['n']] = tmp.copy()
    
    for case in case_data:
        for part in ['ref', '45', '85']:
            tmp = pd.concat([case_data[case]['hw_days'][part][x] \
                             for x in case_data[case]['hw_days'][part]]).\
                              groupby(pd.Grouper(freq=group_p)).mean()
            case_data[case]['hw_days'][part]['all'] = tmp.copy()
            
                

#Then some plottings:
for case in case_data:    
    plt.figure()
    plt.title(case)
    decadals = case_data[case]['decadals']
    heatwaves = case_data[case]['heatwaves']
    hw_days = case_data[case]['hw_days']
#    for s in decadals['45'].keys(): #drop the last ones, as they have too little data
#        plt.plot(decadals['ref'][s].index[:-1],decadals['ref'][s].length[:-1],\
#                 color = 'k', marker = '*', alpha = 0.2, linewidth = 0.0)
#        plt.plot(decadals['45'][s].index[:-1],decadals['45'][s].length[:-1],\
#                 color = 'b', marker = '*', alpha = 0.2, linewidth = 0.0)
#        plt.plot(decadals['85'][s].index[:-1],decadals['85'][s].length[:-1],\
#                 color = 'r', marker = '*', alpha = 0.2, linewidth = 0.0)    
    
#    plt.plot(decadals['ref']['all'].index[:-1],decadals['ref']['all'].length[:-1],\
#             color = 'k', marker = '*', alpha = 1.0)
#    plt.plot(decadals['45']['all'].index[:-1],decadals['45']['all'].length[:-1],\
#             color = 'b', marker = '*', alpha = 1.0)
#    plt.plot(decadals['85']['all'].index[:-1],decadals['85']['all'].length[:-1],\
#             color = 'r', marker = '*', alpha = 1.0)    
    
    for i,c in zip(['ref','45','85'],['k','b','r']):
        plt.errorbar(\
                     decadals[i]['all'].index[:-1],\
                     decadals[i]['all'].length[:-1],\
                     yerr =  decadals[i]['all'].len_std[:-1],\
                     elinewidth = 4, linewidth = 0, color = c)
        minmax = np.vstack((decadals[i]['all'].len_high[:-1]-decadals[i]['all'].length[:-1],\
                            decadals[i]['all'].len_low[:-1]-decadals[i]['all'].length[:-1]))        
        minmax = np.abs(minmax)
        plt.errorbar(\
                     decadals[i]['all'].index[:-1],\
                     decadals[i]['all'].length[:-1],\
                     yerr =  minmax,\
                     elinewidth = 1, linewidth = 0, capsize = 5, color = c)

    #Days of heatwave total
    plt.figure()
    plt.title("Days of Heatwaves in 10 years. {}".format(case))
    for i,c in zip(['ref','45','85'],['k','b','r']):
        plt.errorbar(\
                     hw_days[i]['all'].index[:-1],\
                     hw_days[i]['all'].isinHW[:-1],\
                     yerr =  hw_days[i]['all'].days_std[:-1],\
                     elinewidth = 4, linewidth = 0, color = c)
        minmax = np.vstack((hw_days[i]['all'].days_high[:-1]-hw_days[i]['all'].isinHW[:-1],\
                            hw_days[i]['all'].days_low[:-1]-hw_days[i]['all'].isinHW[:-1]))        
        minmax = np.abs(minmax)
        plt.errorbar(\
                     hw_days[i]['all'].index[:-1],\
                     hw_days[i]['all'].isinHW[:-1],\
                     yerr =  minmax,\
                     elinewidth = 1, linewidth = 0, capsize = 5, color = c)
        
        