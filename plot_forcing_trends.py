# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import datetime as dt
import matplotlib as mp
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
import re
output_dir = "D:\\Data\\Figures\\SmartSea\\Forcings\\"
main_data_dir = "D:\\Data\\svnfmi_merimallit\\smartsea\\derived_data\\"

fig_size = (8,5)
fig_dpi = 300
analyze_inflow = True
analyze_atmosphere = True
analyze_boundary = False
plot_trends = False
plot_smoothed = False
plot_original = False
plot_cloud = True
b_val = 'vosaline' # 'avg_temp'
period={'min':dt.datetime(1980,1,1), 'max':dt.datetime(2060,1,1)}
change_time = dt.datetime(2006,1,1) # used to cut the forecasts before this
def set_style(set_name,alpha=1.0):
    scenario = ""
    if( '1' in set_name or 'hindcast' in set_name):
        scenario = "history"
    if('2' in set_name):
        scenario = "rcp45"
    if('5' in set_name):
        scenario = "rcp85"
    model_type = re.search("^.",set_name).group()
    colors = {
        'A':'b',
        'B':'r',
        'C':'c',
        'D':'g',
        'h':'k'
    }
    scen_styles = {
        'history':'--',
        'rcp45':'-',
        'rcp85':'-'
    }
    line_width = 1.0
    if(scenario == 'rcp85'):
        line_width=2.0
    return {'color':colors[model_type],
            'linestyle':scen_styles[scenario],
            'linewidth':line_width,
            'alpha':alpha}
    
    
        
if analyze_inflow:
    # water inflow is originaly in forcings kg/(m^2 sec)
    # csv's have summed them over and multiplied by each cell's area,
    # so number is kg/sec of water for the whole area.
    inflow_numbers = []
    data_dir = main_data_dir + '\\inflow\\'
    files = os.listdir(data_dir)
    files = [x for x in files if x.endswith('csv')]
    dat={}
    for f in files:
        set_name=re.search('_([^_]*)\.csv',f).groups()[0]
        dat[set_name]=pd.read_csv(data_dir+f,\
                             parse_dates=[0])
        dat[set_name]['inflow'] = dat[set_name]['inflow']*\
                                1000000\
                                *60*60*24*365\
                                *0.0001*0.0001*0.0001  
                                #fixes one eror in csv crations, then
                                # changes unit from kg per second
                                # into km^3/year
        dat[set_name]=dat[set_name].set_index('time')
        multiplier=1.0
        if set_name == 'hindcast':
            multiplier = 30.5
        print(set_name,dat[set_name]['inflow'].sum()*multiplier)
    fig = plt.figure(figsize=fig_size)
    plt.title("River inflow")
    for s in dat:
        d=dat[s]
        if('2' in s or '5' in s):
            period_min = change_time
            period_max = period['max']
        else:
            period_min = period['min']
            period_max = change_time
        d = d[(d.index>period_min) & (d.index<period_max)]
        if(plot_original):
            plt.plot(d['time'],d['inflow'], label='_nolegend_', zorder=11,**set_style(s,0.05))
        smooth_window = 6000
        if(s=='hindcast'):
            smooth_window /= 30.5 # as it is monthly, not daily
            smooth_window = max(1,smooth_window)
        smoothed = d['inflow'].ewm(span = smooth_window,min_periods=smooth_window).mean()
        fitting_time = mp.dates.date2num(d.index)
        fitting = np.polyfit(fitting_time,d['inflow'],1)
        print("{} change: {:.3} unit/year".format(s,fitting[0]*365.15))
        if(plot_smoothed):
            plt.plot(d.index,smoothed,label=s, zorder=15,**set_style(s))
#        plt.gcf().autofmt_xdate()
        s_cloud = set_style(s)
        s_cloud['alpha'] = 0.1
        d_tmp = d.groupby(pd.Grouper(freq='1Y')).mean()
        mean = d_tmp.groupby(pd.Grouper(freq='10Y')).mean()
        median = d_tmp.groupby(pd.Grouper(freq='10Y')).median()
        std = d_tmp.groupby(pd.Grouper(freq='10Y')).std()
        print("Mean std for {}: {}".format(s,std.mean()))
        if(plot_cloud):
            plt.plot(median.index,median['inflow'], label=s,zorder=16,**set_style(s))
            plt.fill_between(median.index,\
                             mean['inflow']-std['inflow'],\
                             mean['inflow']+std['inflow'],
                             **s_cloud)
            
        if(plot_trends):
            plt.plot(mp.dates.num2date(fitting_time),\
                     fitting[0]*fitting_time+fitting[1],\
                     label='_nolegend_', zorder=15,**set_style(s))
        inflow_numbers.append([s,float(mean.mean()),fitting[0]*365.15,float(std.mean())])
#    plt.ylim(4,10)
    plt.xlabel('Year')
    plt.ylabel('Average kg^3 Yearly inflow')
    plt.legend()
    plt.savefig(output_dir+"inflow_{}_{}-{}.png".format(\
                set_name,period['min'].year,period['max'].year),dpi = fig_dpi)
    
print(pd.DataFrame(inflow_numbers,columns=['set','mean','trend','std']).to_latex())

if analyze_boundary:
#    data_dir ='D:\\Data\\svnfmi_merimallit\\smartsea\\derived_data\\boundaryjusthistory\\'
    data_dir = main_data_dir + '\\boundary\\'
#    for subset in ['boundary_mean','5meter','20meter','80meter','120meter']:
    for subset in ['boundary_mean','5meter','80meter']:
        files = os.listdir(data_dir)
        files = [f for f in files if subset in f]
        dat={}
        for f in files:
            set_name=re.search('_([^_]*)\.csv',f).groups()[0]
            dat[set_name]=pd.read_csv(data_dir+f,\
                                 parse_dates=[0])
        plt.figure(figsize=fig_size)
        plt.title("Boundary Salinity, {}".format(subset))
        for s in dat:
            d=dat[s]
            d = d[(d['time']>period['min']) & (d['time']<period['max'])]
            plt.plot(d['time'],d[b_val], label='_nolegend_', zorder=11,**set_style(s,0.2))
    
            smooth_window = 366*2
            smoothed = d[b_val].ewm(span = smooth_window, min_periods=smooth_window).mean()
            plt.plot(d['time'],smoothed, zorder=15,label='_nolegend_',**set_style(s))

            fitting_time = mp.dates.date2num(d['time'])
            fitting = np.polyfit(fitting_time,d[b_val],1)
            print("{} {} change: {:.3} (g/g)/year".format(s,subset,fitting[0]*365.15))
            plt.plot(d['time'],smoothed,label=s, zorder=15,**set_style(s))
            if(plot_trends):
                plt.plot(mp.dates.num2date(fitting_time),fitting[0]*fitting_time+fitting[1],label='_nolegend_', zorder=15,**set_style(s,0.5))

        plt.ylim([4.0,9.0])
        plt.legend()
        plt.savefig(output_dir+"boundary_{}_{}_{}-{}.png".format(\
                    subset,set_name,period['min'].year,period['max'].year),dpi = fig_dpi)
    
