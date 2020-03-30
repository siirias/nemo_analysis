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
out_dir = "D:\\Data\\svnfmi_merimallit\\smartsea\\"
fig_size = (15,10)
analyze_inflow = True
analyze_atmosphere = True
analyze_boundary = True
plot_trends = False
b_val = 'vosaline' # 'avg_temp'
period={'min':dt.datetime(1980,1,1), 'max':dt.datetime(2060,1,1)}
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
    in_dir ='D:\\Data\\svnfmi_merimallit\\smartsea\\derived_data\\inflow\\'
    files = os.listdir(in_dir)
    dat={}
    for f in files:
        set_name=re.search('_([^_]*)\.csv',f).groups()[0]
        dat[set_name]=pd.read_csv(in_dir+f,\
                             parse_dates=[0])
        multiplier=1.0
        if set_name == 'hindcast':
            multiplier = 30.5
        print(set_name,dat[set_name]['inflow'].sum()*multiplier)
    plt.figure(figsize=fig_size)
    plt.title("River inflow")
    for s in dat:
        d=dat[s]
        d = d[(d['time']>period['min']) & (d['time']<period['max'])]
        plt.plot(d['time'],d['inflow'], label='_nolegend_', zorder=11,**set_style(s,0.05))
    for s in dat:
        d=dat[s]
        d = d[(d['time']>period['min']) & (d['time']<period['max'])]
        smooth_window = 1000
        if(s=='hindcast'):
            smooth_window /= 30.5 # as it is monthly, not daily
            smooth_window = max(1,smooth_window)
        smoothed = d['inflow'].ewm(span = smooth_window,min_periods=smooth_window).mean()
        fitting_time = mp.dates.date2num(d['time'])
        fitting = np.polyfit(fitting_time,d['inflow'],1)
        print("{} change: {:.3} unit/year".format(s,fitting[0]*365.15))
        plt.plot(d['time'],smoothed,label=s, zorder=15,**set_style(s))
        if(plot_trends):
            plt.plot(mp.dates.num2date(fitting_time),fitting[0]*fitting_time+fitting[1],label='_nolegend_', zorder=15,**set_style(s))
    plt.legend()
    plt.savefig(out_dir+"inflow_{}_{}-{}.png".format(\
                set_name,period['min'].year,period['max'].year))
    


if analyze_boundary:
#    in_dir ='D:\\Data\\svnfmi_merimallit\\smartsea\\derived_data\\boundaryjusthistory\\'
    in_dir ='D:\\Data\\svnfmi_merimallit\\smartsea\\derived_data\\boundary\\'
#    for subset in ['boundary_mean','5meter','20meter','80meter','120meter']:
    for subset in ['boundary_mean','5meter','80meter']:
        files = os.listdir(in_dir)
        files = [f for f in files if subset in f]
        dat={}
        for f in files:
            set_name=re.search('_([^_]*)\.csv',f).groups()[0]
            dat[set_name]=pd.read_csv(in_dir+f,\
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
        plt.savefig(out_dir+"boundary_{}_{}_{}-{}.png".format(\
                    subset,set_name,period['min'].year,period['max'].year))
    
