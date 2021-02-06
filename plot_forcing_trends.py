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
import smartseahelper
sm = smartseahelper.smh()


output_dir = "D:\\Data\\Figures\\SmartSea\\Forcings\\"
main_data_dir = "D:\\Data\\svnfmi_merimallit\\smartsea\\derived_data\\"

fig_factor = 0.8  #1.5
#fig_size = (10*fig_factor,5*fig_factor)
fig_size = (10*fig_factor,7*fig_factor)
fig_dpi = 300
analyze_inflow = True
analyze_atmosphere = True
analyze_boundary = False

plot_single_models = False
plot_combinations = True

plot_trends = True
plot_smoothed = False
plot_yearly_mean = True
plot_original = False
plot_cloud = False
plot_scatter = True
show_grid = True
fix_inflow_ylims = False # (100,350)# False
show_trends_in_Legend = False

b_val = 'vosaline' # 'avg_temp'
#period={'min':dt.datetime(1980,1,1), 'max':dt.datetime(2060,1,1)}
period={'min':dt.datetime(1976,1,1), 'max':dt.datetime(2060,1,1)}
change_time = dt.datetime(2006,1,1) # used to cut the forecasts before this
   
plot_shift = dt.timedelta(5*365)  # how much decadal errorbars are shifted to middle of the decade
extra_shift_step = dt.timedelta(0.2*365) # keep the errorbars from overlapping (too much)
    
        
if analyze_inflow:
    # water inflow is originaly in forcings kg/(m^2 sec)
    # csv's have summed them over and multiplied by each cell's area,
    # so number is kg/sec of water for the whole area.
    variable = 'inflow'
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
    #calculate the means for History, RCP4.5 and RCP8.5
    if(plot_combinations):
        dat["Control"] = pd.concat([dat['A001'],dat['B001'],dat['D001']])
        dat["RCP45"] = pd.concat([dat['A002'],dat['B002'],dat['D002']])
        dat["RCP85"] = pd.concat([dat['A005'],dat['B005'],dat['D005']])
        if(not plot_single_models): # remove the A,B,D thingies from the list
            for i in list(dat.keys()):
                if(i.startswith('A') or i.startswith('B') or i.startswith('D')):
                    dat.pop(i)
    
    fig = plt.figure(figsize=fig_size)
    plt.title("River inflow")
    extra_shift = -extra_shift_step*3.0  # used to shift whisker plots a bit
    for s in dat:
        d=dat[s]
        if('2' in s or '5' in s):
            period_min = change_time
            period_max = period['max']
        else:
            period_min = period['min']
            period_max = change_time
        d = d[(d.index>period_min) & (d.index<period_max)]
        

        smooth_window = 6000
        if(s=='hindcast'):
            smooth_window /= 30.5 # as it is monthly, not daily
            smooth_window = max(1,smooth_window)
        smoothed = d[variable].ewm(span = smooth_window,min_periods=smooth_window).mean()
        fitting_time = mp.dates.date2num(d.index)
        fitting = np.polyfit(fitting_time,d[variable],1)
        print("{} change: {:.3} km^3/yr".format(s,fitting[0]*365.15))
        s_temp = s
        if(s == 'hindcast'):
            s_temp = 'Hindcast'
        if(show_trends_in_Legend):
            label_text = "{}:{:0.3} km^3/yr".format(s_temp,fitting[0]*365.15)
        else:
            label_text = "{}".format(s_temp)

        if(plot_original):
            plt.plot(d.index,d['inflow'], label='_nolegend_', zorder=11,**sm.set_style(s,0.05))
        if(plot_smoothed):
            plt.plot(d.index,smoothed,label=label_text, zorder=15,**sm.set_style(s))
            label_text = None
#        plt.gcf().autofmt_xdate()
        s_cloud = sm.set_style(s)
        s_cloud['alpha'] = 0.1
        s_cloud.pop('marker') # get rid of the marker, as fill_between doesn't have one
        d_tmp = d.groupby(pd.Grouper(freq='1AS')).mean()
        mean = d_tmp.groupby(pd.Grouper(freq='10AS')).mean()
        median = d_tmp.groupby(pd.Grouper(freq='10AS')).median()
        std = d_tmp.groupby(pd.Grouper(freq='10AS')).std()
        maximum = d_tmp.groupby(pd.Grouper(freq='10AS')).max()
        minimum = d_tmp.groupby(pd.Grouper(freq='10AS')).min()
        quant_min = d_tmp.groupby(pd.Grouper(freq='10AS')).quantile(0.75)
        quant_max = d_tmp.groupby(pd.Grouper(freq='10AS')).quantile(0.25)
        print("Mean std for {}: {}".format(s,std.mean()))
        
        if(plot_scatter):
            plot_shift_plus = plot_shift + extra_shift
            scatter_style = sm.set_style(s)
            scatter_style['marker'] = 'D'
            scatter_style['s'] = scatter_style['linewidth']*30
            scatter_style['linewidth'] = 0.0
            plt.scatter(median.index+plot_shift_plus,median[variable], \
                        label=label_text, zorder=16,**scatter_style)
            label_text = None # to prevent plotting the label more than once
            scatter_style.pop('s')
            scatter_style['marker'] = ''
            scatter_style['linestyle'] = ' '
            scatter_style['elinewidth'] = 3
#                    scatter_style['capsize'] = 5
            
            minmax = np.vstack((mean[variable]-quant_max[variable],\
                                quant_min[variable]- mean[variable]))
            plt.errorbar(median.index+plot_shift_plus,mean[variable], \
                         yerr = minmax,\
                        label=label_text, zorder=16,**scatter_style)

            minmax = np.vstack((mean[variable]-minimum[variable],\
                                maximum[variable]- mean[variable]))
            scatter_style['elinewidth'] = 1
            scatter_style['capsize'] = 3
            plt.errorbar(median.index+plot_shift_plus,mean[variable], \
                         yerr = minmax,\
                        label=label_text, zorder=16,**scatter_style)
            extra_shift += extra_shift_step
        
        if(plot_cloud):
            plt.plot(median.index,median[variable], label=label_text,zorder=16,**sm.set_style(s))
            plt.fill_between(median.index,\
                             mean[variable]-std[variable],\
                             mean[variable]+std[variable],
                             **s_cloud)
            label_text = None
            
        if(plot_trends):
            plt.plot(mp.dates.num2date(fitting_time),\
                     fitting[0]*fitting_time+fitting[1],\
                     label='_nolegend_', zorder=15,**sm.set_style(s))
        if(plot_yearly_mean):
            mean_style = sm.set_style(s)
            mean_style['alpha'] = 0.15
            plt.plot(d_tmp.index,d_tmp[variable], label=label_text,zorder=16,**mean_style)
            label_text = None # to prevent plotting the label more than once
            
        inflow_numbers.append([s,float(mean.mean()),fitting[0]*365.15,float(std.mean())])
    if(fix_inflow_ylims):
        plt.ylim(*fix_inflow_ylims)
    plt.xlim([period['min'],period['max']])
    plt.xlabel('Year')
    plt.ylabel('Average km^3 Yearly inflow')
    if(show_grid):
        plt.grid('on')
    plt.legend()
    extra = ""
    if(plot_combinations):
        extra+="comb"
    plt.savefig(output_dir+"inflow_{}_{}-{}{}.png".format(\
                set_name,\
                period['min'].year,period['max'].year,\
                extra),dpi = fig_dpi)
    
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
            d = d[(d.index>period['min']) & (d.index<period['max'])]
            plt.plot(d.index,d[b_val], label='_nolegend_', zorder=11,**sm.sm.set_style(s,0.2))
    
            smooth_window = 366*2
            smoothed = d[b_val].ewm(span = smooth_window, min_periods=smooth_window).mean()
            plt.plot(d.index,smoothed, zorder=15,label='_nolegend_',**sm.set_style(s))

            fitting_time = mp.dates.date2num(d.index)
            fitting = np.polyfit(fitting_time,d[b_val],1)
            print("{} {} change: {:.3} (g/g)/year".format(s,subset,fitting[0]*365.15))
            plt.plot(d.index,smoothed,label=s, zorder=15,**sm.set_style(s))
            if(plot_trends):
                plt.plot(mp.dates.num2date(fitting_time),fitting[0]*fitting_time+fitting[1],label='_nolegend_', zorder=15,**sm.set_style(s,0.5))

        plt.ylim([4.0,9.0])
        plt.legend()
        plt.savefig(output_dir+"boundary_{}_{}_{}-{}.png".format(\
                    subset,set_name,period['min'].year,period['max'].year),dpi = fig_dpi)
    
