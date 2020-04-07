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
#import seaborn as sns
import os
import re
import netCDF4
from netCDF4 import Dataset

out_dir = "D:\\Data\\SmartSeaModeling\\Images\\"
fig_size = (15,10)
analyze_salt_content = False
analyze_salt_profiles = False
analyze_salt_trends = True
plot_trends = True

trend_file_name = out_dir+'point_trends_salinity.csv'
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

class ValueSet():
    def __init__(self):
        self.data = {}
    def add(self, point, lat, lon, depth, set_name, value):
        if(not point in self.data.keys()):
            self.data[point] = {}
            self.data[point]['lat'] = lat
            self.data[point]['lon'] = lon
        if(not depth in self.data[point].keys()):
            self.data[point][depth] = {}
        if(not set_name in self.data[point][depth].keys()):
            self.data[point][depth][set_name] = value
        return True

    def give_values(self,point,depth,filter_str=".*"):
        all_sets = self.data[point][depth].keys()
        the_sets = [i for i in all_sets if re.match(filter_str,i)]
        return [self.data[point][depth][i] for i in the_sets]

    def mean(self,point,depth,filter_str=".*"):
        return np.mean(self.give_values(point,depth,filter_str))
    def max(self,point,depth,filter_str=".*"):
        return np.max(self.give_values(point,depth,filter_str))
    def min(self,point,depth,filter_str=".*"):
        return np.min(self.give_values(point,depth,filter_str))
        
if analyze_salt_content:
    variable = 'sea_water_absolute_salinity'
    in_dir ='D:\\Data\\SmartSeaModeling\\'
    name_format = 'reserve_salinity_(.*)_total.nc'
    files = os.listdir(in_dir)
    dat={}
    for f in files:
        
        set_name=re.search(name_format,f)
        if(set_name):
            set_name = set_name.groups()[0]
#        dat[set_name]=pd.read_csv(in_dir+f,\
#                             parse_dates=[0])
            print(set_name)
            D = Dataset(in_dir+f)
            values = D['sea_water_absolute_salinity']
            times = D['time']
            times = netCDF4.num2date(times[:],times.units)
            dat[set_name] = pd.DataFrame({'time':times, 'value':values})
    plt.figure(figsize=fig_size)
    plt.title("Total amount of salt in GoB (GT)")
    for s in dat:
        d=dat[s]
        d = d[(d['time']>period['min']) & (d['time']<period['max'])]
        plt.plot(d['time'],d['value'], label='_nolegend_', zorder=11,**set_style(s,0.2))
    for s in dat:
        d=dat[s]
        d = d[(d['time']>period['min']) & (d['time']<period['max'])]
        smooth_window = 12 #yearly 
        smoothed = d['value'].ewm(span = smooth_window,min_periods=smooth_window).mean()
        fitting_time = mp.dates.date2num(d['time'])
        fitting = np.polyfit(fitting_time,d['value'],1)
        print("{} change: {:.3} unit/year".format(s,fitting[0]*365.15))
        label_text = "{}:{:0.3} unit/year".format(s,fitting[0]*365.15)
        plt.plot(d['time'],smoothed,label=label_text, zorder=15,**set_style(s))
        if(plot_trends):
            plt.plot(mp.dates.num2date(fitting_time),fitting[0]*fitting_time+fitting[1],label='_nolegend_', zorder=15,**set_style(s,0.4))
    plt.legend()
    plt.savefig(out_dir+"total_salt_{}-{}.png".format(\
                period['min'].year,period['max'].year))
gathered_profile_trends = ValueSet()
if analyze_salt_profiles:
    variable = 'vosaline'
    all_depths = [0.0,50.0,2000.0] #depth, if under the bottom, the lowest with number is accepted.
    points = ['F64', 'SR5', 'MS4', 'C3', 'US5B', 'F16', 'BO3', 'F3', 'F9', 'BO5']
    fixed_axis= None #[2.0,9.0] #None or [min, max]
    for point in points:
        for depth_in in all_depths:
            in_dir ='D:\\Data\\SmartSeaModeling\\Extracted_profiles\\'
            name_format = 'profile_{}_(.*)_vosaline.nc'.format(point)
            files = os.listdir(in_dir)
            files = [i for i in files if re.match(name_format,i)]
            dat={}
            for f in files:
                set_name=re.search(name_format,f)
                if(set_name):
                    set_name = set_name.groups()[0]
                    print(set_name)
                    D = Dataset(in_dir+f)
                    values = D[variable]
                    times = D['date']
                    times_orig = times[:]
                    lat = float(D['latitude'].getValue())
                    lon = float(D['longitude'].getValue())
                    depths = D['deptht']
                    max_depth = depths[values[0,:][values[0,:].mask == False]\
                                       .shape[0]-1]
                    depth = float(depths[np.abs((depths[:]-depth_in)).argmin()])
                    depth = min(depth,max_depth)
                    depth_layer = np.abs(np.array(depths)-depth).argmin()
                    times = netCDF4.num2date(times[:],times.units)
                    dat[set_name] = pd.DataFrame({'time':times,\
                                   'value':values[:,depth_layer],\
                                   'lat':lat,
                                   'lon':lon})
            plt.figure(figsize=fig_size)
            plt.title("Salinity on {} depth {:0.1f} m (Max Depth {:0.0f} m)"\
                          .format(point,depth,max_depth))
            for s in dat:
                d=dat[s]
                d = d[(d['time']>period['min']) & (d['time']<period['max'])]
                plt.plot(d['time'],d['value'], label='_nolegend_',\
                                 zorder=11,**set_style(s,0.2))
        
                smooth_window = 12 #yearly 
                smoothed = d['value'].ewm(span = smooth_window,\
                                min_periods=smooth_window).mean()
                fitting_time = mp.dates.date2num(d['time'])
                fitting = np.polyfit(fitting_time,d['value'],1)
                print("{} change: {:.3} unit/year".format(s,fitting[0]*365.15))
                label_text = "{}:{:0.3f} u/dec".format(s,fitting[0]*3651.5)
                plt.plot(d['time'],smoothed,label=label_text, \
                             zorder=15,**set_style(s))
                gathered_profile_trends.add(\
                        point,\
                        d['lat'].iloc[0],\
                        d['lon'].iloc[0],\
                        "{:0.1f}".format(depth),\
                        s,\
                        fitting[0]*365.15)
                if(plot_trends):
                    plt.plot(mp.dates.num2date(fitting_time),\
                             fitting[0]*fitting_time+fitting[1],\
                             label='_nolegend_', zorder=15,**set_style(s,0.4))
            plt.legend()
            if(fixed_axis):
                plt.ylim(fixed_axis[0],fixed_axis[1])
            print("saving",depth,point)
            plt.savefig(out_dir+"\\profiles\\"+\
                        "Salinity_profile_{}_{:0.1f}m_{}-{}.png".format(\
                        point,\
                        depth,\
                        period['min'].year,\
                        period['max'].year))
    #write trend analysis
    with open(trend_file_name,'w') as out_f:
        out_f.write("Point\tlat\tlon\tdepth\tscenario\tmean\tmin\tmax\n")
        for fil,tag in zip(['.*1','.*2','.*5'],\
                           ['HISTORY','RCP4.5','RCP8.5']):
            for point in gathered_profile_trends.data.keys():
                for depth in gathered_profile_trends.data[point].keys():
                    try:
                        depth_f=float(depth)
                        ok = True
                    except:
                        ok = False
                    if(ok):
                        mean_val = gathered_profile_trends.mean(point,depth,fil)
                        max_val = gathered_profile_trends.max(point,depth,fil)
                        min_val = gathered_profile_trends.min(point,depth,fil)
                        lat = gathered_profile_trends.data[point]['lat']
                        lon = gathered_profile_trends.data[point]['lon']
                        print(\
                        "{}, {} m {}: mean {:0.3f} (min {:0.3f}, max {:0.3f})".format(\
                         point, depth_f, tag,  mean_val, min_val, max_val))
                        out_f.write("{}\t{:0.2f}\t{:0.2f}\t{}\t{}\t{:0.3f}\t{:0.03f}\t{:0.03f}\n".format(\
                                  point, lat, lon, depth_f, tag, \
                                  mean_val, min_val, max_val))
if analyze_salt_trends:
    #open just saved file as pandas, and do some plotting
    scenarios = ['HISTORY','RCP4.5','RCP8.5']
    for scenario in scenarios:
        
        depths = [1.5, 50.0]
        depth_vars = [0.5, 10.0]
        for depth, depth_var in zip(depths,depth_vars):
            dataf = pd.read_csv(trend_file_name,sep='\t')
            d = dataf[dataf['scenario'] == scenario]
            d = d[d['depth'] > depth - depth_var]
            d = d[d['depth'] < depth + depth_var]
            d = d.sort_values('lat')
            figure = plt.figure(figsize=fig_size)
            plt.title("Salinity trend {} depth {:0.1f} m".format(scenario,depth))
            plt.plot(d['lat'],d['mean'],'b*')
            plt.plot(d['lat'],[0]*len(d['lat']),'k',alpha=0.3)
            axis = figure.axes[0]
            axis.fill_between(d['lat'],d['max'],d['min'],\
                              facecolor = 'b', alpha=0.2)
            for point,lat,val in zip(d['Point'],d['lat'],d['mean']):
                plt.text(lat,val,point)
            #plt.ylim(-0.02,0.04)
            plt.savefig(out_dir+\
                        "Salinity_trends_{}_{:0.1f}m.png".format(\
                        scenario, depth))