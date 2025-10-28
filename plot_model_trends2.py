# -*- coding: utf-8 -*-
"""
Spyder Editor

Modification of plot_model_trends.
new dataset, with slightly different names.
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
import smartseahelper
from netCDF4 import Dataset
import xarray as xr
from pandas.plotting import register_matplotlib_converters
register_matplotlib_converters()

import warnings
warnings.filterwarnings("ignore")

#out_dir = "D:\\Data\\SmartSeaModeling\\Images\\"
sm = smartseahelper.smh()
sm.root_data_in = "D:\\SmartSea\\new_dataset\\"
sm.root_data_out = "C:\\Data\\"
#sm.root_data_in = "D:\\Data\\svnfmi_merimallit\\smartsea\\"
out_dir = sm.root_data_out+"figures\\SmartSeaNEW\\test\\"
fig_factor = 1.5#1.5 #0.8  #1.5
fig_size = (10*fig_factor,5*fig_factor)
#analyze_salt_content = True
#analyze_heat_content = True
#content_types = {"analyze_salt_content":True, "analyze_heat_content":True}
content_types = {"analyze_salt_content":True,\
                 "analyze_heat_content":True}

analyze_profiles = True
profile_types = ["vosaline", "votemper"]
#profile_types = ["vosaline"]
analyze_salt_trends = True
analyze_sbs_changes = True
analyze_correlations = True
plot_single_models = True
plot_combinations = not plot_single_models
model_area = 5959.7#6286  #km^3

plot_original = True
plot_yearly_mean = True
plot_smoothed = False
plot_trends = False
plot_cloud = False
plot_scatter = True
show_grid = True
use_total_salt_amount = False # Total amount, or average salinity.
use_total_heat_energy = False # Heat energy, or average temperature.

plot_shift = dt.timedelta(5*365)  # how much decadal errorbars are shifted to middle of the decade
extra_shift_step = dt.timedelta(0.2*365) # keep the errorbars from overlapping (too much)

create_ensembles = True
ensemble_filters = {'RCP45':'002','RCP85':'005','HISTORY':'001'}

drop_hindcast = False
#period={'min':dt.datetime(2006,1,1), 'max':dt.datetime(2100,1,1)}
period={'min':dt.datetime(1980,1,1), 'max':dt.datetime(2100,1,1)}
#period={'min':dt.datetime(1980,1,1), 'max':dt.datetime(2060,1,1)}
#period={'min':dt.datetime(2006,1,1), 'max':dt.datetime(2060,1,1)}
#period={'min':dt.datetime(1980,1,1), 'max':dt.datetime(2006,1,1)}

if(period['min'] >= dt.datetime(2006,1,1)):
    drop_hindcast = True

def make_ensemble(data_sets, ensemble_string, param = 'value'):
    keys = [x for x in data_sets.keys() if ensemble_string in x]    
    ensemble_vals = np.mean([data_sets[x][param] for x in keys],0)
    ensemble = data_sets[keys[0]].copy()
    ensemble[param] = ensemble_vals
    return ensemble


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
        return_value = [self.data[point][depth][i] for i in the_sets]
        if len(return_value) == 0:
            return [0]
        else:
            return return_value

    def mean(self,point,depth,filter_str=".*"):
        return np.mean(self.give_values(point,depth,filter_str))
    def max(self,point,depth,filter_str=".*"):
        return np.max(self.give_values(point,depth,filter_str))
    def min(self,point,depth,filter_str=".*"):
        return np.min(self.give_values(point,depth,filter_str))


if(analyze_correlations): 
    boundary_data = sm.load_boundary_data()

#    
#    
# Plots conserning the whole Model area    
#        
#    
#    
    
for a in content_types:
    data_multiplier = 1.0 # gludge to change from total salt to salinity.
    data_zero_point = 0.0 # gludge to deal with Celcius/Kelvin
    if content_types[a]:
        if(a == "analyze_salt_content"):
            variable = 'total_salt'
            name_format = 'reserve_vosaline_(.*)\.nc'
            if use_total_salt_amount:
                title_text = "Total amount of salt in GoB (GT)"
                trend_unit = "GT/decade"
            else: #calculate average salinity
                title_text = "Average salinity over model area (g/kg)"
                trend_unit = "(g/kg)/decade"
                data_multiplier = 1./model_area
            
        elif(a == "analyze_heat_content"):
            variable = 'thermal_energy'
            name_format = 'reserve_votemper_(.*)\.nc'            
            if use_total_salt_amount:
                title_text = "Total heat energy (J)"
                trend_unit = "J/decade"
            else: #calculate average salinity
                title_text = "Average temperature over model area (°C)"
                trend_unit = "°C/decade"
                data_multiplier = 1./model_area
                data_zero_point = 273.15 # gludge to deal with Celcius/Kelvin
    #    data_dir ='D:\\Data\\SmartSeaModeling\\'
        data_dir = sm.root_data_in+'derived_data\\figure_data_new\\'
        files = os.listdir(data_dir)
        dat={}
        extra_shift = -extra_shift_step*3.0  # used to shift whisker plots a bit
        for f in files:
            skip_this = True
            set_name=re.search(name_format,f)
            if(set_name):
                skip_this = False
                set_name = set_name.groups()[0]
                if(set_name == "REANALYSIS"):
                    set_name = "hindcast"
                    if(drop_hindcast):
                        skip_this = True
            if(not skip_this):
    #        dat[set_name]=pd.read_csv(data_dir+f,\
    #                             parse_dates=[0])
                D = xr.open_dataset(data_dir+f)
#                print(set_name)
    #            D = Dataset(data_dir+f)
                values =  np.array(D[variable])
                values = np.sum(values,1)
                values =  values*data_multiplier
                values = values - data_zero_point
                times =  np.array(D['time'])
                dat[set_name] = pd.DataFrame(list(zip(times,values)),\
                                   columns=['time',variable])
    #            times = D['time']
    #            times = netCDF4.num2date(times[:],times.units)
    #            dat[set_name] = pd.DataFrame({'time':times, 'value':values})
    #            dat[set_name]=pd.read_csv(data_dir+f,\
    #                                 parse_dates=[0])
                dat[set_name] = dat[set_name].set_index('time')
                D.close()
                
        plt.figure(figsize=fig_size)
        plt.title(title_text)
        #calculate the means for History, RCP4.5 and RCP8.5
        if(plot_combinations):
            dat["Control"] = pd.concat([dat['A001'],dat['B001'],dat['D001']])
            dat["RCP45"] = pd.concat([dat['A002'],dat['B002'],dat['D002']])
            dat["RCP85"] = pd.concat([dat['A005'],dat['B005'],dat['D005']])
            if(not plot_single_models): # remove the A,B,D thingies from the list
                for i in list(dat.keys()):
                    if(i.startswith('A') or i.startswith('B') or i.startswith('D')):
                        dat.pop(i)
        for s in dat:
            d=dat[s]
            d = d[(d.index>period['min']) & (d.index<period['max'])]
            if(plot_original):
                plt.plot(d.index,d[variable], label='_nolegend_', zorder=11,**sm.set_style(s,0.2))
        for s in dat:
            d=dat[s]
            d = d[(d.index>period['min']) & (d.index<period['max'])]
            if(len(d)>1):
                smooth_window = 12 #yearly 
                smoothed = d[variable].ewm(span = smooth_window,min_periods=smooth_window).mean()
                fitting_time = mp.dates.date2num(d.index)
                fitting = np.polyfit(fitting_time,d[variable],1)
                print("{} change: {:.3} {}".format(s,fitting[0]*3651.5, trend_unit))
                label_text = "{}:{:0.2} {}".format(s,fitting[0]*3651.5, trend_unit)
                if(plot_smoothed):
                    plt.plot(d.index,smoothed,label=label_text, zorder=15,**sm.set_style(s))
                    label_text = None # to prevent plotting the label more than once
                if(plot_trends):
    #                plt.plot(mp.dates.num2date(fitting_time),fitting[0]*fitting_time+fitting[1],label='_nolegend_', zorder=15,**set_style(s,0.4))
                    plt.plot(d.index,fitting[0]*fitting_time+fitting[1],label=label_text, zorder=15,**sm.set_style(s,0.4))
                    label_text = None # to prevent plotting the label more than once
                s_cloud = sm.set_style(s)
                s_cloud['alpha'] = 0.1
                s_cloud.pop('marker') # fill_betwen doesn't revognize marker, so this key must be ejected.
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
                    plt.plot(median.index+plot_shift,median[variable], label=label_text, zorder=16,**sm.set_style(s))
                    plt.fill_between(median.index+plot_shift,\
                                     mean[variable]-std[variable],\
                                     mean[variable]+std[variable],
                                     **s_cloud)
                    label_text = None # to prevent plotting the label more than once

                if(plot_yearly_mean):
                    mean_style = sm.set_style(s)
                    mean_style['alpha'] = 0.15
                    plt.plot(d_tmp.index,d_tmp[variable], label=label_text,zorder=16,**mean_style)
                    label_text = None # to prevent plotting the label more than once
                    
        plt.legend()
        plt.xlim([period['min'],period['max']])
        if(show_grid):
            plt.grid('on')
        extra = ""
        if(plot_combinations):
            extra+="comb"
        out_filename = "total_{}_{}-{}{}.png".format(\
                    variable, period['min'].year,period['max'].year, extra)
        plt.savefig(out_dir+ out_filename)
        print("saved figure: {} {}".format(out_dir, out_filename))
gathered_profile_trends = ValueSet()

#

#
# Plots conserning specific measurement points
#
#
#


if analyze_profiles:
#    variable = 'votemper'
#    variable = 'vosaline'
    yearly_means = {}
    full_point_data = {}
    for variable in profile_types:
    #    all_depths = [0.0,50.0, 100.0, 2000.0] #depth, if under the bottom, the lowest with number is accepted.
        all_depths = [0.0,'bottom_sample'] #depth, if under the bottom, the lowest with number is accepted.
    #    points = ['F64', 'SR5', 'MS4', 'C3', 'US5B', 'F16', 'BO3', 'F3', 'F9', 'BO5']
        points = ['F64', 'SR5', 'US5B', 'BO3']
        bottom_sample = {'F64':245.0, 'SR5':110.0, 'US5B':120.0, 'BO3':100.0}
        fixed_axis= None #[2.0,9.0] #None or [min, max]
        if(variable in ['vosaline']):
            variable_name = "Salinity"
        if(variable in ['votemper']):
            variable_name = "Temperature"
        full_point_data[variable] = {}
        for point in points:
            full_point_data[variable][point] = {}
            yearly_means[point] = {}
            for depth_in_list in all_depths:
                full_point_data[variable][point][depth_in_list] = {}
                data_dir = sm.root_data_in+'derived_data\\extracted_profiles\\'
                name_format = 'profile_{}_(.*)_{}.nc'.format(point,variable)
                files = os.listdir(data_dir)
                files = [i for i in files if re.match(name_format,i)]
                depth = 0.0  # default if no other defined
                dat={}
                if(depth_in_list == 'bottom_sample'):
                    depth_in = bottom_sample[point]
                else:
                    depth_in = depth_in_list
                for f in files:
                    set_name=re.search(name_format,f)
                    skip_this = True
                    if(set_name):
                        skip_this = False
                        set_name = set_name.groups()[0]
                        if(set_name == "REANALYSIS"):
                            set_name = "hindcast"
                            if(drop_hindcast):
                                skip_this = True
                    if(not skip_this):
#                        print(set_name)
                        if(not set_name in yearly_means[point].keys()):
                            yearly_means[point][set_name] = {}
                        D = Dataset(data_dir+f)
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
                        # upper gives cftime, convert to datetime
                        times = map(\
                                    lambda x: \
                                        dt.datetime.strptime(str(x),x.format),\
                                        times)
                        dat[set_name] = pd.DataFrame({'time':times,\
                                       variable:values[:,depth_layer],\
                                       'lat':lat,
                                       'lon':lon})
                        dat[set_name] = dat[set_name].set_index('time')
                        yearly_means[point][set_name][depth_in] = \
                            dat[set_name].groupby(pd.Grouper(freq='1AS')).mean()
                full_point_data[variable][point][depth_in_list] = dat.copy()
                #calculate the means for History, RCP4.5 and RCP8.5
                if(plot_combinations):
                    dat["Control"] = pd.concat([dat['A001'],dat['B001'],dat['D001']])
                    dat["RCP45"] = pd.concat([dat['A002'],dat['B002'],dat['D002']])
                    dat["RCP85"] = pd.concat([dat['A005'],dat['B005'],dat['D005']])
                    if(not plot_single_models): # remove the A,B,D thingies from the list
                        for i in list(dat.keys()):
                            if(i.startswith('A') or i.startswith('B') or i.startswith('D')):
                                dat.pop(i)

                
                extra_shift = -extra_shift_step*3.0  # used to shift whisker plots a bit
                plt.figure(figsize=fig_size)
                plt.title("{} on {} depth {:0.1f} m (Max Depth {:0.0f} m)"\
                              .format(variable_name, point,depth,max_depth))
                for s in dat:
                    d=dat[s]
                    d = d[(d.index>period['min']) & (d.index<period['max'])]
                    if(len(d)>0):
                        smooth_window = 12*3 #yearly 
                        smoothed = d[variable].ewm(span = smooth_window,\
                                        min_periods=smooth_window).mean()
                        fitting_time = mp.dates.date2num(d.index)
                        fitting = np.polyfit(fitting_time,d[variable],1)
                        print("{} change: {:.3} unit/year".format(s,fitting[0]*365.15))
                        label_text = "{}:{:0.3f} u/dec".format(s,fitting[0]*3651.5)
                        if(plot_original):
                            plt.plot(d.index,d[variable], label='_nolegend_',\
                                             zorder=11,**sm.set_style(s,0.2))
                        if(plot_smoothed):
                            plt.plot(d.index,smoothed,label=label_text, \
                                         zorder=15,**sm.set_style(s))
                            label_text = None
                        gathered_profile_trends.add(\
                                point,\
                                d['lat'].iloc[0],\
                                d['lon'].iloc[0],\
                                "{:0.1f}".format(depth),\
                                s,\
                                fitting[0]*365.15)
                        if(plot_trends):
                            plt.plot(d.index,\
                                     fitting[0]*fitting_time+fitting[1],\
                                     label=label_text, zorder=15,**sm.set_style(s,0.4))
                            label_text = None
                    ## Handle the yearly, decadal, etc.
                        s_cloud = sm.set_style(s)
                        s_cloud['alpha'] = 0.1
                        s_cloud.pop('marker') # fill_betwen doesn't revognize marker, so this key must be ejected.
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
                            plt.plot(median.index+plot_shift,median[variable], label=label_text, zorder=16,**sm.set_style(s))
                            plt.fill_between(median.index+plot_shift,\
                                             mean[variable]-std[variable],\
                                             mean[variable]+std[variable],
                                             **s_cloud)
                            label_text = None # to prevent plotting the label more than once
        
                        if(plot_yearly_mean):
                            mean_style = sm.set_style(s)
                            mean_style['alpha'] = 0.15
                            plt.plot(d_tmp.index,d_tmp[variable], label=label_text,zorder=16,**mean_style)
                            label_text = None # to prevent plotting the label more than once
                        
                        
                plt.legend()
                if(fixed_axis):
                    plt.ylim(fixed_axis[0],fixed_axis[1])
                plt.xlim([period['min'],period['max']])
                if(show_grid):
                    plt.grid('on')
#                print("saving",depth,point)
                if(depth_in_list == "bottom_sample"):
                    depth_str = "bottom"
                else:
                    depth_str = "{:.1f}m".format(depth)
                extra = ""
                if(plot_combinations):
                    extra += "comb"
                out_filename = "{}_profile_{}_{}_{}-{}{}.png".format(\
                            variable_name,\
                            point,\
                            depth_str,\
                            period['min'].year,\
                            period['max'].year, extra)
                
                plt.savefig(out_dir+"Profiles\\"+ out_filename)
                print("Saved: {} {}".format(out_dir+"Profiles\\",out_filename))
        #write trend analysis
        trend_file_name = \
            out_dir+'point_trends_{}.csv'.format(variable_name.lower())
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
        print(pd.read_csv(trend_file_name,'\t')\
              .to_latex(caption = variable_name, index = False))

#
#
#The trend plots 
#
#
#

if analyze_salt_trends:
    #open just saved file as pandas, and do some plotting
    scenarios = ['HISTORY','RCP4.5','RCP8.5']
    for scenario in scenarios:
        shade_color = 'b'
        if(variable_name == "Temperature"):
            shade_color = 'r'
        depths = [1.5, 50.0, 100.0]
        depth_vars = [0.5, 10.0, 15.0]
        for depth, depth_var in zip(depths,depth_vars):
            dataf = pd.read_csv(trend_file_name,sep='\t')
            d = dataf[dataf['scenario'] == scenario]
            d = d[d['depth'] > depth - depth_var]
            d = d[d['depth'] < depth + depth_var]
            d = d.sort_values('lat')
            figure = plt.figure(figsize=fig_size)
            plt.title("{} trend {} depth {:0.1f} m".format(\
                      variable_name, scenario,depth))
            plt.plot(d['lat'],d['mean'],'b*')
            plt.plot(d['lat'],[0]*len(d['lat']),'k',alpha=0.3)
            axis = figure.axes[0]
            axis.fill_between(d['lat'],d['max'],d['min'],\
                              facecolor = shade_color, alpha=0.2)
            for point,lat,val in zip(d['Point'],d['lat'],d['mean']):
                plt.text(lat,val,point)
            #plt.ylim(-0.02,0.04)
            out_filename = "{}_trends_{}_{:0.1f}m.png".format(\
                        variable_name, scenario, depth)
            plt.savefig(out_dir+out_filename)
            print("saved: {} {}".format(out_dir, out_filename))

if analyze_correlations:
    # correlate the river inflows
    inflow_numbers = []
    data_dir = sm.root_data_in + '\\derived_data\\inflow\\'
    files = os.listdir(data_dir)
    files = [x for x in files if x.endswith('csv')]
    inflow_dat={}
    for f in files:
        set_name=re.search('_([^_]*)\.csv',f).groups()[0]
        inflow_dat[set_name]=pd.read_csv(data_dir+f,\
                             parse_dates=[0])
        inflow_dat[set_name]['inflow'] = inflow_dat[set_name]['inflow']*\
                                1000000\
                                *60*60*24*365\
                                *0.0001*0.0001*0.0001  
                                #fixes one eror in csv creations, then
                                # changes unit from kg per second
                                # into km^3/year
        inflow_dat[set_name]=inflow_dat[set_name].set_index('time')
        multiplier=1.0
        if set_name == 'hindcast':
            multiplier = 30.5
        print(set_name,inflow_dat[set_name]['inflow'].sum()*multiplier)
    #calculate the means for History, RCP4.5 and RCP8.5
    if(plot_combinations):
        inflow_dat["Control"] = pd.concat([inflow_dat['A001'],inflow_dat['B001'],inflow_dat['D001']])
        inflow_dat["Control"].sort_index(inplace = True)
        inflow_dat["RCP45"] = pd.concat([inflow_dat['A002'],inflow_dat['B002'],inflow_dat['D002']])
        inflow_dat["RCP45"].sort_index(inplace = True)
        inflow_dat["RCP85"] = pd.concat([inflow_dat['A005'],inflow_dat['B005'],inflow_dat['D005']])
        inflow_dat["RCP85"].sort_index(inplace = True)
                    
    # calculate correlations
    correlation_set = '5meter'
    for correlation_set in list(boundary_data.keys()) + ['inflow']:
        print("####{}####".format(correlation_set))
        if(correlation_set == 'inflow'):
            corr_set = inflow_dat
        else:
            corr_set = boundary_data[correlation_set]
        for variable in full_point_data.keys():
            if(correlation_set == 'inflow'):
                variable2 = 'inflow'
                correlation_type = 'river'
            else:
                variable2 = variable
                correlation_type = 'boundary'
                
            for point in full_point_data[variable].keys():
                for depth in full_point_data[variable][point].keys():
                    correlation_values = []
                    dat = full_point_data[variable][point][depth]
                    print("=={},{},{}==".format(variable, point, depth))
                    for serie in dat.keys():
                        dat[serie] = dat[serie].sort_index()
                        dat[serie] = dat[serie][dat[serie].index>=period['min']] #to trim some 50's values off first.
                        max_lim = pd.DatetimeIndex([dat[serie].index.max(),\
                                                    corr_set[serie].index.max(),
                                                    period['max']]).min()
                        min_lim = pd.DatetimeIndex([dat[serie].index.min(),\
                                                    corr_set[serie].index.min(),
                                                    period['min']]).max()
                        #print(dat[serie].corr(boundary_data['5meter'][serie]))
                        dat[serie] = dat[serie][(dat[serie].index>=min_lim) \
                                               & (dat[serie].index<=max_lim)]
                        # dat[serie]['vosaline'].plot()
                        corr_set[serie] = \
                            corr_set[serie][(corr_set[serie].index>=min_lim) \
                                              & (corr_set[serie].index<=max_lim)]
                        # make sure both sets have the minimum value to get the bins right
                        if(not min_lim in dat[serie]):
                            dat[serie] = dat[serie].append(pd.DataFrame(None,[min_lim]))
                            dat[serie] = dat[serie].sort_index()
                        if(not min_lim in corr_set[serie]):
                            corr_set[serie] = corr_set[serie].append(pd.DataFrame(None,[min_lim]))
                            corr_set[serie] = corr_set[serie].sort_index()
                        # make sure both sets have the maximum value to get the bins right
                        if(not max_lim in dat[serie]):
                            dat[serie] = dat[serie].append(pd.DataFrame(None,[max_lim]))
                            dat[serie] = dat[serie].sort_index()
                        if(not max_lim in corr_set[serie]):
                            corr_set[serie] = corr_set[serie].append(pd.DataFrame(None,[max_lim]))
                            corr_set[serie] = corr_set[serie].sort_index()
                            
                        #corr_set[serie].plot()        
                        d1 = pd.DataFrame(dat[serie][variable])
                        d2 = corr_set[serie]
                        # let's take monthly means for the comparison
                        d1 = d1.groupby(pd.Grouper(freq='12M', offset = min_lim - d1.index.min())).mean()
                        d2 = d2.groupby(pd.Grouper(freq='12M', offset = min_lim - d2.index.min())).mean()
                        the_correlation = d1[variable].corr(d2[variable2])
                        correlation_values.append(the_correlation)
                        print("Correlation with {} {} in {} is {}".format(\
                                                    correlation_type, correlation_set,\
                                                    serie,\
                                                    the_correlation))        
                        plt.figure()
                        plt.plot(np.array(d1),np.array(d2[variable2]),'.',\
                                 label = "{:.4f}".format(the_correlation))
                        plt.legend()
                        plt.title("{},{},{}, {} {}\n{}".format(\
                                            variable,point,depth,
                                            correlation_type, correlation_set,\
                                            serie))
                        out_filename = "Correlation_{}_{}_{}_{}_{}_{}.png".format(\
                                            variable,point,depth,\
                                            correlation_type, correlation_set,\
                                            serie)
                        out_dir_plus = "\\{}\\".format(correlation_set)
                        if(not os.path.exists(out_dir+out_dir_plus)):
                            os.makedirs(out_dir+out_dir_plus)
                        plt.savefig(out_dir+out_dir_plus+out_filename)
                        print("saved: {} {}".format(out_dir+out_dir_plus, out_filename))
                        plt.close()
                    print("On average: {}\n\n".format(np.array(correlation_values).mean()))
