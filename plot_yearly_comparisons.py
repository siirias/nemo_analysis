# -*- coding: utf-8 -*-
"""
Created on Mon Jul 16 15:52:51 2018

@author: siirias
"""
import datetime as dt
import calendar
import matplotlib as mp
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from smartseahelper import smh
from matplotlib import rc
from netCDF4 import Dataset

ss=smh()
datadir=ss.root_data_out+"/tmp/Images_for_video/"
ss.main_data_folder= ss.root_data_out+"/derived_data/"
datadir = ss.root_data_out+"/tmp/Images_for_video/" #where everyt output is stored
name_markers=['A001','A002','B001','B002']
name_markers=['A002']
other_name_markers=['A001','B001','C001']
plot_types=['normal','long','circular']
plot_types=['long']
substract_climatology=True
#plot_types='circular' #normal, circular, long
compare_two=False
y_label=""
animation=False
plot_param='mean_SST'#'mean_SST' #ice, mean_SST, mean_SSS
climat_param='SST'
name_markers=['A001','A002']
for plot_type in plot_types:
    for name_marker in name_markers:
        if(substract_climatology):
            climatology_file=ss.root_data_out+'derived_data/climatology_{}_d_{}.nc'.format(name_marker,climat_param)  #this one tells which dataseries is handled.
            climatology_data_nc=Dataset(climatology_file)
            climatology_data=climatology_data_nc.variables['{}_mean'.format(climat_param)][:]
            clim_lats=climatology_data_nc.variables['nav_lat'][:]
            clim_lons=climatology_data_nc.variables['nav_lon'][:]
            clim_lats,clim_lons=ss.fix_latslons(clim_lats,clim_lons)
            areas=ss.give_areas(clim_lats,clim_lons)
            climat_mean_variable=np.zeros(climatology_data.shape[0])
            for t_f in range(climatology_data.shape[0]):
                variable_times_area=np.sum(climatology_data[t_f,:,:]*areas)
                if(np.sum(~climatology_data[t_f,:,:].mask)>0): #which means there is some actual number here
                    climat_mean_variable[t_f]=variable_times_area/np.sum(~climatology_data[t_f,:,:].mask*areas)            
                else:
                    climat_mean_variable[t_f]=0.0
            climatology_data_nc.close()
        if(not compare_two):
            other_name_markers=['empty']
        for other_name_marker in [x for x in other_name_markers if x is not name_marker]:
            
            
            if plot_param=='ice':
                x_name='time'
                data_name='ice_extent'
                in_filename='ice_extent_{}.csv'.format(name_marker)
                y_label=r"ice surface ($km^2$)"
                file_identifier=data_name
                val_min=0.
                grid_min=val_min
                if(plot_type=="circular"):
                    val_min=-60000.
                    grid_min=0
                    
                val_max=120000.
                val_step=10000.
                first_day=0
                if(compare_two):
                    other_in_filename='ice_extent_{}.csv'.format(other_name_marker)
                    val_min=-100000.
                    val_max=100000.
                    val_step=10000.
                    grid_min=-100000  #stop grid lines
                    y_label=r"ice surface difference ($km^2$)"
                
    
            if plot_param=='mean_SST':
                x_name='time'
                data_name='mean_SST'
                in_filename='mean_SST_{}.csv'.format(name_marker)
                file_identifier=data_name
                y_label=r"mean SST ($^o$C)"
                val_min=-20.
                val_max=20.
                val_step=2.
                grid_min=-2  #stop grid lines
                first_day=0
                x_step=30
                if(compare_two):
                    other_in_filename='mean_SST_{}.csv'.format(other_name_marker)
                    val_min=-8.
                    val_max=8.
                    val_step=1.
                    grid_min=-8  #stop grid lines
                    y_label=r"mean SST difference ($^o$C)"
            
            if plot_param=='mean_SSS':
                x_name='time'
                data_name='mean_SSS'
                in_filename='mean_{}_{}.csv'.format('SSS',name_marker)
                file_identifier=data_name
                y_label=r"mean SSS "
                val_min=3.8
                val_max=5.6
                val_step=0.2
                grid_min=3  #stop grid lines
                first_day=0
                x_step=30
                if(compare_two):
                    other_in_filename='mean_{}_{}.csv'.format('SSS',other_name_marker)
                    val_min=-8.
                    val_max=8.
                    val_step=1.
                    grid_min=0  #stop grid lines
                    y_label=r"mean SSS difference"
            
            data=pd.read_csv(ss.main_data_folder+in_filename,parse_dates=[x_name])
            if(compare_two):
                other_data=pd.read_csv(ss.main_data_folder+other_in_filename,parse_dates=[x_name])
                file_identifier+="_comp"
            time_data=list(data['time'])
            clim_fix=np.zeros(len(time_data))
            if(substract_climatology):
                clim_fix=map(lambda x:climat_mean_variable[x.timetuple().tm_yday-1],time_data) 
            fig=plt.figure(figsize=(6,6))
            plt.clf();
            
            if plot_type=='circular':
                rotation=0./360.*2.*np.pi
                x_axis=list(map(lambda x:2.*np.pi*x.dayofyear/366.+ rotation,list(data[x_name])))
                ax=plt.subplot(111,projection='polar');
                ax.set_ylim(val_min,val_max)
                ax.set_yticks(np.arange(grid_min,val_max,val_step))
                 #ax.set_rmax(val_max)
            #    ax.set_rmin(val_min)
                #ax.set_rticks(np.arange(val_min,val_max,val_step))  # less radial ticks
                ax.set_theta_zero_location("N")
                ax.set_theta_direction(-1)
                ax.set_rlabel_position(180.0)  # get radial labels away from plotted line
                ax.set_xticks(np.arange(0,np.pi*2.,np.pi*2./12.))    
                ax.set_xticklabels(map(lambda x:x[:3],calendar.month_name[1:]))
                y_label_distance=40
                x_label=''
                print( "CIRC")
            if plot_type=='normal':
                x_axis=list(map(lambda x:(x.dayofyear+first_day)%366.,list(data[x_name])))
            #    x_axis=data[x_name]
                ax=plt.subplot(111);
                ax.set_ylim(grid_min,val_max)
                ax.set_yticks(np.arange(grid_min,val_max,val_step))
                ax.set_xlim(0,365)
                ax.set_xticks(np.arange(0,365,x_step))
                y_label_distance=10
                x_label='day of year'
                print ("NORM")
            if plot_type=='long':
                x_axis=data[x_name]
            #    x_axis=data[x_name]
                ax=plt.subplot(111);
                ax.set_ylim(grid_min,val_max)
                ax.set_yticks(np.arange(grid_min,val_max,val_step))
                y_label_distance=10
                x_label='time'
                print ("long")
            
            segments=len(x_axis)-1
            at_once=5
            #segments=20
            if(compare_two):
                    tmp_x=np.arange(0,367,1)
                    if(plot_type in ['circular']):
                        tmp_x=np.arange(0,2.*np.pi,0.001)
                    plt.plot(tmp_x,np.zeros((len(tmp_x),)),color='#000000',zorder=10,alpha=0.6,linewidth=2)
            
            time_string="{} -- {}".format(data[x_name][0],data[x_name][len(data[x_name])-1])
            plt.title("{}, {}\n{} ".format(name_marker,plot_type,time_string))
            if(compare_two):
                plt.title("difference {}-{}, {}\n{} ".format(name_marker,other_name_marker,data_name,time_string))
            plt.xlabel(x_label)
            plt.ylabel(y_label,labelpad=y_label_distance)
            
            running_number=0
            for i in range(0,segments,at_once):
                color=(float(i)/segments,0.,1.-float(i)/segments)
                zorder=10+300*(float(i)/segments)
                skip_this_segment=False
                if(plot_type in ['normal'] and np.diff(x_axis[i:i+at_once+1]).min()<0.0):  #skip drawing the ones which overlap border
                    skip_this_segment=True
                if(not skip_this_segment):  #skip drawing the ones which overlap border
                    if(compare_two):
                        plt.plot(x_axis[i:i+at_once+1],data[data_name][i:i+at_once+1]-other_data[data_name][i:i+at_once+1],color=color,zorder=zorder,alpha=0.3)
                    else:
                        plt.plot(x_axis[i:i+at_once+1],data[data_name][i:i+at_once+1]-clim_fix[i:i+at_once+1],color=color,zorder=zorder,alpha=0.3) #0.3
                if(plot_type in ['normal','long']):
                    plt.grid(True,alpha=0.4)
                if animation:
                    annotation=plt.annotate(time_data[i].strftime("%Y-%m-%d"),xy=(0.25, 0.95), xycoords='axes fraction',zorder=100)
                    plt.savefig("{}{}{:05d}.png".format(datadir,file_identifier,running_number),facecolor='w',dpi=300)
                    annotation.remove()
                running_number+=1
            print(time_string)
            if(not animation):
                    running_number=1
                    names_str=name_marker
                    if(compare_two):
                        names_str="{}vs{}".format(name_marker,other_name_marker)
                        
                    plt.savefig("{}{}{}_{}.png".format(datadir,file_identifier,names_str,plot_type),facecolor='w',dpi=300)
                
