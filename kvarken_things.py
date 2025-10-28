# -*- coding: utf-8 -*-
"""
Created on Mon Feb 22 17:32:35 2021

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
import datetime as dt

output_dir = "D:\\Data\\Figures\\SmartSea\\kvarken\\"
data_dir = "E:\\SmartSea\\kvarken\\"  
plot_animation = False
plot_analysis  = True
datasets = ['A001']
figure_size = (10,5)
fig_dpi = 300
mod_min_lat = 62.5
mod_max_lat = 64.0
mod_min_lon = 17.8
mod_max_lon = 24.0
mod_shape_lat = 90
mod_shape_lon = 204
lat = np.linspace(mod_min_lat,mod_max_lat,mod_shape_lat)  # This is a gludge, as the original lat lon ar bit weird.
lon = np.linspace(mod_min_lon,mod_max_lon,mod_shape_lon)
selected_points={
        'north':(21.780,63.595),
        'south':(20.612,63.448)
        }

def moving_average(a, n=3):
    ret = np.cumsum(a, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    ret = ret[n - 1:] / n
    ret = np.concatenate((np.tile(ret[0],int(np.floor(n/2.0))),
                          ret,
                          np.tile(ret[-1],int(np.floor(n/2.0))-(n+1)%2)))
    return ret

for dataset in datasets:
    file = 'kvarken_{}.nc'.format(dataset)
    data = xr.open_dataset(data_dir+file)
    time = pd.to_datetime(np.array(data['time_counter']))
    data = np.array(data['SBS'])  # Take the bottom salinity
    data[data == 0.0] = np.nan # remove the ones marked as natural zeros away
    if(plot_analysis):
        plt.figure()
        for i in selected_points:
            lon_point = np.abs(lon-selected_points[i][0]).argmin()
            lat_point = np.abs(lat-selected_points[i][1]).argmin()
            mean_value = np.mean(data[:,lat_point,lon_point])
#            plt.plot(time,data[:,lat_point,lon_point]-mean_value,label=i)
#            plt.fill_between(time,data[:,lat_point,lon_point]-mean_value,alpha=0.5)
            d = data[:,lat_point,lon_point] - mean_value
            d = moving_average(d,7) # one week average
            perc = 80.0
            d_filt = d>np.percentile(d,perc)
            plt.plot(time,d,label=i)
            plt.fill_between(time,d,np.percentile(d,perc),
                             where = d>np.percentile(d,perc),alpha=0.5)
        plt.legend()
        plt.figure(figsize=figure_size)
        plt.pcolor(lon,lat,data[0],vmin=1.0,vmax=7.,cmap=cmo.cm.haline)
        for i in selected_points:
            plt.plot(selected_points[i][0],selected_points[i][1],\
                     'o',markersize=10,fillstyle='none')

        plt.title("{} {}".format(dataset,'sample'))
        plt.colorbar()
    if(plot_animation):
        for t in range(len(time)):
            plt.figure(figsize=figure_size)
            plt.pcolor(lon,lat,data[t],vmin=1.0,vmax=7.,cmap=cmo.cm.haline)
            for i in selected_points:
                plt.plot(selected_points[i][0],selected_points[i][1],\
                         'o',markersize=10,fillstyle='none')
            
            
            plt.title("{} {}".format(dataset,time[t].strftime("%Y-%m-%d")))
            plt.colorbar()
            filename = "kvark_{}_{:04}.png".format(dataset,t)
            plt.savefig(output_dir+filename,\
                            facecolor='w',dpi=fig_dpi,bbox_inches='tight')
    
            print("Saved: {}".format(output_dir+filename))
            
            plt.close()
#    flattened = np.nanmean(data,2)
#    plt.figure()
#    plt.pcolor(time,lat,flattened.transpose())
