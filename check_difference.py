#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  4 10:44:59 2018

@author: siirias
"""


import datetime
import matplotlib as mp
import matplotlib.pyplot as plt
import numpy as np
from scipy.io import netcdf
from mpl_toolkits.basemap import Basemap
from mpl_toolkits.basemap import Basemap, shiftgrid, cm
from netCDF4 import Dataset
from smartseahelper import smh
import os
import cmocean
import re

name_marker='A002'  #this one tells which dataseries is handled.
other_name_marker='B001'  #this one tells which dataseries is handled.
compare_two=False

startdate=datetime.datetime(1985,1,1)
enddate=datetime.datetime(1986,12,31)
ss=smh()
ss.grid_type='T'
ss.interwall='d'
ss.save_interval='year'
ss.main_data_folder= ss.root_data_in+"/tmp/ice_test/"
ss.file_name_format="SS-GOB_1{}_{}_{}_grid_{}.nc"
datadir = ss.root_data_out+"/tmp/ice_test/" #where everyt output is stored

filenames=ss.filenames_between(startdate,enddate)
ok_files=0
files_working=[]

#plotting
font_size=10.0
resolution='h'
projection='laea'
lon_min=16;lat_min=60;lon_max=26.01;lat_max=66.01;
#first setup teh main map:
fig=plt.figure(figsize=(6,6))
if projection in ['laea']:
    lat_0=0.5*(lat_min+lat_max)
    lon_0=0.5*(lon_min+lon_max)
    bmap = Basemap(llcrnrlon=lon_min,llcrnrlat=lat_min,urcrnrlon=lon_max,urcrnrlat=lat_max, \
                   lat_0=lat_0, lon_0=lon_0,resolution = resolution,
                   projection=projection,fix_aspect=False)  #resolution c,l,i,h,f
elif projection in ['merc','cyl']:
    bmap = Basemap(llcrnrlon=lon_min,llcrnrlat=lat_min,urcrnrlon=lon_max,urcrnrlat=lat_max, \
                   resolution = resolution,
                   projection=projection,fix_aspect=False)  #resolution c,l,i,h,f
    
bmap.drawcoastlines(zorder=21,linewidth=0.5,color='gray')
#bmap.fillcontinents([0.9,0.9,0.9],lake_color=[0.85,0.85,0.85],zorder=20) 
bmap.drawparallels(np.arange(lat_min,lat_max,1.),linewidth=1,zorder=50,labels=[True,False,False,False],dashes=[1,0],color="#00000020",fontsize=10)        
bmap.drawmeridians(np.arange(lon_min,lon_max,2.),linewidth=1,zorder=50,labels=[False,False,False,True],dashes=[1,0],color="#00000020",fontsize=10)        
is_first=True
for f in filenames:
    if(os.path.isfile(ss.main_data_folder+f)):
        ok_files+=1
        files_working.append(f)
    else:
        print(f)
    time_frame=0
   
    data=Dataset(ss.main_data_folder+f)
    d=data.variables['icevolume'][:]
    d=np.ma.masked_where(d==0.0,d)
    times=data.variables['time_counter'][:].data
    lons = data.variables['nav_lon'][:]
    lats = data.variables['nav_lat'][:]
    lats,lons=ss.fix_latslons(lats,lons)
    var_min=0.
    var_max=5.
    var1_cm=cmocean.cm.thermal
    
    time=ss.nemo_time_to_datetime(times[time_frame])
    tmp_lon,tmp_lat=bmap(lons,lats)
    colors_fig=bmap.pcolormesh(tmp_lon,tmp_lat,d[time_frame,:,:],vmin=var_min,vmax=var_max,zorder=-3,cmap=var1_cm)
    cb=plt.colorbar()
    cb.set_clim(vmin=var_min,vmax=var_max)
    cb.ax.tick_params(labelsize=font_size)