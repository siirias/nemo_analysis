# -*- coding: utf-8 -*-
"""
Created on Tue Jul  3 14:11:45 2018

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
compare_two=False
other_name_marker='B002'
var1 = 'SST'  #'VEL', 'SST', 'SSS'

ss=smh()
ss.grid_type='T'
ss.interwall='d'
ss.main_data_folder= ss.root_data_in+"/OUTPUT{}/".format(name_marker)
datadir = ss.root_data_out+"/tmp/Images_for_video2/" #where everyt output is stored

climatology_file=ss.root_data_out+'derived_data/climatology_{}_d_{}.nc'.format(name_marker,var1)  #this one tells which dataseries is handled.
series_name='climatology{}'.format(name_marker)

if(compare_two):
    other_climatology_file=ss.root_data_out+'derived_data/climatology_{}_d_{}.nc'.format(other_name_marker,var1)  #this one tells which dataseries is handled.
    series_name='climatology{}vs{}'.format(name_marker,other_name_marker)
    


#Let's try to plot soemthing to start with:
lon_min=16;lat_min=60;lon_max=26.01;lat_max=66.01;
running_number=0
show_contours=True
plotted_depth=0. #meters

var2 = None  #None, 'icecon'
var2_min=0.0
var2_max=1.0
var2_cm='gray'


if var1 in ['VEL']:
    var_min=0.0
    var_max=0.15
    var1_cm=cmocean.cm.speed
    var1_cm=cmocean.cm.speed
    ss.grid_type='U'  #U,V,T
    ss.interval='m'
    plotted_depth=-10000. #meters
    var2=None
    show_contours=False
if var1 in ['SST']:
    var_min=-0.5
    var_max=20.0
    var1_cm=cmocean.cm.thermal
    ss.grid_type='T'  #U,V,T
    ss.interval='d'
    if(compare_two):
        var_min=-5.
        var_max=5.
if var1 in ['SSH']:
    var_min=1.0 
    var_max=6.0 
    var1_cm=cmocean.cm.haline 
    ss.grid_type='T'  #U,V,T
    ss.interval='d'



font_size=10.0
resolution='h'
projection='laea'
just_one=False

    


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
data=Dataset(climatology_file)
if(compare_two):
    other_data=Dataset(other_climatology_file)
if(var1 in ['SST','SSH','SSS']):
    d=data.variables['{}_mean'.format(var1)][:]
    
    if(compare_two):
        other_d=other_data.variables['{}_mean'.format(var1)][:]
        d=d-other_d
#    d=np.ma.masked_where(d==0.0,d)
if(var1 in ['VEL']):
    data2=Dataset(ss.main_data_folder+re.sub("_grid_.","_grid_V",f))
    if(plotted_depth>=0):
        depth_index=np.abs(data.variables['depthu'][:]-plotted_depth).argmin()
        d=data.variables['uos'][:,depth_index,:,:]
        d2=data2.variables['vos'][:,depth_index,:,:]
    else: #negative plotted depth measn bottom.
        d=ss.give_bottom_values(data.variables['uos'][:,:,:,:])
        d2=data2.variables['vos'][:,:,:]
    d=np.ma.masked_where(d==0.0,d)
    #this one requires that we take another file too:
    d2=np.ma.masked_where(d==0.0,d)
    d=np.sqrt(d*d+d2*d2)
    data2.close()
    
if(var2 is not None):
    ice_d=data.variables[var2][:]
    ice_d=np.ma.masked_where(ice_d<0.2,ice_d)
lons = data.variables['nav_lon'][:]
lats = data.variables['nav_lat'][:]
lats,lons=ss.fix_latslons(lats,lons)
times=data.variables['day_of_year'][:].data
data.close()
if compare_two:
    other_data.close()
if(just_one):
    times=[times[0]]
for time_frame in range(len(times)):
    time=time_frame
    tmp_lon,tmp_lat=bmap(lons,lats)
    #cb=plt.colorbar(fraction=0.027, pad=0.01)
    colors_fig=bmap.pcolormesh(tmp_lon,tmp_lat,d[time_frame,:,:],vmin=var_min,vmax=var_max,zorder=-3,cmap=var1_cm)
    if is_first:
        cb=plt.colorbar()
        cb.set_clim(vmin=var_min,vmax=var_max)
        cb.ax.tick_params(labelsize=font_size)
        is_first=False
    if(var2 is not None):
        ice_fig=bmap.pcolormesh(tmp_lon,tmp_lat,ice_d[time_frame,:,:],vmin=var2_min,vmax=var2_max,zorder=16,cmap=var2_cm)

    if show_contours:
        cont_fig=bmap.contour(tmp_lon,tmp_lat,d[time_frame,:,:],vmin=var_min,vmax=var_max,zorder=15,colors='black',linewidths=0.5,alpha=0.5)
        cont_labels=plt.clabel(cont_fig,inline=1,fontsize=5,fmt="%1.1f")
    annotation=plt.annotate("Year day: {}".format(time_frame+1),xy=(0.25, 0.95), xycoords='axes fraction',zorder=100)
    plt.savefig("{}{}{:05d}.png".format(datadir,series_name,running_number),facecolor='w',dpi=300)
    if(not just_one):
        #clean up the changing things, so we don't have to do everything again, just these:
        colors_fig.remove()
        if(var2 is not None):
            ice_fig.remove()
        annotation.remove()
        if show_contours:
            for i in cont_labels:
                i.remove()
            for i in cont_fig.collections:
                i.remove()
#            plt.close(fig)
    running_number+=1
    print("Image {} done.".format(running_number))
    