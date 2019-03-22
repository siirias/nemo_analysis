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
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from netCDF4 import Dataset
from smartseahelper import smh
import os
import cmocean

startdate=datetime.datetime(1975,1,1)
enddate=datetime.datetime(2005,12,31)
#startdate=datetime.datetime(2006,1,1)
#enddate=datetime.datetime(2058,12,31)
datadir = ".//Images_for_video//"
series_name='Hindcast_A001_SSS'
ss=smh()
ss.grid_type='T'
ss.interwall='d'
ss.main_data_folder= "d:\\OUTPUTA001\\"

#just_one=False
just_one=True

filenames=ss.filenames_between(startdate,enddate)
ok_files=0
files_working=[]
for f in filenames:
    if(os.path.isfile(ss.main_data_folder+f)):
        ok_files+=1
        files_working.append(f)
    else:
        print f
print
print "ok {} out of {}".format(ok_files,len(filenames))
    
#Let's try to plot soemthing to start with:
lon_min=16;lat_min=60;lon_max=26.01;lat_max=66.01;
running_number=0
var1 = 'SST'  #'SST', 'SSS'
var_min=1.0 #-0.5
var_max=6.0 #20.0
var1_cm=cmocean.cm.haline #cmocean.cm.thermal

var2 = 'icecon'  #None, 'icecon'
var2_min=0.0
var2_max=1.0
var2_cm='gray'
font_size=10.0
resolution='h'
projection='laea'
if(just_one):
    files_working=[files_working[0]]
    
#first setup teh main map:
fig=plt.figure(figsize=(10,6))
if projection in ['laea']:
    lat_0=0.5*(lat_min+lat_max)
    lon_0=0.5*(lon_min+lon_max)
    projection=ccrs.LambertAzimuthalEqualArea(lat_0,lon_0)
#    projection=ccrs.PlateCarree()
    cp_ax = plt.axes(projection=projection)
    #Basemap(llcrnrlon=lon_min,llcrnrlat=lat_min,urcrnrlon=lon_max,urcrnrlat=lat_max, \
   #                lat_0=lat_0, lon_0=lon_0,resolution = resolution,projection=projection,fix_aspect=False)  #resolution c,l,i,h,f
elif projection in ['merc','cyl']:
    projection=ccrs.Mercator()
    cp_ax = plt.axes(projection=projection)
#    bmap = Basemap(llcrnrlon=lon_min,llcrnrlat=lat_min,urcrnrlon=lon_max,urcrnrlat=lat_max, \
#                   resolution = resolution,projection=projection,fix_aspect=False)  #resolution c,l,i,h,f
    

#bmap.drawcoastlines(zorder=21,linewidth=0.5,color='gray')
#cp_ax.coastlines()
cp_ax.add_feature(cfeature.LAND)
##bmap.fillcontinents([0.9,0.9,0.9],lake_color=[0.85,0.85,0.85],zorder=20) 

#BM# bmap.drawparallels(np.arange(lat_min,lat_max,1.),linewidth=1,zorder=50,labels=[True,False,False,False],dashes=[1,0],color="#00000020",fontsize=10)        

#BM# #bmap.drawmeridians(np.arange(lon_min,lon_max,2.),linewidth=1,zorder=50,labels=[False,False,False,True],dashes=[1,0],color="#00000020",fontsize=10)        

is_first=True
for f in files_working:
    data=Dataset(ss.main_data_folder+f)
    d=data.variables[var1][:]
    d=np.ma.masked_where(d==0.0,d)
    if(var2 is not None):
        ice_d=data.variables[var2][:]
        ice_d=np.ma.masked_where(ice_d<0.2,ice_d)
    lons = data.variables['nav_lon'][:]
    lats = data.variables['nav_lat'][:]
    lats,lons=ss.fix_latslons(lats,lons)
    times=data.variables['time_counter'][:].data
    if(just_one):
        times=[times[0]]
    for time_frame in range(len(times)):
        time=ss.nemo_time_to_datetime(times[time_frame])
        #cb=plt.colorbar(fraction=0.027, pad=0.01)
        colors_fig=plt.pcolormesh(lons,lats,d[time_frame,:,:],vmin=var_min,vmax=var_max,zorder=0,cmap=var1_cm,transform=projection)
        if is_first:
            cb=plt.colorbar()
            cb.set_clim(vmin=var_min,vmax=var_max)
            cb.ax.tick_params(labelsize=font_size)
            is_first=False
        if(var2 is not None):
            ice_fig=plt.pcolormesh(lons,lats,ice_d[time_frame,:,:],vmin=var2_min,vmax=var2_max,zorder=18,cmap=var2_cm,transform=projection)

        cont_fig=plt.contour(lons,lats,d[time_frame,:,:],vmin=var_min,vmax=var_max,zorder=15,colors='black',linewidths=0.5,alpha=0.5,transform=projection)
        cont_labels=plt.clabel(cont_fig,inline=1,fontsize=7,fmt="%1.1f",zorder=15,transform=projection)
        annotation=plt.annotate(time.strftime("%Y-%m-%d"),xy=(0.25, 0.95), xycoords='axes fraction',zorder=100)
        plt.savefig("{}{}{:05d}.png".format(datadir,series_name,running_number),facecolor='w',dpi=300)
        if(not just_one):
            #clean up the changing things, so we don't have to do everything again, just these:
            colors_fig.remove()
            if(var2 is not None):
                ice_fig.remove()
            annotation.remove()
            for i in cont_labels:
                i.remove()
            for i in cont_fig.collections:
                i.remove()
#           plt.close(fig)
        running_number+=1
        print "Image {} of (approx){} done.".format(running_number,len(files_working)*30)
        