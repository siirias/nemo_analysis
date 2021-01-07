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
#from mpl_toolkits.basemap import Basemap
#from mpl_toolkits.basemap import Basemap, shiftgrid, cm
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from netCDF4 import Dataset
from smartseahelper import smh
import os
import cmocean
import re


font_size=10.0
resolution='i'  #h
projection='laea'
data_dir= "D:\\Data\\SmartSeaModeling\\SharkExamples\\"
out_dir = "D:\\Data\\SmartSeaModeling\\Images\\"
plot_area = [17.0, 26.0, 60.0, 66.0]
figure_size = (10,10)
the_proj = ccrs.PlateCarree()

def create_main_map(the_proj):
        fig=plt.figure(figsize=figure_size)
        plt.clf()
        ax = plt.axes(projection=the_proj)
        ax.set_extent(plot_area)
        ax.set_aspect('auto')
        ax.coastlines('10m',zorder=4, alpha = 0.5)
        ax.add_feature(cfeature.NaturalEarthFeature('physical', 'land', '10m',\
                                                edgecolor='face', \
                                                facecolor='#555560', alpha = 0.3))
        gl = ax.gridlines(crs=the_proj, draw_labels=True,
                  linewidth=2, color='gray', alpha=0.1, linestyle='-')
        gl.xlabels_top = False
        gl.ylabels_right = False
        return fig

files = os.listdir(data_dir)
dat={}
name_format = '\d*_shark(.*)\.nc'
for f in files:
    point_name=re.search(name_format,f)
    if(point_name):
        point_name = point_name.groups()[0]
        print(point_name)
        with Dataset(data_dir+f) as D:
            lat = D['nav_lat'][0].data[0]
            lon = D['nav_lon'][0].data[0]
            dat[point_name] = {'lat':lat, 'lon':lon}

##Let's try to plot soemthing to start with:
#lon_min=16;lat_min=60;lon_max=26.01;lat_max=66.01;
#first setup teh main map:

#fig=plt.figure(figsize=(10,10))
#if projection in ['laea']:
#    lat_0=0.5*(lat_min+lat_max)
#    lon_0=0.5*(lon_min+lon_max)
#    bmap = Basemap(llcrnrlon=lon_min,llcrnrlat=lat_min,urcrnrlon=lon_max,urcrnrlat=lat_max, \
#                   lat_0=lat_0, lon_0=lon_0,resolution = resolution,
#                   projection=projection,fix_aspect=False)  #resolution c,l,i,h,f
#elif projection in ['merc','cyl']:
#    bmap = Basemap(llcrnrlon=lon_min,llcrnrlat=lat_min,urcrnrlon=lon_max,urcrnrlat=lat_max, \
#                   resolution = resolution,
#                   projection=projection,fix_aspect=False)  #resolution c,l,i,h,f
#    
#bmap.drawcoastlines(zorder=21,linewidth=0.5,color='gray')
#bmap.fillcontinents([0.9,0.9,0.9],lake_color=[0.85,0.85,0.85],zorder=20) 
#bmap.drawparallels(np.arange(lat_min,lat_max,1.),linewidth=1,zorder=50,labels=[True,False,False,False],dashes=[1,0],color="#00000020",fontsize=10)        
#bmap.drawmeridians(np.arange(lon_min,lon_max,2.),linewidth=1,zorder=50,labels=[False,False,False,True],dashes=[1,0],color="#00000020",fontsize=10)        
fig = create_main_map(the_proj)
for name in dat:
    lon = dat[name]['lon']
    lat = dat[name]['lat']
#    mlon,mlat = bmap(lon,lat)
#    plt.plot(mlon,mlat,'r.',zorder = 100)
    plt.plot(lon,lat,'r.',zorder = 100, transform = the_proj)
    if(name == 'B7'):
        lon += 0.2
    if(name == 'NB1'):
        lon -= 0.2
#    txtlon,txtlat = bmap(lon,lat+0.07)
    plt.text(lon, lat + 0.07, name, \
             fontsize = 'medium', horizontalalignment = 'center',\
             zorder = 120, transform = the_proj)
#    txtlon,txtlat = bmap(lon,lat-0.20)
    plt.text(lon, lat - 0.2,\
             u"{:.4}° E\n {:.4}° N".format(dat[name]['lon'],dat[name]['lat']), \
             horizontalalignment = 'center', fontsize = 'x-small', \
             zorder = 120, transform = the_proj)
plt.title('SHARK Points stored in SmartSea')                   
plt.savefig(out_dir+"Shark_points.png",facecolor='w',dpi=300)
        