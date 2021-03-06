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


name_marker='A005'  #this one tells which dataseries is handled.
other_name_marker='A005_d_SST'  #this one tells which dataseries is handled.
font_size=10.0
resolution='h'
projection='laea'
just_one=False
climatology_dir='/derived_data//'
climatology_name='climatology_{}.nc'.format(other_name_marker)
if '1' in name_marker: #the 001 series are hindcasts, all other scenarios
    startdate=datetime.datetime(1975,1,1)
    enddate=datetime.datetime(2005,12,31)
else:
    startdate=datetime.datetime(2006,1,1)
    enddate=datetime.datetime(2058,12,31)
ss=smh()
ss.grid_type='T'
ss.interwall='d'
ss.main_data_folder= ss.root_data_in+"/OUTPUT{}/".format(name_marker)
datadir = ss.root_data_out+"/tmp/Images_for_video2/" #where everyt output is stored

image_per_month={'h':24*30.5,'d':30.5,'m':1}


var1 = 'SST'  #'VEL', 'SST', 'SSS'
clim_var1='SST_mean'
combined_name_markers=name_marker
if '1' in name_marker: #the 001 series are hindcasts, all other scenarios
    series_name='Hindcast_{}_{}'.format(combined_name_markers,var1)
else:
    series_name='Scenario_{}_{}'.format(combined_name_markers,var1)
    

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
    var_min=-10.
    var_max=10.0
    var1_cm=cmocean.cm.thermal
    ss.grid_type='T'  #U,V,T
    ss.interval='d'
if var1 in ['SSH']:
    var_min=1.0 
    var_max=6.0 
    var1_cm=cmocean.cm.haline 
    ss.grid_type='T'  #U,V,T
    ss.interval='d'
var1_cm=cmocean.cm.balance




filenames=ss.filenames_between(startdate,enddate)
ok_files=0
files_working=[]
for f in filenames:
    if(os.path.isfile(ss.main_data_folder+f)):
        ok_files+=1
        files_working.append(f)
    else:
        print(f)
print()
print("ok {} out of {}".format(ok_files,len(filenames)))

def give_bottom_values(array4d):
    bottom_index=(array4d[0,:,:,:]!=0.0).sum(axis=0)-1 #first axis time, then depth, x,y
    bottom_index[bottom_index<0]=0
    values=array4d[:,0,:,:].copy()
    for i in range(values.shape[1]):
        for j in range(values.shape[2]):
            if ~bottom_index.mask[i,j]:
                values[:,i,j]=array4d[:,bottom_index[i,j],i,j]
    return values
    



if(just_one):
    files_working=[files_working[0]]
    
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
    
climatology_data=Dataset(ss.root_data_in+climatology_dir+climatology_name)
for f in files_working:
    data=Dataset(ss.main_data_folder+f)
        
    if(var1 in ['SST','SSH']):
        d=data.variables[var1][:]
        d=np.ma.masked_where(d==0.0,d)
    if(var1 in ['VEL']):
        data2=Dataset(ss.main_data_folder+re.sub("_grid_.","_grid_V",f))
        if(plotted_depth>=0):
            depth_index=np.abs(data.variables['depthu'][:]-plotted_depth).argmin()
            d=data.variables['uos'][:,depth_index,:,:]
            d2=data2.variables['vos'][:,depth_index,:,:]
        else: #negative plotted depth measn bottom.
            d=give_bottom_values(data.variables['uos'][:,:,:,:])
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
    times=data.variables['time_counter'][:].data
    #TODO: Tässä on oikeasti vähennettävä oikea vuodenpäivä kustakin!
    for time_index in range(len(times)):
        print(time_index)
        day_of_year=ss.nemo_time_to_datetime(times[time_index]).timetuple().tm_yday
        d[time_index,:,:]=d[time_index,:,:]-climatology_data.variables[clim_var1][day_of_year-1,:,:]
    print("----Passed")
    data.close()
    if(just_one):
        times=[times[0]]
    for time_frame in range(len(times)):
        time=ss.nemo_time_to_datetime(times[time_frame])
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
        annotation=plt.annotate(time.strftime("%Y-%m-%d"),xy=(0.25, 0.95), xycoords='axes fraction',zorder=100)
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
        print("Image {} of (approx){} done.".format(running_number,len(files_working)*image_per_month[ss.interval]))
climatology_data.close()
        