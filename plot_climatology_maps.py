# -*- coding: utf-8 -*-
"""
Created on Sun Nov 15 14:29:42 2020

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
import cmocean as cmo

output_dir = "D:\\Data\\Figures\\SmartSea\\"
output_dir_plus = ""
fig_dpi = 300
data_dir = "E:\\SmartSea\\climatologies\\"  
bathymetric_file = "D:\\Data\\ArgoData\\iowtopo2_rev03.nc"
all_variables = {"Temperature_monthly":"votemper",\
                 "Salinity_monthly":"vosaline",
                 "SSS":"SSS",
                 "SBS":"SBS",
                 "SST":"SST"
                 }

color_maps = {"Temperature_monthly":cmo.cm.thermal,\
              "Salinity_monthly":cmo.cm.haline,
             "SSS":cmo.cm.haline,
             "SBS":cmo.cm.haline,
             "SST":cmo.cm.thermal
              }
#serie_types = ["SBS_2vs1_diff", "SBS_5vs1_diff", "SBS_5vs2_diff"]
#serie_types = [ "SBS_2vs1_diff", "SBS_5vs2_diff","SSS_2vs1_diff", "SSS_5vs1_diff", "SSS_5vs2_diff"]
serie_types = [ "SST_2vs1_diff", "SST_5vs2_diff","SST_5vs1_diff"]
the_proj = ccrs.PlateCarree()
set_configurations = open('climatology_dataset_types.txt').readlines()

# all these are deined later in config file
file0 = None
file = None
output_dir_plus = None
output_dir_plus_means = None
set_name = None
set_name0 = None
color_map = None
var_lims= None
#var_name = "Temperature_monthly"
var_name = "SST"
var = all_variables[var_name]

def create_main_map(the_proj):
        fig=plt.figure(figsize=figure_size)
        plt.clf()
        ax = plt.axes(projection=the_proj)
        ax.set_extent(plot_area)
        ax.set_aspect('auto')
        ax.coastlines('10m',zorder=4)
        ax.add_feature(cfeature.NaturalEarthFeature('physical', 'land', '10m',\
                                                edgecolor='face', facecolor='#555560'))
        gl = ax.gridlines(crs=the_proj, draw_labels=True,
                  linewidth=2, color='gray', alpha=0.1, linestyle='-')
        gl.xlabels_top = False
        gl.ylabels_right = False
        
        # Then possible bathymetry or contours
        if plot_bathymetry or plot_bathy_contours:
            topodata = Dataset(bathymetric_file)
            topoin = topodata.variables['Z_WATER'][:]
            lons = topodata.variables['XT_I'][:]
            lats = topodata.variables['YT_J'][:]
            x=np.tile(lons,(lats.shape[0],1))
            y=np.tile(lats,(lons.shape[0],1)).T
            if(plot_bathy_contours):
                cn = plt.contour(x,y,-1*topoin,colors='k',vmin=0,vmax=bathy_max,\
                                 alpha=0.3,levels = bathy_levels, zorder = 5,\
                                 transform = the_proj)
                plt.clabel(cn,fmt='%1.0f')
            if(plot_bathymetry):
                plt.pcolor(x,y,-1*topoin,cmap=cmo.cm.deep,vmin=0,\
                           vmax=bathy_max, transform = the_proj)
                cb=plt.colorbar()
                cb.ax.invert_yaxis()
                cb.set_label('Depth (m)')
        return fig

for serie_type in serie_types:
    figure_size = (10,10)
    plot_area = [17.0, 26.0, 60.0, 66.0]
    plot_bathymetry = False
    plot_bathy_contours = True
    plot_yearly_average = True
    plot_daily_figures = False
    bathy_max = 250.0
    bathy_levels = list(np.arange(0,bathy_max,50.0))  # number or list of numbers
    mod_min_lat = 59.92485
    mod_max_lat = 65.9080876
    mod_min_lon = 16.40257
    mod_max_lon = 25.8191425
    mod_shape_lat = 360
    mod_shape_lon = 340
    data_set = "ABD"
    data_set = "D"
    # read configuration from a file.
    # note: xonfig files are python statements,
    # this doesn't check for malignant code or anything like that,
    # so be vary.
    in_set = False
    for i in set_configurations:
        r = re.search("#SET (.*)",i.strip())
        if(r):
            in_set = r.groups()[0]
        elif(in_set == serie_type):
            print(i)
            exec(i)    
#    # load the actual data:
#    if(serie_type == "SBS_5vs1_diff"):
#        file0 = "daily_average_{}001.nc".format(data_set)
#        file = "daily_average_{}005.nc".format(data_set)
#        output_dir_plus = "\\SBS_{}5vs1\\".format(data_set)
#        set_name = "{}005".format(data_set)
#        set_name0 = "{}001".format(data_set)
#        color_map = cmo.cm.diff
#        var_lims=[-0.5,0.5]
#        #var_name = "Temperature_monthly"
#        var_name = "SBS"
#        var = all_variables[var_name]
#    
#    if(serie_type == "SBS_2vs1_diff"):
#        file0 = "daily_average_{}001.nc".format(data_set)
#        file = "daily_average_{}002.nc".format(data_set)
#        output_dir_plus = "\\SBS_{}2vs1\\".format(data_set)
#        set_name = "{}002".format(data_set)
#        set_name0 = "{}001".format(data_set)
#        color_map = cmo.cm.diff
#        var_lims=[-0.5,0.5]
#        #var_name = "Temperature_monthly"
#        var_name = "SBS"
#        var = all_variables[var_name]
#
#    if(serie_type == "SBS_5vs2_diff"):
#        file0 = "daily_average_{}002.nc".format(data_set)
#        file = "daily_average_{}005.nc".format(data_set)
#        output_dir_plus = "\\SBS_{}5vs2\\".format(data_set)
#        set_name = "{}005".format(data_set)
#        set_name0 = "{}002".format(data_set)
#        color_map = cmo.cm.diff
#        var_lims=[-0.5,0.5]
#        #var_name = "Temperature_monthly"
#        var_name = "SBS"
#        var = all_variables[var_name]
#    if(serie_type == "SSS_5vs1_diff"):
#        file0 = "daily_average_{}001.nc".format(data_set)
#        file = "daily_average_{}005.nc".format(data_set)
#        output_dir_plus = "\\SSS_{}5vs1\\".format(data_set)
#        set_name = "{}005".format(data_set)
#        set_name0 = "{}001".format(data_set)
#        color_map = cmo.cm.diff
#        var_lims=[-0.5,0.5]
#        #var_name = "Temperature_monthly"
#        var_name = "SSS"
#        var = all_variables[var_name]
#    
#    if(serie_type == "SSS_2vs1_diff"):
#        file0 = "daily_average_{}001.nc".format(data_set)
#        file = "daily_average_{}002.nc".format(data_set)
#        output_dir_plus = "\\SSS_{}2vs1\\".format(data_set)
#        set_name = "{}002".format(data_set)
#        set_name0 = "{}001".format(data_set)
#        color_map = cmo.cm.diff
#        var_lims=[-0.5,0.5]
#        #var_name = "Temperature_monthly"
#        var_name = "SSS"
#        var = all_variables[var_name]
#
#    if(serie_type == "SSS_5vs2_diff"):
#        file0 = "daily_average_{}002.nc".format(data_set)
#        file = "daily_average_{}005.nc".format(data_set)
#        output_dir_plus = "\\SSS_{}5vs2\\".format(data_set)
#        set_name = "{}005".format(data_set)
#        set_name0 = "{}002".format(data_set)
#        color_map = cmo.cm.diff
#        var_lims=[-0.5,0.5]
#        #var_name = "Temperature_monthly"
#        var_name = "SSS"
#        var = all_variables[var_name]
#    if(serie_type == "SST_5vs1_diff"):
#        file0 = "daily_average_{}001.nc".format(data_set)
#        file = "daily_average_{}005.nc".format(data_set)
#        output_dir_plus = "\\SST_{}5vs1\\".format(data_set)
#        set_name = "{}005".format(data_set)
#        set_name0 = "{}001".format(data_set)
#        color_map = 'RdBu_r'
#        var_lims=[-2.5,2.5]
#        #var_name = "Temperature_monthly"
#        var_name = "SST"
#        var = all_variables[var_name]
#    if(serie_type == "SST_2vs1_diff"):
#        file0 = "daily_average_{}001.nc".format(data_set)
#        file = "daily_average_{}002.nc".format(data_set)
#        output_dir_plus = "\\SST_{}2vs1\\".format(data_set)
#        set_name = "{}002".format(data_set)
#        set_name0 = "{}001".format(data_set)
#        color_map = 'RdBu_r'
#        var_lims=[-2.5,2.5]
#        #var_name = "Temperature_monthly"
#        var_name = "SST"
#        var = all_variables[var_name]
#    if(serie_type == "SST_5vs2_diff"):
#        file0 = "daily_average_{}002.nc".format(data_set)
#        file = "daily_average_{}005.nc".format(data_set)
#        output_dir_plus = "\\SST_{}5vs2\\".format(data_set)
#        set_name = "{}005".format(data_set)
#        set_name0 = "{}002".format(data_set)
#        color_map = 'RdBu_r'
#        var_lims=[-2.5,2.5]
#        #var_name = "Temperature_monthly"
#        var_name = "SST"
#        var = all_variables[var_name]
#
#
#
    data = xr.open_dataset(data_dir+file)
    data0 = xr.open_dataset(data_dir+file0)
    try:
        os.mkdir(output_dir+output_dir_plus)
    except:
        pass
    #first plot the yearly average:
    if(plot_yearly_average):
        day_filters = {
                "Year":slice(0,-1),
                "DJF":[slice(0,58),slice(337,368)],
                "MAM":slice(59,151),
                "JJA":slice(152,244),
                "SON":slice(245,336)}
        for i in day_filters:
            fig = create_main_map(the_proj)             
            lat = np.linspace(mod_min_lat,mod_max_lat,mod_shape_lat)
            lon = np.linspace(mod_min_lon,mod_max_lon,mod_shape_lon)
            lon,lat = np.meshgrid(lon,lat)
            if(type(day_filters[i]) == list):
                d = np.mean(np.concatenate(
                        (data[var][day_filters[i][0],:,:],
                         data[var][day_filters[i][1],:,:])),0)
                d0 = np.mean(np.concatenate(
                        (data0[var][day_filters[i][0],:,:],
                         data0[var][day_filters[i][1],:,:])),0)
            else:
                d = np.mean(data[var][day_filters[i],:,:],0)
                d0 = np.mean(data0[var][day_filters[i],:,:],0)
            #d = np.ones(d.shape)
            plt.pcolor(lon,lat,d-d0,transform = the_proj, cmap = color_map, \
                       vmin = var_lims[0], vmax = var_lims[1])
            cb=plt.colorbar()
            cb.set_label('Difference')
            plt.title("{} Diff, Whole Year \n {}-{}".format(var_name,file,file0))
            filename = "{}_{}vs{}_{}.png".format(var,set_name,set_name0,i)
            plt.savefig(output_dir+output_dir_plus_means+filename,\
                            facecolor='w',dpi=fig_dpi,bbox_inches='tight')
    
            plt.close()
    if(plot_daily_figures):
        for time_step in range(365):
            #First let's plot the general map:
            fig = create_main_map(the_proj)             
            #color_map = color_maps[var_name]
            #lat_orig = xr.where(data["nav_lat"]!=0.0,data["nav_lat"],None)
            #lon_orig = xr.where(data["nav_lon"]!=0.0,data["nav_lon"],None)
            # the original ones are just broken for this purpose,
            # have to do this manually for now:
            lat = np.linspace(mod_min_lat,mod_max_lat,mod_shape_lat)
            lon = np.linspace(mod_min_lon,mod_max_lon,mod_shape_lon)
            lon,lat = np.meshgrid(lon,lat)
            #d = xr.where(not np.isnan(data[var][0,0,:,:]),data[var][0,0,:,:],None)
            #d = xr.where(data[var][0,0,:,:] != 0.0, data[var][0,0,:,:], None)
            d = data[var][time_step,:,:]
            d0 = data0[var][time_step,:,:]
            #d = np.ones(d.shape)
            plt.pcolor(lon,lat,d-d0,transform = the_proj, cmap = color_map, \
                       vmin = var_lims[0], vmax = var_lims[1])
            cb=plt.colorbar()
            cb.set_label('Difference')
            plt.title("{} Diff, day: {:03d} \n {}-{}".format(var_name,time_step,file,file0))
            filename = "{}_{}vs{}_{:03d}.png".format(var,set_name,set_name0,time_step)
            plt.savefig(output_dir+output_dir_plus+filename,\
                            facecolor='w',dpi=fig_dpi,bbox_inches='tight')
        
            plt.close()
    data.close()
    data0.close()
    
    
