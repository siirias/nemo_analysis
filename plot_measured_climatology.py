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
from smartseahelper import smh

output_dir = "D:\\Data\\Figures\\SmartSea\\"
output_dir_plus = ""
measurement_dir = "D:\\Data\\SmartSeaModeling\\Climatologies\\"
data_dir = "E:\\SmartSea\\climatologies\\"  
bathymetric_file = "D:\\Data\\ArgoData\\iowtopo2_rev03.nc"

fig_dpi = 300

mod_min_lat = 59.92485
mod_max_lat = 65.9080876
mod_min_lon = 16.40257
mod_max_lon = 25.8191425
mod_shape_lat = 360
mod_shape_lon = 340


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
serie_types = [ "SSS_5vs1_diff", "SSS_2vs1_diff",  "SBS_5vs1_diff", "SBS_2vs1_diff"]
serie_types = [ "SST_1vsABD1_diff", "SBS_1vsABD1_diff", "SSS_1vsABD1_diff"]
serie_types = [ "SBS_1vsABD1_diff", "SSS_1vsABD1_diff"]
serie_types = [ "SSS_1", "SBS_1", "SST_1", "SST_2", "SST_5"]

serie_types = [ "SST_2vs1_diff_special", "SST_5vs1_diff_special"]
serie_types = [ "SBS_1vsABD1_diff_test", "SSS_1vsABD1_diff_test", "SST_1vsABD1_diff_test"]
serie_types = [ "SBS_1vsABD1_diff_test"]
serie_types = [ "SSS_1", "SBS_1", "SST_1"]

data_sets = ["ABD", "A", "B", "D"]
data_sets = ["ABD"]
#data_sets = ["A", "B", "D"]
plot_bathymetry = False
plot_bathy_contours = True
plot_yearly_average = True
plot_daily_figures = False
comparison = True   # This one is set depending on do 
                    # the setup give name for another dataset
comparison_climatology = 'SDC'  # None, 'BNSC', 'SDC' 'TSO50' # if not none, overrides configuration comparison

the_proj = ccrs.PlateCarree()
# read configuration from a file.
# note: xonfig files are python statements,
# this doesn't check for malignant code or anything like that,
# so be vary.
set_configurations = open('climatology_dataset_types.txt').readlines()

# all these are defined later in config file
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

def interp_for_climatology(data, clim_data):
    interp_data = data.copy()
    lat = np.linspace(mod_min_lat,mod_max_lat,mod_shape_lat)  # This is a gludge, as the original lat lon ar bit weird.
    lon = np.linspace(mod_min_lon,mod_max_lon,mod_shape_lon)
    interp_data['y'] = lat
    interp_data['x'] = lon
    interp_data = interp_data.interp(y=clim_data['lat'],x=clim_data['lon'])
    return interp_data

for serie_type in serie_types:
    figure_size = (10,10)
    plot_area = [17.0, 26.0, 60.0, 66.0]
    bathy_max = 250.0
    bathy_levels = list(np.arange(0,bathy_max,50.0))  # number or list of numbers
    for data_set in data_sets:
        in_set = False
        for i in set_configurations:
            r = re.search("#SET (.*)",i.strip())
            if(r):
                in_set = r.groups()[0]
            elif(in_set == serie_type):
                print(i)
                exec(i)    
        data = xr.open_dataset(data_dir+file)
        if(len(file0)>0 or comparison_climatology):
            var0 = var  # change if comparing to measurements which have other names
            if(comparison_climatology == 'BNSC'):
                plot_yearly_average = True
                plot_daily_figures = False
                data0 = xr.open_dataset(measurement_dir+'BNSC_BothnianSea_inter20062015_TS.nc')
                data = interp_for_climatology(data,data0)
                set_name0 = 'BNSC'
                file0 = 'BNSC'
                if(var == 'SST'):
                    var0 = 'temperature_oan'
                if(var == 'SSS'):
                    var0 = 'salinity_oan'
                if(var == 'SBS'):
                    var0 = 'salinity_oan'
            elif(comparison_climatology == 'SDC'):
                plot_yearly_average = True
                plot_daily_figures = False
                data0 = xr.open_dataset(measurement_dir+'SDC_BAL_CLIM_T_1955_2014_00625_m.nc')
                data = interp_for_climatology(data,data0)
                set_name0 = 'SDC'
                file0 = 'SDC'
                if(var == 'SST'):
                    var0 = 'ITS-90 water temperature'
                if(var == 'SSS'):
                    var0 = 'Water body salinity'
                    data0.close()
                    data0 = xr.open_dataset(measurement_dir+\
                            'SDC_BAL_CLIM_S_1955_2014_00625_m.nc')
                if(var == 'SBS'):
                    var0 = 'Water body salinity'
                    data0.close()
                    data0 = xr.open_dataset(measurement_dir+\
                            'SDC_BAL_CLIM_S_1955_2014_00625_m.nc')
            elif(comparison_climatology == 'TSO50'):
                plot_yearly_average = True
                plot_daily_figures = False
                data0 = xr.open_dataset(measurement_dir+'tso50.nc')
                data0 = data0.rename({'latitude':'lat','longitude':'lon'})
                data0 = data0.transpose('time','depth','lat','lon')
                data = interp_for_climatology(data,data0)
                set_name0 = 'TSO50'
                file0 = 'TSO50'
                if(var == 'SST'):
                    var0 = 'temperature'
                if(var == 'SSS'):
                    var0 = 'salinity'
                if(var == 'SBS'):
                    var0 = 'salinity'
                    
            else:
                data0 = xr.open_dataset(data_dir+file0)
            comparison = True
        else:
            comparison = False
        try:
            os.mkdir(output_dir+output_dir_plus)
        except:
            pass
        #first plot the yearly average:
        if(plot_yearly_average):
            day_filters = {
                    "Year":slice(0,None),
                    "DJF":[slice(0,58),slice(337,368)],
                    "MAM":slice(59,151),
                    "JJA":slice(152,244),
                    "SON":slice(245,336)}
            day_filters0 = day_filters.copy() 
            for i in day_filters:
                fig = create_main_map(the_proj)
                if(comparison_climatology == None):
                    lat = np.linspace(mod_min_lat,mod_max_lat,mod_shape_lat)
                    lon = np.linspace(mod_min_lon,mod_max_lon,mod_shape_lon)
                    lon,lat = np.meshgrid(lon,lat)
                    d0_dat = data0[var0][:,:,:]
                else: # data has been modified, so get the lat lon there
                    lat = np.array(data0['lat'])
                    lon = np.array(data0['lon'])
                    d0_dat = data0[var0][:,0,:,:]
                    if(var == 'SBS'): # have to gather bottom layer
                        d0_dat = smh.get_bottom(None,\
                                np.ma.array(data0[var0],\
                                mask = np.isnan(data0[var0])))
                    day_filters0 = {
                            "Year":slice(0,None),
                            "DJF":[slice(0,2),slice(11,None)],
                            "MAM":slice(2,5),
                            "JJA":slice(5,8),
                            "SON":slice(8,11)}
                if(type(day_filters[i]) == list):
                    d = np.mean(np.concatenate(
                            (data[var][day_filters[i][0],:,:],
                             data[var][day_filters[i][1],:,:])),0)
                    d0_dat = np.mean(np.concatenate(
                            (d0_dat[day_filters0[i][0],:,:],
                             d0_dat[day_filters0[i][1],:,:])),0)
                else:
                    d = np.mean(data[var][day_filters[i],:,:],0)
                    d0_dat = np.mean(d0_dat[day_filters0[i],:,:],0)
#                if(not comparison):
#                    d0_dat = np.zeros(d.shape)
                #d = np.ones(d.shape)
                plt.pcolor(lon,lat,d0_dat,transform = the_proj, cmap = color_map, \
                           vmin = var_lims[0], vmax = var_lims[1])
                cb=plt.colorbar()
                cb.set_label(var0)
                plt.title("{}  {} \n {}".format(var_name,i,file0))
                filename = "Measured_{}_{}_{}.png".format(var,set_name0,i)
                plt.savefig(output_dir+output_dir_plus_means+filename,\
                                facecolor='w',dpi=fig_dpi,bbox_inches='tight')
        
                print("Saved: {}".format(output_dir+output_dir_plus_means+filename))
                plt.close()
        data.close()
#        if(comparison):
#            data0.close()
#    
    
