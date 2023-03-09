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

#output_dir = "D:\\Data\\Figures\\SmartSea\\"
#output_dir_plus = ""
#measurement_dir = "D:\\Data\\SmartSeaModeling\\Climatologies\\"
#data_dir = "E:\\SmartSea\\climatologies\\"  
#bathymetric_file = "D:\\Data\\ArgoData\\iowtopo2_rev03.nc"
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
        gl.top_labels = False
        gl.right_labels = False
        
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

output_dir = "C:\\Data\\Figures\\SmartSeaNew\\"
output_dir_plus = ""
data_dir = "D:\\SmartSea\\new_dataset\\A001\\"  
data_dir0 = "D:\\SmartSea\\Final\\SmartSea_NEMO_A001\\A001\\"  
bathymetric_file = "D:\\Data\\ArgoData\\iowtopo2_rev03.nc"

all_variables = {"Temperature_monthly":"votemper",\
                 "Salinity_monthly":"vosaline",
                 "SSS":"SSS",
                 "SBS":"SBS",
                 "SST":"SST",
                 "SBT":"votemper",
                 "ICE_C":"icecon",
                 "ICE_T":"icethic",
                 "ICE_V":"icethic",
                 "OXY_B":"oxy_bottom",
                 "MixedLayerD":"mldr10_1"
                 }

color_maps = {"Temperature_monthly":cmo.cm.thermal,\
              "Salinity_monthly":cmo.cm.haline,
             "SSS":cmo.cm.haline,
             "SBS":cmo.cm.haline,
             "SST":cmo.cm.thermal,
             "ICE_C":cmo.cm.ice,
             "ICE_T":cmo.cm.ice,
             "ICE_V":cmo.cm.ice,
             "OXY_B":cmo.cm.oxy,
             "MixedLayerD":cmo.cm.deep
              }

shown_units = {"SSS":"-",
               "SBS":"-",
               "SST":"°C",
               "SBT":"°C",
               "ICE_C":"%",
               "ICE_T":"m",
               "ICE_V":"m",
               "OXY_B":"ml/l",   # reasonable limits perhaps 0..15?
               "MixedLayerD":"m"
               }

fig_dpi = 300
close_figures = True  # set to False to keep figures open.
mod_min_lat = 59.92485
mod_max_lat = 65.9080876
mod_min_lon = 16.40257
mod_max_lon = 25.8191425
mod_shape_lat = 360
mod_shape_lon = 340

plot_bathymetry = False
plot_bathy_contours = False
plot_daily_figures = False
comparison = False   # This one is set depending on do 
                    # the setup give name for another dataset
#comparison_climatology = None
specific_depth = 0.0    # set target depth in meters. 0.0 is ignored. 
                        # Note: this does something only on SBT or SBS.

# Plot Empty map is for plotting some measuring points etc in the same format than rest of the figures.
plot_empty_map = False
tmp_col = '#a00000'
points_to_plot = {"F16":[63.52,21.10,tmp_col],"BO5":[64.19, 22.90,tmp_col],
                  "BO3":[64.31,22.35,tmp_col], "F9":[64.71, 22.07,tmp_col],
                  "F3":[65.17, 23.24,tmp_col] }
tmp_col = '#00a000'
#F64", "SR5", "MS4", "C3", "US5B" 
points_to_plot.update( {"F64":[60.19,19.15,tmp_col],"SR5":[61.09, 19.60,tmp_col],
                  "MS4":[62.09,18.54,tmp_col], "C3":[62.66, 18.96,tmp_col],
                  "US5B":[62.59, 19.99,tmp_col] })


the_proj = ccrs.PlateCarree()
years = list(range(1978,1984+1))
offset=0
for the_year in years:
    file0 = "NORDIC-GOB_1d_{y}0101_{y}1231_grid_T.nc".format(y=the_year)
    file = "NORDIC-GOB_1d_{y}0101_{y}1231_grid_T.nc".format(y=the_year)
    name_num_plus = 365*offset
    offset+=1
#    output_dir_plus = "\\test_diff\\"
#    output_dir_plus = "\\SBS_nocont\\"
    output_dir_plus = "\\ICE_C\\"
    set_name = "A001"
    set_name0 = "A001_old"
#    var_lims=[-2.0,2.0]
    #var_lims=[0.0,7.0]
#    var_lims=[0.0,80.0]
    var_lims = [0.0,1.0]
    var_name = "ICE_C"
    var = all_variables[var_name]
    if(comparison):
        color_map = cmo.cm.diff
    else:
        color_map = color_maps[var_name]
    
    
    
    
    figure_size = (10,10)
    plot_area = [17.0, 26.0, 60.0, 66.0]
    bathy_max = 250.0
    bathy_levels = list(np.arange(0,bathy_max,50.0))  # number or list of numbers
    in_set = False
    # slight glugdge next. SBS is in daily files, but other salinity
    # depths are not. if SBS is used to get specific depth, the file
    # mus be switchedfrom daily, to monthly file:
    if(specific_depth != 0.0 and var_name == 'SBS'):
        file = re.sub('daily','monthly',file)
        file0 = re.sub('daily','monthly',file0)
        var = 'vosaline'
        
    data = xr.open_dataset(data_dir+file)
    try:
        depth_axis = np.array(data['deptht'])
    except:
        depth_axis = None
    if(len(file0)>0 and comparison):
        var0 = var  # change if comparing to measurements which have other names
        data0 = xr.open_dataset(data_dir0+file0)
        try:
            depth_axis0 = np.array(data0['deptht'])
        except:
            depth_axis0 = None
        comparison = True
    else:
        comparison = False
    try:
        os.mkdir(output_dir+output_dir_plus)
    except:
        pass
    if(plot_empty_map):
            fig = create_main_map(the_proj)             
            cb=plt.colorbar()
            cb.set_ticks([0.0,0.5,1.0],True)  # just to get similar numbers
                                              # than comparable filled plot
                                              # to keep same aspects
            cb.set_label("none")
            plt.title("Empty \n measuring points")
            for p in points_to_plot:
                p_x = [points_to_plot[p][1]]
                p_y = [points_to_plot[p][0]]
                p_name = p
                p_color = points_to_plot[p][2]
                plt.plot(p_x,p_y, 
                         marker = '.', color = p_color, markersize = 20,\
                         transform = the_proj)


            filename = "empty_{}_{}.png".format(var,set_name)
            plt.savefig(output_dir+output_dir_plus+filename,\
                            facecolor='w',dpi=fig_dpi,bbox_inches='tight')
            print("Saved: {}, from {} ".format(\
                  output_dir+output_dir_plus+filename,
                  file))
            if(close_figures):
                plt.close()
            plot_empty_map = False
        
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
            if(var_name != "ICE_V"):  # This has to be calculated from consentration and thickness
                d = data[var][time_step,:,:]
                if(comparison):
                    d0 = data0[var][time_step,:,:]
                else:
                    d0 = np.zeros(d.shape)
            else:
                d = data[var][time_step,:,:]
                d = d*data[all_variables['ICE_C']][time_step,:,:]
                if(comparison):
                    d0 = data0[var][time_step,:,:]
                    d0 = d0*data0[all_variables['ICE_C']][time_step,:,:]
                else:
                    d0 = np.zeros(d.shape)
            #d = np.ones(d.shape)
            plt.pcolor(lon,lat,d-d0,transform = the_proj, cmap = color_map, \
                       vmin = var_lims[0], vmax = var_lims[1],shading='auto')
            cb=plt.colorbar()
            if(comparison):
                cb.set_label('Difference,({})'.format(shown_units[var_name]))
                plt.title("{} Diff, day: {:03d} \n {}-{}".format(var_name,time_step,file,file0))
                filename = "{}_{}vs{}_{:03d}.png".format(var,set_name,set_name0,time_step+name_num_plus)
            else:
                cb.set_label('{}({})'.format(var, shown_units[var_name]))
                plt.title("{} ,day: {:03d} \n {}".format(var_name,time_step,file))
                filename = "{}_{}_{:03d}.png".format(var,set_name,time_step+name_num_plus)
                
            plt.savefig(output_dir+output_dir_plus+filename,\
                            facecolor='w',dpi=fig_dpi,bbox_inches='tight')
            print("Saved: {}, from {} (day {})".format(\
                  output_dir+output_dir_plus+filename,
                  file,
                  time_step))
            if(close_figures):
                plt.close()
    data.close()
    if(comparison):
        data0.close()
    
    
