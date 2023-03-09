# -*- coding: utf-8 -*-
"""
Created on Thu Jul 29 19:01:31 2021

@author: siirias
"""

import os
import matplotlib.pyplot as plt
import numpy as np
import cartopy.crs as ccrs
import cartopy.feature as cfeature

from netCDF4 import Dataset
import xarray as xr
import cmocean as cmo

output_dir = "C:\\Data\\Figures\\SmartSeaNew\\climatologies\\"
main_data_dir = "D:\\SmartSea\\new_dataset\\derived_data\\test\\"  
bathymetric_file = "C:\\Data\\ArgoData\\iowtopo2_rev03.nc"

color_maps = {"SBT":cmo.cm.thermal,
             "SSS":cmo.cm.haline,
             "SBS":cmo.cm.haline,
             "SST":cmo.cm.thermal,
             "dSBT":cmo.cm.balance,
             "dSSS":cmo.cm.delta,
             "dSBS":cmo.cm.delta,
             "dSST":cmo.cm.balance,
              }
shown_units = {"SSS":"g/kg",
               "SBS":"g/kg",
               "SST":"°C",
               "SBT":"°C",
               "ICE_C":""}
var_name = {"SST":"SST_mean",
            "SBT":"SBT_mean",
            "SSS":"SSS_mean",
            "SBS":"SBS_mean",
            "ICE_C":""}
model_names = {'A':'model 1',
               'B':'model 2',
               'D':'model 3',
               'ABD':'all models'}
var_lims_list = {'SST':(None,None),
                 'dSST':(-5.0,5.0),
                 'SBT':(None,None),
                 'dSBT':(-5.0,5.0),
                 'SSS':(None,None),
                 'dSSS':(-1.2,1.2),
                 'SBS':(None,None),
                 'dSBS':(-1.2,1.2),}

var_lims_list = {'SST':(3.0,15.0),
                 'dSST':(-5.0,5.0),
                 'SBT':(3.0,15.0),
                 'dSBT':(-5.0,5.0),
                 'SSS':(0.0,8.0),
                 'dSSS':(-1.2,1.2),
                 'SBS':(0.0,8.0),
                 'dSBS':(-1.2,1.2),}

figure_size = (10,10)
fig_dpi = 300
mod_min_lat = 59.92485
mod_max_lat = 65.9080876
mod_min_lon = 16.40257
mod_max_lon = 25.8191425
mod_shape_lat = 360
mod_shape_lon = 340
plot_area = [17.0, 26.0, 60.0, 66.0]
bathy_max = 250.0
bathy_levels = list(np.arange(0,bathy_max,50.0))  # number or list of numbers
lat = np.linspace(mod_min_lat,mod_max_lat,mod_shape_lat)
lon = np.linspace(mod_min_lon,mod_max_lon,mod_shape_lon)
lon,lat = np.meshgrid(lon,lat)


the_proj = ccrs.PlateCarree()

#clim_sets = ['c30v', 'c70v']
clim_sets = ['c30v', 'c70v']
var_list = ['SSS', 'SBS', 'SST', 'SBT']
#var_list = ['ICE_C']
#var_list = ['SST', 'SBT']
var_lims = [None, None]
model_sets = ['ABD','A','B','D'] # 'A','B','D','ABD'
#model_sets = ['ABD'] # 'A','B','D','ABD'
#data_sets = ['RCP45', 'RCP85'] #'reference', 'RCP45', 'RCP85'
data_sets = ['reference','RCP45', 'RCP85'] #'reference', 'RCP45', 'RCP85'
time_scale = 'd' # monthly or daily data
mean_types = ['Year', 'DJF', 'MAM','JJA','SON'] #'Year', 'DJF', 'MAM','JJA','SON'
mean_types= ['Year']
compare_to_reference = False
plot_bathymetry = False
plot_bathy_contours = True
color_map = cmo.cm.haline
close_windows = True
def create_main_map(the_proj):
        details = '10m'
        fig=plt.figure(figsize=figure_size)
        plt.clf()
        ax = plt.axes(projection=the_proj)
        ax.set_extent(plot_area)
        ax.set_aspect('auto')
        ax.coastlines(details,zorder=4)
        ax.add_feature(cfeature.NaturalEarthFeature('physical', 'land', details,\
                                                edgecolor='face', facecolor='#e0e0f0'))
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

def get_dataset(model_set, data_set, variable, clim_set):
    setnames = {'reference':'001',
                'RCP45':'002',
                'RCP85':'005'}
    if model_set =='ABD':
        models = ['A','B','D']
    else:
        models = [model_set]
    tmp_dat = 0
    # if(clim_set == 'reference'):
    #     clim_set = 'c30v' # doesn't really matter, 
    #                       # now we just want the history data, either is fine.
    #     data_set = 'reference'
    for m in models:
        data_dir = "{}\\{}\\".format(main_data_dir, clim_set)
        data_filename = "{}_climatology_{}_{}_{}.nc".format(\
                                    clim_set, 
                                    m + setnames[data_set],
                                    'd',
                                    variable)
        if(type(tmp_dat) == int): #Just to check if it is empty
            tmp_dat = xr.open_dataset(data_dir+data_filename)\
                [var_name[variable]][:,:,:]
        else:
            tmp_dat += xr.open_dataset(data_dir+data_filename)\
                [var_name[variable]][:,:,:]
    tmp_dat = tmp_dat/len(models)
    return tmp_dat


    
def get_mean(data, mean_type='Year'):
    #mean_type = 'Year', 'DJF', 'MAM','JJA','SON'
        if(data.shape[0] == 12): # this is monthly values
            day_filters = {
                    "Year":slice(0,None),
                    "DJF":[slice(0,2),slice(11,None)],
                    "MAM":slice(2,5),
                    "JJA":slice(5,8),
                    "SON":slice(8,11)}
        else: # let's assume daily values
            day_filters = {
                    "Year":slice(0,None),
                    "DJF":[slice(0,58),slice(337,368)],
                    "MAM":slice(59,151),
                    "JJA":slice(152,244),
                    "SON":slice(245,336)}    
        if(type(day_filters[mean_type]) == list):
            the_data = np.concatenate(
                    (data[day_filters[mean_type][0],:,:],
                     data[day_filters[mean_type][1],:,:]))
        else:
            the_data = data[day_filters[mean_type],:,:]
        return np.mean(the_data,0)
    

## ACTUAL SCRIPT
for clim_set in clim_sets: # which climatology period
    for var in var_list:  # which variable
        output_dir_plus = "\\"+"values"+"\\"
        if(compare_to_reference):
            output_dir_plus = "\\"+"vs_reference"+"\\"
        output_dir_plus += "\\"+clim_set+"\\"+var+'\\'
        if(not os.path.exists(output_dir+output_dir_plus)):
            os.makedirs(output_dir+output_dir_plus)
        for model_set in model_sets: # which A,B, D, ABD
            for data_set in data_sets: # RCP
                for mean_type in mean_types: # time period
                    extra_type = ''
                    var_lims = var_lims_list[var]
                    create_main_map(the_proj)
                    data = get_dataset(model_set,data_set,var, clim_set)
                    d = get_mean(data,mean_type)
                    c_map = color_maps[var]
                    if(compare_to_reference):
                        extra_type = 'diff'
                        r_data = get_dataset(model_set,'reference',var, clim_set)
                        r_d = get_mean(r_data, mean_type)
                        d = d -r_d
                        var_lims = var_lims_list['d'+var]
                        c_map = color_maps['d'+var]
                    plt.pcolor(lon,lat,d,transform = the_proj, cmap = c_map, \
                               vmin = var_lims[0], vmax = var_lims[1])
                    cb = plt.colorbar()
                    cb.set_label('{} {}({})'.format(extra_type, var, shown_units[var]))
                    plt.title("{} {} {}, {}".format(
                        var,
                        data_set,
                        mean_type,
                        model_names[model_set]))
                    filename = "{}_{}_{}_{}_{}_{}.png".format(var,
                                                         clim_set,
                                                         data_set,
                                                         model_set,
                                                         mean_type,
                                                         extra_type)
                    plt.savefig(output_dir+output_dir_plus+filename,\
                                    facecolor='w',dpi=fig_dpi,bbox_inches='tight')
                    print("Saved: {}".format(output_dir+output_dir_plus+filename))
                    if(close_windows):
                        plt.close()
