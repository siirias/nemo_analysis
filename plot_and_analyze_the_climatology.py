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
import importlib
import param_sets_for_plotting
months = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']

param_sets =[\
             {  'var1':'icecon',\
                'name_marker':'A005',\
                'other_name_marker':'A001',\
                'slice_wanted':-1,\
                'compare_two':True,\
                'output_dir_end':"/tmp/images/ice/",\
                'cludge':True\
             }
            ]  # actually parameters are set in next line. this is just an example.
importlib.reload(param_sets_for_plotting)
param_sets = param_sets_for_plotting.temperature_param_sets

for param_set in param_sets:
#    param_set['name_marker'] = param_set['name_marker'].replace('A','C')  # a cludge to quickly change the set.
#    param_set['other_name_marker'] = param_set['other_name_marker'].replace('A','C')  # a cludge to quickly change the set.
    print(param_set)
    name_marker=param_set['name_marker']  #this one tells which dataseries is handled.
    compare_two=param_set['compare_two']
    other_name_marker=param_set['other_name_marker']
    land_color = (0.8, 0.83, 0.8)
    sea_color = (0.5, 0.5, 0.7)  # meaningful only on sea-area where no data present
    var1 = param_set['var1']  #'VEL', 'SST', 'SSS', 'votemper', 'vosaline'
    slice_wanted = param_set['slice_wanted']  # if 3d grid, which layer to take. if -1, bottommost.
    ss=smh()
    ss.grid_type='T'
    ss.interval='m'
    ss.main_data_folder= ss.root_data_in+"/OUTPUT{}/".format(name_marker)
    datadir = ss.root_data_out+param_set['output_dir_end'] #where everyt output is stored
    extra =""
    other_extra = ""
    cludge = param_set['cludge']
    
    if(not '1' in name_marker):
        extra = 'c20v_'
    if(not '1' in other_name_marker):
        other_extra = 'c20v_'
    
    climatology_file=ss.root_data_out+'derived_data/test/{}climatology_{}_{}_{}.nc'.format(\
                extra,name_marker,ss.interval,var1)  #this one tells which dataseries is handled.
    
    
    
    
    
    
    series_name='climatology{}'.format(name_marker)
    if(compare_two):
        other_climatology_file=ss.root_data_out+\
            'derived_data/test/{}climatology_{}_{}_{}.nc'.\
            format(other_extra,other_name_marker,ss.interval,var1)  #this one tells which dataseries is handled.
        series_name='climatology{}vs{}'.format(name_marker,other_name_marker)
    if cludge:    
        climatology_file=ss.root_data_out+\
                'derived_data/test/monthly_climatology_{}_d_icecon.nc'.format(name_marker) #GLUDGE!
        other_climatology_file=ss.root_data_out+\
                'derived_data/test/monthly_climatology_{}_d_icecon.nc'.format(other_name_marker) #GLUDGE!
    
    def get_bottom(grid):
        # grid is supposed to be masked array, Time, D,Lat,Lon
        # The idea in this is to shifht the mask one layer up,
        # and find the values which are masked in one (and only one) of
        # these masks. 
        full_shape = grid.shape
        bottom_layers = np.zeros((  full_shape[0],\
                                    full_shape[2],
                                    full_shape[3]))
        bottom_layers = np.ma.masked_array(bottom_layers,False)
        mask_roll = np.roll(grid.mask,-1,1) # move mask values one up.
        mask_roll[:,-1,:,:] = True  # And mark bottom most mask as True.
                                   # This to get bottom values if there are no mask at end
        grid.mask = ~(grid.mask ^ mask_roll)
        bottom_layers = np.sum(grid,1)
        values = np.array(np.sum(~grid.mask,1),bool)  # used to get the mask
        bottom_layers.mask = ~values
        return bottom_layers
    
    #Let's try to plot soemthing to start with:
    lon_min=16;lat_min=60;lon_max=26.01;lat_max=66.01;
    running_number=0
    show_contours=True
    contour_step = 0.5
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
        plotted_depth=-10000. #meters
        var2=None
        show_contours=False
    if var1 in ['SST','votemper']:
        var_min=-0.5
        var_max=18.0
        var1_cm=cmocean.cm.thermal
        ss.grid_type='T'  #U,V,T
        if(compare_two):
            var_min=-2.
            var_max=2.
            var1_cm='coolwarm'
    if var1 in ['SSS','vosaline']:
        var_min=0.0
        var_max=8.0
        var1_cm=cmocean.cm.haline
        ss.grid_type='T'  #U,V,T
        if(compare_two):
            var_min=-2.
            var_max=2.
            var1_cm='RdGy'
    if var1 in ['SSH']:
        var_min=1.0 
        var_max=6.0 
        var1_cm=cmocean.cm.haline 
        ss.grid_type='T'  #U,V,T
    
    if var1 in ['icecon']:
        var_min=0.0 
        var_max=1.0 
        var1_cm=cmocean.cm.ice 
        ss.grid_type='T'  #U,V,T
        if(compare_two):
            var_min=-1.
            var_max=1.
            var1_cm='RdGy'
    
    
    
    
    
    font_size=10.0
    resolution='h'
    projection='laea'
    just_one=False
    
        
    
    
    #first setup the main map:
    fig=plt.figure(figsize=(6,6))
    if projection in ['laea']:
        lat_0=0.5*(lat_min+lat_max)
        lon_0=0.5*(lon_min+lon_max)
        bmap = Basemap(llcrnrlon=lon_min,llcrnrlat=lat_min,urcrnrlon=lon_max,urcrnrlat=lat_max, \
                       lat_0=lat_0, lon_0=lon_0,resolution = resolution, \
                       projection=projection,fix_aspect=False)  #resolution c,l,i,h,f
    elif projection in ['merc','cyl']:
        bmap = Basemap(llcrnrlon=lon_min,llcrnrlat=lat_min,urcrnrlon=lon_max,urcrnrlat=lat_max, \
                       resolution = resolution,
                       projection=projection,fix_aspect=False)  #resolution c,l,i,h,f
        
    bmap.drawmapboundary(fill_color=sea_color, zorder = -10)
    bmap.fillcontinents(color=land_color, lake_color=sea_color,zorder = -9)
    bmap.drawcoastlines(zorder=21,linewidth=0.5,color='gray')
    #bmap.fillcontinents([0.9,0.9,0.9],lake_color=[0.85,0.85,0.85],zorder=20) 
    bmap.drawparallels(np.arange(lat_min,lat_max,1.),linewidth=1,zorder=50,\
            labels=[True,False,False,False],dashes=[1,0],color="#00000020",fontsize=10)        
    bmap.drawmeridians(np.arange(lon_min,lon_max,2.),linewidth=1,zorder=50,\
            labels=[False,False,False,True],dashes=[1,0],color="#00000020",fontsize=10)        
    is_first=True
    data=Dataset(climatology_file)
    if(compare_two):
        other_data=Dataset(other_climatology_file)
    if(var1 in ['SST','SSH','SSS','votemper','vosaline','icecon','ivevolume']):
        d=data.variables['{}_mean'.format(var1)][:]
        if(len(d.shape)>3): # means this is 3d data, we only want one slice
            if(slice_wanted>=0):  # Just a regular slice
                d = d[:,slice_wanted,:,:]
            else:  #Now we want the bottom-most
                d = get_bottom(d)            
        if(compare_two):
            other_d=other_data.variables['{}_mean'.format(var1)][:]
            if(len(other_d.shape)>3): # means this is 3d data, we only want one slice
                if(slice_wanted>=0):  # Just a regular slice
                    other_d = other_d[:,slice_wanted,:,:]
                else:
                    other_d = get_bottom(other_d)
            #fix to get values fromplaces where one is masked, and the other is not
            mask_changed = d.mask ^ other_d.mask
            other_d.data[mask_changed & other_d.mask == True] = 0.0
            d.data[mask_changed & d.mask == True] = 0.0
            d.mask[mask_changed] = False  # set the mask value of these off.
            other_d.mask[mask_changed] = False  # set the mask value of these off.
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
    if(ss.interval == 'd'):
        times=data.variables['day_of_year'][:].data
    else:
        times=data.variables['month_of_year'][:].data
        times = times[:12]  # gludge to get over a bug. figure out later.
    data.close()
    if compare_two:
        other_data.close()
    if(just_one):
        times=[times[0]]
    for time_frame in range(len(times)):
        time=time_frame
        tmp_lon,tmp_lat=bmap(lons,lats)
        #cb=plt.colorbar(fraction=0.027, pad=0.01)
        value_range = np.max(d[time_frame,:,:])-np.min(d[time_frame,:,:])
        print(time_frame, value_range)
        colors_fig=bmap.pcolormesh(tmp_lon,tmp_lat,d[time_frame,:,:],vmin=var_min,vmax=var_max,zorder=-3,cmap=var1_cm)
        if is_first:
            cb=plt.colorbar()
            cb.set_clim(vmin=var_min,vmax=var_max)
            cb.ax.tick_params(labelsize=font_size)
            is_first=False
        if(var2 is not None):
            ice_fig=bmap.pcolormesh(tmp_lon,tmp_lat,ice_d[time_frame,:,:],vmin=var2_min,vmax=var2_max,zorder=16,cmap=var2_cm)
    
        if show_contours and value_range>0.0001:
            cont_fig=bmap.contour(tmp_lon,tmp_lat,d[time_frame,:,:],\
                                levels = np.arange(var_min, var_max, contour_step),\
                                vmin=var_min,vmax=var_max,\
                                zorder=15,colors='black',linewidths=0.5,alpha=0.5)
            cont_labels=plt.clabel(cont_fig,inline=1,fontsize=5,fmt="%1.1f")
    
        if(ss.interval == 'd'):
            annotation=plt.annotate("Year day: {}".format(time_frame+1),\
                            xy=(0.25, 0.95), xycoords='axes fraction',zorder=100)
        if(ss.interval == 'm'):
            annotation=plt.annotate("{}".format(months[time_frame%12]),\
                            xy=(0.25, 0.95), xycoords='axes fraction',zorder=100)
        layer_info = ""
        if(slice_wanted > 0):
            layer_info = "_D{}_".format(slice_wanted)
        if(slice_wanted<0):
            layer_info = "_bottom_"
        plt.savefig("{}{}{}{}{:05d}.png".format(  datadir,\
                                                var1,\
                                                series_name,\
                                                layer_info,\
                                                running_number),\
                                                facecolor='w',dpi=300)
        if(not just_one):
            #clean up the changing things, so we don't have to do everything again, just these:
            colors_fig.remove()
            if(var2 is not None):
                ice_fig.remove()
            annotation.remove()
            if show_contours and value_range > 0.0001:
                for i in cont_labels:
                    i.remove()
                for i in cont_fig.collections:
                    i.remove()
    #            plt.close(fig)
        running_number+=1
        print("Image {} done.".format(running_number))
    plt.close('all')    
