#!/usr/bin/env ipython
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
import netCDF4 as nc4
from netCDF4 import Dataset
from smartseahelper import smh
import os
import cmocean
import pandas as pd
import xarray as xr
import iris
import cf_units

ss=smh()
ss.grid_type='T'
ss.interval='m'
#folder_start='OUTPUT'
#name_markers=['new_REANALYSIS']
#name_markers=['REANALYSIS','A001','D001','C001']
name_markers=['A001','B001','D001','A002','B002','D002', 'A005','B005','D005']
#variables = ['SSS', 'SST', 'SSH_inst']
variables= ['vosaline']
#variables= ['vosaline','votemper']
make_climatology = True
#name_marker='A001'  #this one tells which dataseries is handled.
climatology_time_slots=366
if(ss.interval=='m'):
    climatology_time_slots=12
    
for name_marker in name_markers:
    folder_start=''
    ss.save_interval='year'   # 'month'
    ss.file_name_format='NORDIC-GOB_1{}_{}_{}_grid_{}.nc'
    print("Processing {}".format(name_marker))
    for var1 in variables:
        
        if '1' in name_marker: #the 001 series are hindcasts, all other scenarios
            startdate=datetime.datetime(1975,1,1)
            enddate=datetime.datetime(2005,12,31)
        elif 'REANALYSIS' in name_marker:
            startdate=datetime.datetime(1980,1,1)
            enddate=datetime.datetime(2012,12,31)
            ss.save_interval='year'
            folder_start=''
            if 'new_' in name_marker:
              ss.file_name_format='SS-GOB_1{}_{}_{}_grid_{}.nc'
        else:
            startdate=datetime.datetime(2006,1,1)
            enddate=datetime.datetime(2058,12,31)
        if 'D' in name_marker:
            ss.file_name_format='SS-GOB_1{}_{}_{}_grid_{}.nc'
            
        datadir = ss.root_data_out+"/derived_data/" #where everyt output is stored
        
        ss.main_data_folder= ss.root_data_in+"/{}{}/".format(folder_start,name_marker)
        depth_ax='deptht'
        if var1 in ['SST','SSS','SSH_inst']:
            ss.grid_type='T'
            ss.interval='d'
        if var1 in ['vosaline','votemper']:
            ss.grid_type='T'
            ss.interval='m'
            depth_ax='deptht'
        
        
        if '1' in name_marker: #the 001 series are hindcasts, all other scenarios
            series_name='Hindcast_{}_{}'.format(name_marker,var1)
        elif 'REANALYSIS' in name_marker:
            series_name='{}_{}'.format(name_marker,var1)
        else:
            series_name='Scenario_{}_{}'.format(name_marker,var1)
        
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
        print(ss.main_data_folder,f)
            
        running_number=0
        
        just_one=False
        if(just_one):
            files_working=[files_working[0]]
        mean_variables=None
        time_axis=None    
        is_first=True
        yearly_variable_data={}
        for f in files_working:
            data=Dataset(ss.main_data_folder+f)
            d=data.variables[var1][:]
            d=np.ma.masked_where(d==0.0,d)
            if(data.variables[var1].units=='degC'):
              d+=273.15
            #dshape should be time,depth,lats,lons
            lons = data.variables['nav_lon'][:]
            lats = data.variables['nav_lat'][:]
            lats,lons=ss.fix_latslons(lats,lons)
            times=data.variables['time_counter'][:].data
            time_units=data.variables['time_counter'].units
            time_calendar=data.variables['time_counter'].calendar
            areas=ss.give_areas(lats,lons)
            if(len(d.shape)==4):
                depths=data.variables['deptht'][:]
                volumes=np.repeat(areas[np.newaxis,:,:],len(depths),axis=0)
                for i in range(len(depths)): #Fix to use depth_bounds
                    if i>0:
                        layer_depth=depths[i]-depths[i-1]
                    else:
                        layer_depth=depths[i]
                    volumes[:,:,i]*=layer_depth*0.001  #because km
            else:
                raise ValueError("data doesn't seem to be 3D") 
            if(just_one):
                times=[times[0]]
            for time_frame in range(len(times)):
                time=ss.nemo_time_to_datetime(times[time_frame])
                if(np.sum(~d[time_frame,:,:].mask)>0): #which means there is some actual number here
                    variable_times_volume=np.sum(d[time_frame,:,:,:]*volumes,axis=(1,2)) #sum over lat,lon
                    mean_variable=variable_times_volume/np.sum(~d[time_frame,:,:,:].mask*volumes,axis=(1,2))
                    if time.year in yearly_variable_data.keys():
                        yearly_variable_data[time.year]['samples']+=1
                        yearly_variable_data[time.year]['mean']+=mean_variable
#                        if(yearly_variable_data[time.year]['min']>mean_variable):
#                            yearly_variable_data[time.year]['min']=mean_variable
#                        if(yearly_variable_data[time.year]['max']<mean_variable):
#                            yearly_variable_data[time.year]['max']=mean_variable
                    else:
                        yearly_variable_data[time.year]={'year':time.year,
#                                                                  'min':mean_variable,
#                                                                  'max':mean_variable,
                                                                  'mean':mean_variable,
                                                                  'samples':1}
                        
                    day_year=time.timetuple().tm_yday
                    month_year=time.month-1
                    if is_first:
                        mean_variables=np.array([variable_times_volume])
                        time_axis=np.array([time])
                        is_first=False
                    else:
                        mean_variables=np.concatenate((mean_variables,[variable_times_volume]))
                        time_axis=np.concatenate((time_axis,[time]))
                running_number+=1
                print("Analysing {} ({} of approx {})".format(time,running_number,int(len(files_working)*30.5)))
            data.close()
        #iris system
        #first define time axis:
        iris_time_unit=cf_units.Unit(time_units,calendar=time_calendar)
        iris_time_coord=iris.coords.DimCoord(iris_time_unit.date2num(time_axis),standard_name='time',units=iris_time_unit)
        iris_depth_coord=iris.coords.DimCoord(depths,standard_name='depth',units='m')

        cube_J=iris.cube.Cube(mean_variables,long_name='thermal_energy',units='J')
        cube_J.add_dim_coord(iris_time_coord,0)
        cube_J.add_dim_coord(iris_depth_coord,1)
        iris.save(cube_J,datadir+'reserve_{}_{}.nc'.format(var1,name_marker))        
#        #fix the averages:
#        for i in yearly_variable_data:
#            yearly_variable_data[i]['mean']/=float(yearly_variable_data[i]['samples'])
#
#        full_data=pd.DataFrame({'time':time_axis,'mean_{}'.format(var1):mean_variables})
#        full_data.to_csv(datadir+'reserve_{}_{}.csv'.format(var1,name_marker),index=False)
#        
#        yearly_data=pd.DataFrame(yearly_variable_data).transpose()
#        yearly_data.to_csv(datadir+'reserve_year_{}_{}.csv'.format(var1,name_marker),index=False)
