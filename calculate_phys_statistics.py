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
import netCDF4 as nc4
from netCDF4 import Dataset
from smartseahelper import smh
import os
import cmocean
import pandas as pd

ss=smh()
ss.grid_type='T'
ss.interval='d'
ss.save_interval='year'
ss.file_name_format="NORDIC-GOB_1{}_{}_{}_grid_{}.nc"  
just_bottom = False  # used to get bottom values from 3d grid
#folder_start='OUTPUT'
#name_markers=['A001','A002']
folder_start=''
#name_markers=['REANALYSIS','A001','B001','C001']
#name_markers=['A001','B001','C001', 'D001']
name_markers=['A005','A002','A001','B005','B002','B001','D001','D002','D005']
#name_markers=['B001','B005','B002','D001','D005','D002']
#name_markers=['D001','D005','D002']
#variables = ['SSS', 'SST', 'SSH_inst']
#variables= ['icecon','icevolume','SSS','SST','SBS','vosaline','votemper']
variables= ['SSS','SST','SBS','SBT','vosaline','votemper']
#variables= ['SBT']
#variables= ['vosaline', 'votemper']
#variables= ['SBS','SSS']
make_climatology = True 
extra_definition = 'c30v_'  # c30v_ or c70v_ 
#extra_definition = ''
#name_marker='A001'  #this one tells which dataseries is handled.
climatology_time_slots=366
if(ss.interval=='m'):
    climatology_time_slots=12
    
for name_marker in name_markers:
    for var1 in variables:
        var1_actual = var1  # some variables load actually a different name        
        if '1' in name_marker: #the 001 series are hindcasts, all other scenarios
            startdate=datetime.datetime(1975,1,1)
            enddate=datetime.datetime(2005,12,31)
        elif 'REANALYSIS' in name_marker:
            startdate=datetime.datetime(1980,1,1)
            enddate=datetime.datetime(2012,12,31)
            ss.save_interval='year'
        else:
            if(extra_definition == 'c70v_'):
                startdate=datetime.datetime(2070,1,1)  # should be 2006  #2028 for earlier
                enddate=datetime.datetime(2099,12,31)
            else:
                startdate=datetime.datetime(2030,1,1)  # should be 2006  #2028 for earlier
                enddate=datetime.datetime(2059,12,31)
        datadir = ss.root_data_out+"/derived_data/test/" #where everyt output is stored
        
        ss.main_data_folder= ss.root_data_in+"/{}{}/".format(folder_start,name_marker)
        depth_ax='deptht'
        if var1 in ['SST','SSS','SSH_inst','SBT','icecon','icevolume']:
            ss.grid_type='T'
            ss.interval='d'
            climatology_time_slots=366
        if var1 in ['vosaline','votemper']:
            ss.grid_type='T'
            ss.interval='m'
            climatology_time_slots=12
            depth_ax='deptht'
        
        
        if '1' in name_marker: #the 001 series are hindcasts, all other scenarios
            series_name='Hindcast_{}_{}'.format(name_marker,var1)
        elif 'REANALYSIS' in name_marker:
            series_name='{}_{}'.format(name_marker,var1)
        else:
            series_name='Scenario_{}_{}'.format(name_marker,var1)
        
        print("Analysing files from {}".format(ss.main_data_folder))
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
            
        running_number=0
        
        just_one=False
        if(just_one):
            files_working=[files_working[0]]
        mean_variables=None
        time_axis=None    
        is_first=True
        yearly_variable_data={}
        if(make_climatology):
            climatology_data=nc4.Dataset(datadir+'{}climatology_{}_{}_{}.nc'.format(\
                    extra_definition,name_marker,ss.interval,var1),'w',format='NETCDF4')
            clim_lat=None
            clim_lon=None
            clim_day=None
            clim_variable=None
        for f in files_working:
            data=Dataset(ss.main_data_folder+f)
            d=data.variables[var1_actual][:]
            if(just_bottom): # this is 3d value, where bottom is wanted
                d = ss.give_bottom_values(d)
            d=np.ma.masked_where(d==0.0,d)
            coord_extra = '' # if file has several gfrids, the coordinates files have the _grid_T
            if('nav_lon_grid_T' in data.variables.keys()):
                coord_extra = '_grid_T'
            lons = data.variables['nav_lon'+coord_extra][:]
            lats = data.variables['nav_lat'+coord_extra][:]
            lats,lons=ss.fix_latslons(lats,lons)
            times=data.variables['time_counter'][:].data
            areas=ss.give_areas(lats,lons)
            if(len(d.shape)==4):
                data_3d=True
                depths=data.variables['deptht'][:]
                volumes=np.repeat(areas[np.newaxis,:,:],len(depths),axis=0)
                for i in range(len(depths)):
                    if i>0:
                        layer_depth=depths[i]-depths[i-1]
                    else:
                        layer_depth=depths[i]
                    volumes[:,:,i]*=layer_depth*0.001  #because km
                
                
            else:
                data_3d=False
            if(just_one):
                times=[times[0]]
            for time_frame in range(len(times)):
                time=ss.nemo_time_to_datetime(times[time_frame])
                if(np.sum(~d[time_frame,:,:].mask)>0): #which means there is some actual number here
                    if(data_3d):
                        variable_times_volume=np.sum(d[time_frame,:,:,:]*volumes)
                        mean_variable=variable_times_volume/np.sum(~d[time_frame,:,:,:].mask*volumes)
                    else:
                        variable_times_area=np.sum(d[time_frame,:,:]*areas)
                        mean_variable=variable_times_area/np.sum(~d[time_frame,:,:].mask*areas)
                    if time.year in yearly_variable_data.keys():
                        yearly_variable_data[time.year]['samples']+=1
                        yearly_variable_data[time.year]['mean']+=mean_variable
                        if(yearly_variable_data[time.year]['min']>mean_variable):
                            yearly_variable_data[time.year]['min']=mean_variable
                        if(yearly_variable_data[time.year]['max']<mean_variable):
                            yearly_variable_data[time.year]['max']=mean_variable
                    else:
                        yearly_variable_data[time.year]={'year':time.year,
                                                                  'min':mean_variable,
                                                                  'max':mean_variable,
                                                                  'mean':mean_variable,
                                                                  'samples':1}
                        
                    day_year=time.timetuple().tm_yday
                    month_year=time.month-1
                    if is_first:
                        if(make_climatology):
                            climatology_data.createDimension('x',data.dimensions['x'+coord_extra].size)
                            climatology_data.createDimension('y',data.dimensions['y'+coord_extra].size)
                            climatology_data.createDimension('time',climatology_time_slots)
                            if(data_3d):
                                climatology_data.createDimension('deptht',data.dimensions['deptht'].size)
                                variable_axis=('time','deptht','y','x')
                            else:
                                variable_axis=('time','y','x')
                            clim_lat=climatology_data.createVariable('nav_lat','f4',('y','x'))
                            clim_lon=climatology_data.createVariable('nav_lon','f4',('y','x'))
                            clim_samples=climatology_data.createVariable('sample_number','i4','time')
                            if(data_3d):
                                clim_day=climatology_data.createVariable('month_of_year','i4','time')
                            else:
                                clim_day=climatology_data.createVariable('day_of_year','i4','time')
                            clim_day[:]=np.arange(1,climatology_time_slots+1)
                            clim_samples[:]=np.zeros((climatology_time_slots,))
                            clim_variable_mean=climatology_data.createVariable('{}_mean'.format(var1),'f4',variable_axis)
                            clim_lat[:]=lats
                            clim_lon[:]=lons
                            clim_variable_mean[:]=0.0
                            if(ss.interval=='m' and data_3d):
                                clim_variable_mean[month_year,:,:,:]+=d[time_frame,:,:,:] #-1 to start from day 0
                                clim_samples[month_year]+=1
                            if(ss.interval=='m' and not data_3d):
                                clim_variable_mean[month_year,:,:]+=d[time_frame,:,:] #-1 to start from day 0
                                clim_samples[month_year]+=1
                            if(ss.interval=='d' and not data_3d):
                                clim_variable_mean[day_year-1,:,:]+=d[time_frame,:,:] #-1 to start from day 0
                                clim_samples[day_year-1]+=1
                            if(ss.interval=='d' and data_3d):
                                clim_variable_mean[day_year-1,:,:,:]+=d[time_frame,:,:,:] #-1 to start from day 0
                                clim_samples[day_year-1]+=1
                            clim_lon.units = 'degrees east'
                            clim_lat.units = 'degrees north'
                            clim_day.units = 'day of year'
                            clim_variable_mean.units = data.variables[var1_actual].units
                            clim_variable_mean.coordinates = "time_centered nav_lon nav_lat"                
                        
                        mean_variables=np.array([mean_variable])
                        time_axis=np.array([time])
                        is_first=False
                    else:
                        mean_variables=np.concatenate((mean_variables,[mean_variable]))
                        time_axis=np.concatenate((time_axis,[time]))
                        if(make_climatology):
                            if(ss.interval=='m' and data_3d):
                                clim_variable_mean[month_year,:,:,:]+=d[time_frame,:,:,:] #-1 to start from day 0
                                clim_samples[month_year]+=1
                            if(ss.interval=='m' and not data_3d):
                                clim_variable_mean[month_year,:,:]+=d[time_frame,:,:] #-1 to start from day 0
                                clim_samples[month_year]+=1
                            if(ss.interval=='d' and not data_3d):
                                clim_variable_mean[day_year-1,:,:]+=d[time_frame,:,:] #-1 to start from day 0
                                clim_samples[day_year-1]+=1
                            if(ss.interval=='d' and data_3d):
                                clim_variable_mean[day_year-1,:,:,:]+=d[time_frame,:,:,:] #-1 to start from day 0
                                clim_samples[day_year-1]+=1
#                                clim_variable_mean[day_year-1,:,:]+=d[time_frame,:,:] #-1 to start from day 0
#                                clim_samples[day_year-1]+=1
                                
                running_number+=1
                print("Analysing {} ({} of approx {})".format(time,running_number,int(len(files_working)*30.5)))
            data.close()
            
        #fix the averages:
        for i in yearly_variable_data:
            yearly_variable_data[i]['mean']/=float(yearly_variable_data[i]['samples'])
        if( make_climatology):
            for i in range(climatology_time_slots):
                if(clim_samples[i]>1):
                    clim_variable_mean[i,:,:]=clim_variable_mean[i,:,:]/float(clim_samples[i])
            climatology_data.close()
        
        full_data=pd.DataFrame({'time':time_axis,'mean_{}'.format(var1):mean_variables})
        full_data.to_csv(datadir+'mean_{}_{}.csv'.format(var1,name_marker),index=False)
        
        yearly_data=pd.DataFrame(yearly_variable_data).transpose()
        yearly_data.to_csv(datadir+'yearly_{}_{}.csv'.format(var1,name_marker),index=False)
