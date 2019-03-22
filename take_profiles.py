# -*- coding: utf-8 -*-
"""
Created on Mon Jul  9 10:16:05 2018

@author: siirias
"""

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

ss=smh()

name_markers=['A001']
variables=['votemper','vosaline','SSH_inst','SST','SSS']
dataset_to_substract = ''    
#profiles=[{'name':'koe','lat':65.0,'lon':23.0},
#          {'name':'SR5','lat':61.0520,'lon':19.3555}
#          ]
profiles=open("monitoring_points.txt",'r').readlines()
profiles=map(lambda x:{'name':x.split('\t')[0],'lat':float(x.split('\t')[1]),'lon':float(x.split('\t')[2])},profiles)

for name_marker in name_markers:
    for var1 in variables:
        datadir = ss.root_data_out+"/derived_data/" #where everyt output is stored
        
        ss.main_data_folder= ss.root_data_in+"/OUTPUT{}/".format(name_marker)

        if '1' in name_marker: #the 001 series are hindcasts, all other scenarios
            startdate=datetime.datetime(1975,1,1)
            enddate=datetime.datetime(2005,12,31)
        else:
            startdate=datetime.datetime(2006,1,1)
            enddate=datetime.datetime(2058,12,31)
        
        
        running_number=0
        variable_name=var1
        depth_ax='deptht'
        if var1 in ['SST','SSS','SSH_inst']:
            ss.grid_type='T'
            ss.interval='d'
        if var1 in ['vosaline','votemper']:
            ss.grid_type='T'
            ss.interval='m'
            depth_ax='deptht'
        
        
        
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
            
        
        
        
        just_one=False
        if(just_one):
            files_working=[files_working[0]]
        mean_temperatures=None
        time_axis=None    
        is_first=True
        
        profile_data={}
        for profile in profiles:
            profile_data[profile['name']]=nc4.Dataset(datadir+'profile_{}_{}_{}.nc'.format(profile['name'],name_marker,var1),'w',format='NETCDF4')
        prof_lat=None
        prof_lon=None
        prof_day=None
        prof_variable=None
        for f in files_working:
            data=Dataset(ss.main_data_folder+f)
            d=data.variables[var1][:]
            d=np.ma.masked_where(d==0.0,d)
            lons = data.variables['nav_lon'][:]
            lats = data.variables['nav_lat'][:]
            lats,lons=ss.fix_latslons(lats,lons)
            times=data.variables['time_counter'][:].data
            if ( depth_ax in data.dimensions.keys() ):
                number_of_depth_points=data.dimensions[depth_ax].size
                depths=data.variables[depth_ax][:]
                print("depth data")
            else:
                print("Warning: no depth data")
                number_of_depth_points=1
                depths=[0.0]
                
                    
            if(just_one):
                times=[times[0]]
            for time_frame in range(len(times)):
                time=ss.nemo_time_to_datetime(times[time_frame])
                variable_data=np.sum(d[time_frame,:,:])
                for profile in profiles:
                    prof_lat_index,prof_lon_index = ss.latlon_index(profile['lat'],profile['lon'],lats,lons)
                    if(number_of_depth_points>1):
                        actual_data=d[time_frame,:,prof_lat_index,prof_lon_index]
                    else:
                        actual_data=d[time_frame,prof_lat_index,prof_lon_index]
                    if is_first:
                        profile_data[profile['name']].createDimension(depth_ax,number_of_depth_points)
                        profile_data[profile['name']].createDimension('time',None)
                        prof_day=profile_data[profile['name']].createVariable('date','i4','time')
                        prof_depth=profile_data[profile['name']].createVariable(depth_ax,'f4',depth_ax)
                        prof_variable=profile_data[profile['name']].createVariable(variable_name,'f4',('time',depth_ax))
                        file_latitude=profile_data[profile['name']].createVariable('latitude','f4')
                        file_latitude[:]=profile['lat']
                        file_longitude=profile_data[profile['name']].createVariable('longitude','f4')
                        file_longitude[:]=profile['lon']
                        prof_variable[running_number,:]=actual_data
                        prof_day[running_number]=times[time_frame]
                        prof_day.units = data.variables['time_counter'].units
                        prof_day.time_origin = data.variables['time_counter'].time_origin
                        prof_day.calendar = data.variables['time_counter'].calendar
                        prof_variable.units = data.variables[var1].units
                        
                        prof_depth[:]=depths
                    else:
                        prof_variable=profile_data[profile['name']].variables[variable_name]
                        prof_day=profile_data[profile['name']].variables['date']
                        prof_variable[running_number,:]=actual_data
                        prof_day[running_number]=times[time_frame]
                is_first=False
                        
                running_number+=1
                print("Analysing {} ({} of approx {})".format(time,running_number,int(len(files_working)*30.5)))
            data.close()
        
        for profile in profile_data:    
            profile_data[profile].close()
        
