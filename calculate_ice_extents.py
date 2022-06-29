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
#import matplotlib as mp
#import matplotlib.pyplot as plt
import numpy as np
from scipy.io import netcdf
#from mpl_toolkits.basemap import Basemap
#from mpl_toolkits.basemap import Basemap, shiftgrid, cm
from netCDF4 import Dataset
from smartseahelper import smh
import os
#import cmocean
import pandas as pd

ss=smh()
ss.grid_type='T'
ss.interwall='d'

#name_markers=['A001','A002','A005']
name_markers=['REANALYSIS']
#name_markers=['A001']  #this one tells which dataseries is handled.
for name_marker in name_markers:
    
    if '1' in name_marker: #the 001 series are hindcasts, all other scenarios
        startdate=datetime.datetime(1975,1,1)
        enddate=datetime.datetime(2005,12,31)
    elif 'REANALYSIS' in name_marker:
        startdate=datetime.datetime(1980,1,1)
        enddate=datetime.datetime(2008,12,31)
    else:
        startdate=datetime.datetime(2006,1,1)
        enddate=datetime.datetime(2058,12,31)
    datadir = ss.root_data_out+"/derived_data/" #where everyt output is stored
    ss=smh()
    ss.grid_type='T'
    ss.interwall='d'
    ss.main_data_folder= ss.root_data_in+"/{}/".format(name_marker)
    var1 = 'SST'  #'SST', 'SSS'
    var2 = 'soicecov'    
    
    if '1' in name_marker: #the 001 series are hindcasts, all other scenarios
        series_name='Hindcast_{}_{}'.format(name_marker,var1)
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
        
    running_number=0
    
    just_one=False
    if(just_one):
        files_working=[files_working[0]]
    ice_extents=0.0
    time_axis=None    
    is_first=True
    yearly_ice_maximum={}
    for f in files_working:
        data=Dataset(ss.main_data_folder+f)
        d=data.variables[var1][:]
        d=np.ma.masked_where(d==0.0,d)
        if(var2 is not None):
            ice_d=data.variables[var2][:]
            ice_d=np.ma.masked_where(ice_d<0.2,ice_d)
        lons = data.variables['nav_lon'][:]
        lats = data.variables['nav_lat'][:]
        lats,lons=ss.fix_latslons(lats,lons)
        times=data.variables['time_counter'][:].data
        areas=ss.give_areas(lats,lons)
        if(just_one):
            times=[times[0]]
        for time_frame in range(len(times)):
            time=ss.nemo_time_to_datetime(times[time_frame])
            ice_extent=np.sum(ice_d[time_frame,:,:]*areas)
            winter=time+datetime.timedelta(days=6*30) #winter isat the change of year so...
            if(np.sum(~ice_d[time_frame,:,:].mask)==0): #which means we must put zero as ice
                ice_extent=0.
            if winter.year in yearly_ice_maximum.keys():
                if(ice_extent>yearly_ice_maximum[winter.year]):
                    yearly_ice_maximum[winter.year]=ice_extent
            else:
                yearly_ice_maximum[winter.year]=ice_extent
                
            if is_first:
                ice_extents=np.array([ice_extent])
                time_axis=np.array([time])
                is_first=False
            else:
                ice_extents=np.concatenate((ice_extents,[ice_extent]))
                time_axis=np.concatenate((time_axis,[time]))
            running_number+=1
            print("Analysing {} ({} of approx {})".format(time,running_number,int(len(files_working)*30.5)))
    
    full_data=pd.DataFrame({'time':time_axis,'ice_extent':ice_extents,'ice_percentage':ice_extents/areas[d[0].mask == False].sum()})
    full_data.to_csv(datadir+'ice_extent_{}.csv'.format(name_marker),index=False)
    yearly_data=pd.DataFrame({'year':list(yearly_ice_maximum.keys()),'max_ice':list(yearly_ice_maximum.values())})
    yearly_data.to_csv(datadir+'yearly_max_ice_{}.csv'.format(name_marker),index=False)
