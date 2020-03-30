# -*- coding: utf-8 -*-
"""
Created on Mon Mar 30 13:17:42 2020

@author: siirias
"""

# -*- coding: utf-8 -*-
"""
Created on Tue Jul  3 14:11:45 2018

@author: siirias
"""

import datetime
import numpy as np
from scipy.io import netcdf
import netCDF4
from netCDF4 import Dataset
from smartseahelper import smh
import os
import re


data_dir= "D:\\Data\\SmartSeaModeling\\SharkFiles\\A001\\"
out_dir = "D:\\Data\\SmartSeaModeling\\converted_csv\\"

csv_header = "Cruise,Station,Type,yyyy-mm-ddThh:mm,Latitude [degrees_north],Longitude [degrees_east],Bot. Depth [m],Secchi Depth [m]:METAVAR:FLOAT,PRES [db],TEMP [deg C],PSAL [psu],DOXY [ml/l],PHOS [umol/l],TPHS [umol/l],SLCA [umol/l],NTRA [umol/l],NTRI [umol/l],AMON [umol/l],NTOT [umol/l],PHPH [],ALKY [meq/l],CPHL [ug/l]"
csv_format = \
"Nemo,X,{},{:},{:4.2f},{:4.2f},{:4.2f},,{:4.2f},{:4.2f},{:4.2f},{:4.2f},,,,,,,,,,,\n"
#csv_header = csv_header.split(',')
files = os.listdir(data_dir)
dat={}
name_format = '\d*_shark(.*)\.nc'
csv_filename_header = 'A001_'
def date_to_csv_format(dt_arg):
    return dt_arg.strftime("%Y-%m-%dT%H:%M")

for f in files:
    point_name=re.search(name_format,f)
    if(point_name):
        point_name = point_name.groups()[0]
        print(point_name)
        with Dataset(data_dir+f) as D:
            lat = D['nav_lat'][0].data[0]
            lon = D['nav_lon'][0].data[0]
            deptht = D['deptht'][:]
            temperature = D['votemper'][:]
            salinity = D['vosaline'][:]
            oxy = D['oxygen'][:]
            time_stamp = netCDF4.num2date(D['time_counter'][:],\
                                          D['time_counter'].units)
            time_stamp = map(date_to_csv_format,time_stamp)                        
            csv_filename = out_dir+csv_filename_header+\
                        re.search("(.*)\.nc",f).groups()[0]+'.csv'
            dat[point_name] = {'lat':lat, 'lon':lon,
                               'date':time_stamp,
                               'deptht':deptht,
                               'votemper':temperature,
                               'vosaline':salinity,
                               'oxygen':oxy,
                               'csv_filename':csv_filename
                               }
for d_ind in dat:
    d = dat[d_ind]
    with open(d['csv_filename'],'w') as csv_file:
        csv_file.write(csv_header+'\n')
        print(d['csv_filename'])
        print(d['deptht'])
        profiles = d['votemper'].shape[0]
        depths = d['votemper'].shape[1]
        for profile in range(profiles):
            for depth in range(depths):
                try:
                    max_depth = d['deptht']\
                                    [(d['votemper'][profile][:] != 0.0)[:,0,0]].max()
                except:
                    max_depth = 0.0
                    print("Warning, something wrong with data")
                csv_file.write(csv_format.format(\
                               profile+1,
                               d['date'][profile],
                               d['lat'],
                               d['lon'],
                               max_depth,
                               d['deptht'][depth],
                               float(d['votemper'][profile][depth]),
                               float(d['vosaline'][profile][depth]),
                               float(d['oxygen'][profile][depth]),
                               ))