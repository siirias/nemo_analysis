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
csv_header = "Cruise,Station,Type,yyyy-mm-ddThh:mm,Latitude [degrees_north],Longitude [degrees_east],Bot. Depth [m],Secchi Depth [m]:METAVAR:FLOAT,PRES [db],TEMP [deg C],PSAL [psu],DOXY [ml/l],PHOS [umol/l],TPHS [umol/l],SLCA [umol/l],NTRA [umol/l],NTRI [umol/l],AMON [umol/l],NTOT [umol/l],PHPH [],ALKY [meq/l],CPHL [ug/l]"
set_names = ['A001','A002','A005','B001','B002','B005','D002','D005']
for set_name in set_names:
    data_dir= "/scratch/project_2001635/siiriasi/smartsea_data/{}/".format(set_name)
    out_dir = "/scratch/project_2001635/siiriasi/smartsea_data/converted_csv/"
#    data_dir= "D:\\Data\\SmartSeaModeling\\SharkExamples\\"
#    out_dir = "D:\\Data\\SmartSeaModeling\\converted_csv\\"

    csv_format = \
    "Nemo,X,{},{:},{:4.2f},{:4.2f},{:4.2f},,{:4.2f},{:4.2f},{:4.2f},{:4.2f},,,,,,,,,,,\n"
    #csv_header = csv_header.split(',')
    files = os.listdir(data_dir)
    dat={}
    name_format = '(\d*)_shark(.*)\.nc'
    csv_filename_header = '{}_'.format(set_name)
    def date_to_csv_format(dt_arg):
        return dt_arg.strftime("%Y-%m-%dT%H:%M")
    print("gathering data from files..")
    for f in files:
        name_search=re.search(name_format,f)
        if(name_search):
            point_year = name_search.groups()[0]
            point_name = name_search.groups()[1]
            entry_name = point_year+'_'+point_name
            print(".",end="")
            with Dataset(data_dir+f) as D:
                lat = D['nav_lat'][0].data[0]
                lon = D['nav_lon'][0].data[0]
                deptht = D['deptht'][:]
                temperature = D['votemper'][:]
                salinity = D['vosaline'][:]
                oxy = D['oxygen'][:]
                time_stamp = netCDF4.num2date(D['time_counter'][:],\
                                              D['time_counter'].units)
                time_stamp = list(map(date_to_csv_format,time_stamp))
                csv_filename = out_dir+csv_filename_header+\
                            re.search("(.*)\.nc",f).groups()[0]+'.csv'
                dat[entry_name] = {'lat':lat, 'lon':lon,
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
            warnings_on_file = 0
            csv_file.write(csv_header+'\n')
            print(d['csv_filename'])
            profiles = d['votemper'].shape[0]
            depths = d['votemper'].shape[1]
            for profile in range(profiles):
                for depth in range(depths):
                    try:
                        max_depth = d['deptht']\
                                        [(d['votemper'][profile][:] != 0.0)[:,0,0]].max()
                    except:
                        max_depth = 0.0
                        warnings_on_file+=1
                        if(warnings_on_file==1):
                            print("Warning, something wrong with data: {}".format(\
                                                d['csv_filename']))
                    if(d['votemper'][profile][depth] != 0.0): #skip the empty lines
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
            if(warnings_on_file>0):
                print("warning: {}".format(warnings_on_file))
