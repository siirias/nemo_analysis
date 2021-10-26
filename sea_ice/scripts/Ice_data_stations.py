#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 18 09:54:54 2020

@author: oikkonea
"""
# This script extracts sea ice data (concetration, volume) for coastal stations.
# Daily data is saved in txt files for further analysis and plotting
# Stations: Kemi, Oulu (Saapaskari), Kalajoki, Kylmäpihlaja and Sälgrund


from __future__ import division
import matplotlib.pyplot as plt
import numpy as np
import netCDF4 as nc
import datetime as dt
import sys
sys.path.insert(0,'/projappl/project_2001635/siiriasi/nemo_analysis/')
from smartseahelper import smh
ss=smh()
in_dir = "/scratch/project_2002540/siiriasi/smartsea_data/"
out_dir = "/scratch/project_2002540/siiriasi/smartsea_data/derived_data/ice_data/"
reftime=dt.datetime(1950,1,1)
#Kemi
measuring_points = {
'Kemi':{'lat':65.72, 'lon':24.43},
'Saapaskari':{'lat':65.05, 'lon':25.17},
'Kalajoki':{'lat':64.29, 'lon':23.89},
'Kylmapihlaja':{'lat':61.14, 'lon':21.31},
'Salgrund':{'lat':62.33, 'lon':21.21}
}
simulation_sets = ['A001', 'B001', 'D001', 
                   'A002', 'B002', 'D002', 
                   'A005', 'B005', 'D005']
simulation_sets = ['D002', 
                   'A005', 'B005', 'D005']
simulation_sets = ['REANALYSIS']
for simulation_set in simulation_sets:
    start_year = 2006
    end_year = 2100
    if(simulation_set in ['A001', 'B001', 'D001']):
        start_year = 1976 
        end_year = 2006
    if(simulation_set in ['REANALYSIS']):
        start_year = 1981 
        end_year = 2008
    for point_name in measuring_points.keys():
        obs_lat = measuring_points[point_name]['lat']
        obs_lon = measuring_points[point_name]['lon']
        with open(out_dir + "ice_seasons/Ice_season_{}_{}.txt".\
                    format(point_name, simulation_set), "w") as a_file:
            for i in range(start_year,end_year):  
                year=i
                #data=nc.Dataset('run_data/NORDIC-GOB_1d_{:}0101_{:}1231_grid_T.nc'.format(year,year))

                try: 
                    data=nc.Dataset('{}/{}/NORDIC-GOB_1d_{:}0101_{:}1231_grid_T.nc'.format(in_dir, simulation_set, year,year))
                    print(year)
                    lons = data.variables['nav_lon'][:]
                    lats = data.variables['nav_lat'][:]
                    
                  
                    abslat = np.abs(lats-obs_lat)
                    abslon= np.abs(lons-obs_lon)
                    c = np.maximum(abslon,abslat)
                    latlon_idx = np.argmin(c)
                    x, y = np.where(c == np.min(c)) 
                    vol=data.variables['icevolume'][:,:,:]
                    conc = data.variables['icecon'][:,:,:]
                    
                    vol = vol[:,x,y]
                    conc = conc[:,x,y]

                    for j in range(len(vol)):
                        a_file.write("\n")
                        a_file.write("{:},{:},{:},{:}".format(year,j,conc[j][0],vol[j][0]))
                except FileNotFoundError:
                    print("Year {} file not found".format(year))

