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
from smartseahelper import smh
ss=smh()

reftime=dt.datetime(1950,1,1)

#Kemi
obs_lat = 65.72
obs_lon = 24.43


with open("ice_seasons/Ice_season_Kemi_D005.txt", "w") as a_file:

    for i in range(2006,2060):  
        year=i
        print year
        #data=nc.Dataset('run_data/NORDIC-GOB_1d_{:}0101_{:}1231_grid_T.nc'.format(year,year))

        data=nc.Dataset('run_data/SS-GOB_1d_{:}0101_{:}1231_grid_T.nc'.format(year,year))
        
        lons = data.variables['nav_lon'][:]
        lats = data.variables['nav_lat'][:]
        
      
        abslat = np.abs(lats-obs_lat)
        abslon= np.abs(lons-obs_lon)
        c = np.maximum(abslon,abslat)
        latlon_idx = np.argmin(c)
        x, y = np.where(c == np.min(c)) 
        vol=data.variables['icevolume'][:,:,:]
        conc = data.variables['icecon'][:,:,:]
        
        vol_Kemi=vol[:,x,y]
        conc_Kemi=conc[:,x,y]

        for j in range(len(vol_Kemi)):
            a_file.write("\n")
            a_file.write("{:},{:},{:},{:}".format(year,j,conc_Kemi[j][0],vol_Kemi[j][0]))

            
#saapaskari (oulu)
obs_lat = 65.05
obs_lon = 25.17


with open("ice_seasons/Ice_season_Saapaskari_D005.txt", "w") as a_file:

    for i in range(2006,2060):  
        year=i
        print year
#        data=nc.Dataset('run_data/NORDIC-GOB_1d_{:}0101_{:}1231_grid_T.nc'.format(year,year))

        data=nc.Dataset('run_data/SS-GOB_1d_{:}0101_{:}1231_grid_T.nc'.format(year,year))
        
        lons = data.variables['nav_lon'][:]
        lats = data.variables['nav_lat'][:]
        

        abslat = np.abs(lats-obs_lat)
        abslon= np.abs(lons-obs_lon)
        c = np.maximum(abslon,abslat)
        latlon_idx = np.argmin(c)
        x, y = np.where(c == np.min(c)) 
        vol=data.variables['icevolume'][:,:,:]
        conc = data.variables['icecon'][:,:,:]
        
        vol_Saapaskari=vol[:,x,y]
        conc_Saapaskari=conc[:,x,y]

        for j in range(len(vol_Saapaskari)):
            
            a_file.write("\n")
            a_file.write("{:},{:},{:},{:}".format(year,j,conc_Saapaskari[j][0],vol_Saapaskari[j][0]))


#Kalajoki
obs_lat = 64.29
obs_lon = 23.89


with open("ice_seasons/Ice_season_Kalajoki_D005.txt", "w") as a_file:

    for i in range(2006,2060):  
        year=i
        print year
#        data=nc.Dataset('run_data/NORDIC-GOB_1d_{:}0101_{:}1231_grid_T.nc'.format(year,year))

        data=nc.Dataset('run_data/SS-GOB_1d_{:}0101_{:}1231_grid_T.nc'.format(year,year))

        lons = data.variables['nav_lon'][:]
        lats = data.variables['nav_lat'][:]
        
        abslat = np.abs(lats-obs_lat)
        abslon= np.abs(lons-obs_lon)
        c = np.maximum(abslon,abslat)
        latlon_idx = np.argmin(c)
        x, y = np.where(c == np.min(c)) 
        vol=data.variables['icevolume'][:,:,:]
        conc = data.variables['icecon'][:,:,:]
        
        vol_Kalajoki=vol[:,x,y]
        conc_Kalajoki=conc[:,x,y]

        for j in range(len(vol_Kalajoki)):
            
            a_file.write("\n")
            a_file.write("{:},{:},{:},{:}".format(year,j,conc_Kalajoki[j][0],vol_Kalajoki[j][0]))


#Kylmäpihlaja
obs_lat = 61.14
obs_lon = 21.31


with open("ice_seasons/Ice_season_Kylmapihlaja_D005.txt", "w") as a_file:

    for i in range(2006,2060):  
        year=i
        print year
#        data=nc.Dataset('run_data/NORDIC-GOB_1d_{:}0101_{:}1231_grid_T.nc'.format(year,year))

        data=nc.Dataset('run_data/SS-GOB_1d_{:}0101_{:}1231_grid_T.nc'.format(year,year))
        
        lons = data.variables['nav_lon'][:]
        lats = data.variables['nav_lat'][:]
        
        abslat = np.abs(lats-obs_lat)
        abslon= np.abs(lons-obs_lon)
        c = np.maximum(abslon,abslat)
        latlon_idx = np.argmin(c)
        x, y = np.where(c == np.min(c)) 
        vol=data.variables['icevolume'][:,:,:]
        conc = data.variables['icecon'][:,:,:]
        
        vol_Kylmapihlaja=vol[:,x,y]
        conc_Kylmapihlaja=conc[:,x,y]

        for j in range(len(vol_Kylmapihlaja)):
            
            a_file.write("\n")
            a_file.write("{:},{:},{:},{:}".format(year,j,conc_Kylmapihlaja[j][0],vol_Kylmapihlaja[j][0]))



#Sälgrund
obs_lat = 62.33
obs_lon = 21.21


with open("ice_seasons/Ice_season_Salgrund_D005.txt", "w") as a_file:

    for i in range(1980,20092006,2060):  
        year=i
        print year
#        data=nc.Dataset('run_data/NORDIC-GOB_1d_{:}0101_{:}1231_grid_T.nc'.format(year,year))

        data=nc.Dataset('run_data/SS-GOB_1d_{:}0101_{:}1231_grid_T.nc'.format(year,year))

        lons = data.variables['nav_lon'][:]
        lats = data.variables['nav_lat'][:]
         
        abslat = np.abs(lats-obs_lat)
        abslon= np.abs(lons-obs_lon)
        c = np.maximum(abslon,abslat)
        latlon_idx = np.argmin(c)
        x, y = np.where(c == np.min(c)) 
        vol=data.variables['icevolume'][:,:,:]
        conc = data.variables['icecon'][:,:,:]
        
        vol_Salgrund=vol[:,x,y]
        conc_Salgrund=conc[:,x,y]

        for j in range(len(vol_Salgrund)):
            
            a_file.write("\n")
            a_file.write("{:},{:},{:},{:}".format(year,j,conc_Salgrund[j][0],vol_Salgrund[j][0]))
            
