# -*- coding: utf-8 -*-
"""
Created on Thu Mar  9 16:10:05 2023

@author: siirias
"""
from smartseahelper import smh
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt

ss=smh()
dat = xr.open_dataset('D:/SmartSea/new_dataset/B001/NORDIC-GOB_1m_19970101_19971231_grid_T.nc')
lons = dat.variables['nav_lon'][:]
lats = dat.variables['nav_lat'][:]
lats,lons=ss.fix_latslons(lats,lons)
areas = ss.give_areas(lats, lons)
depths=dat.variables['deptht'][:]
d_bounds = np.array(dat['deptht_bounds'])
volumes=np.repeat(areas[np.newaxis,:,:],len(depths),axis=0)
for i in range(len(depths)): #Fix to use depth_bounds
    layer_depth = d_bounds[i,1] - d_bounds[i,0]
    volumes[i,:,:]*=layer_depth*0.001  #because km
    
x = np.array(dat.variables['votemper'][0,:,:,:]) #Just something with values
volumes = np.ma.masked_array(volumes, np.isnan(x) + (x == 0.0))

for i in range(len(depths)):
    tmp = volumes[i,:,:]
    l_depth = d_bounds[i,1] - d_bounds[i,0]
    print("layer {}\t volume {:0.1f} km^3\t layer depth {:0.2f} \t area:{:0.2f}".format(
                    i,
                    np.sum(tmp[~tmp.mask]),
                    l_depth,
                    np.sum(tmp[~tmp.mask])/(0.001*l_depth)
                    ))

print("Total volume of the model area: {}".format(
                    np.sum(volumes[~volumes.mask])))



