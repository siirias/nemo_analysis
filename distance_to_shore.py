# -*- coding: utf-8 -*-
"""
Created on Thu Jan 17 17:47:05 2019

@author: siirias
Script to determine routes for closest land-point for every sea point.
"""

import sys

import datetime
import matplotlib as mp
import matplotlib.pyplot as plt
import numpy as np
import random
from netCDF4 import Dataset

data_dir="C:\\Data\\SmartSeaModeling\\"
input_file="bathy_meter.nc"
output_dir = "c:\\Data\\"

nc_data=Dataset(data_dir+input_file,'r')
bathymetry=nc_data['Bathymetry'][:,:]
lats=nc_data['lat'][:]
lons=nc_data['lon'][:]
make_fake_ice_mask=True

bathymetry.data[bathymetry.data==bathymetry.fill_value]=0.0 #get rid of the fillvalues
bathymetry=bathymetry.data
nc_data.close()

#Now lets make the distance thing, start from every non 0 is high number, and work from there.
high_number=10000.0
small_number=0.01 #number small enough to have no meaning when optimizing closest route
x_min=0
y_min=0
x_max=bathymetry.shape[0]
y_max=bathymetry.shape[1]
diag_extra_cost=np.sqrt(2)-1.0 #diagonal contact is this much further
def neighbour_num(dx,dy):
    result=0
    if(dx==1):
        result+=1
    if(dx==-1):
        result+=2
    if(dy==1):
        result+=4
    if(dy==-1):
        result+=8
    return result
def direction_from_neighbour_num(neighbour_num_in):
    dx=0
    dy=0
    if(neighbour_num_in&1):
        dx+=1
    if(neighbour_num_in&2):
        dx-=1
    if(neighbour_num_in&4):
        dy+=1
    if(neighbour_num_in&8):
        dy-=1
    return (dx,dy)

def give_value_and_index(distance_grid,x,y):
    return(distance_grid[x,y],y*x_max+x)
    
def nearest_neighbor(distance_grid,x,y,ignore=high_number):
    best=high_number
    best_index=y*x_max+x
    if(x>x_min):
        (neighbor,neighbor_index)=give_value_and_index(distance_grid,x-1,y)
        if(neighbor<best):
            (best,best_index)=(neighbor,neighbor_index)
    if(y>y_min):
        (neighbor,neighbor_index)=give_value_and_index(distance_grid,x,y-1)
        if(neighbor<best):
            (best,best_index)=(neighbor,neighbor_index)
    if(x<x_max-1):
        (neighbor,neighbor_index)=give_value_and_index(distance_grid,x+1,y)
        if(neighbor<best):
            (best,best_index)=(neighbor,neighbor_index)
    if(y<y_max-1):
        (neighbor,neighbor_index)=give_value_and_index(distance_grid,x,y+1)
        if(neighbor<best):
            (best,best_index)=(neighbor,neighbor_index)
    #diagonals:
    if(x>x_min):
        if(y>y_min):
            (neighbor,neighbor_index)=give_value_and_index(distance_grid,x-1,y-1)
            neighbor+=diag_extra_cost
            if(neighbor<best):
                (best,best_index)=(neighbor,neighbor_index)
        if(y<y_max-1):
            (neighbor,neighbor_index)=give_value_and_index(distance_grid,x-1,y+1)
            neighbor+=diag_extra_cost
            if(neighbor<best):
                (best,best_index)=(neighbor,neighbor_index)
    if(x<x_max-1):
        if(y>y_min):
            (neighbor,neighbor_index)=give_value_and_index(distance_grid,x+1,y-1)
            neighbor+=diag_extra_cost
            if(neighbor<best):
                (best,best_index)=(neighbor,neighbor_index)
        if(y<y_max-1):
            (neighbor,neighbor_index)=give_value_and_index(distance_grid,x+1,y+1)
            neighbor+=diag_extra_cost
            if(neighbor<best):
                (best,best_index)=(neighbor,neighbor_index)
    new_x=best_index%x_max
    new_y=np.floor(best_index/x_max)
    return (best,best_index,neighbour_num(new_x-x,new_y-y))

distance_grid=bathymetry[:,:].copy()
distance_grid[distance_grid>0.0]=high_number
index_grid=np.zeros(distance_grid.shape,dtype=np.int)
    
changes_needed=[]
top_distance=high_number
while top_distance==high_number or len(changes_needed)>0:
    changes_needed=[]
    for x in range(x_min,x_max):
            for y in range(y_min,y_max):
                if(distance_grid[x,y]>0.0+small_number): #don't touch the land
                    (dist,index,neighbor)=nearest_neighbor(distance_grid,x,y)
                    if(dist<1.0): #This is one next to land then
                        neighbor=0
                    if(dist+1+small_number<distance_grid[x,y]):
                        changes_needed.append((x,y,dist+1,index,neighbor))
    #update the grid as we have went it through:
    for i in changes_needed:
        distance_grid[i[0],i[1]]=i[2]
        last=index_grid[i[0],i[1]]
        index_grid[i[0],i[1]]=i[4]
    top_distance=distance_grid.max()
    print(distance_grid[distance_grid<high_number].max(), top_distance, len(changes_needed))


#Plot what we got
distance_copy=distance_grid[:,:].copy()
distance_copy[distance_copy<small_number]=-15.0   
plt.imshow(distance_copy[-1:0:-1,:])
#plt.figure()
#plt.imshow(bathymetry[-1:0:-1,:])

#Let's write out the indices netcdf
index_netcdf=Dataset(output_dir+'neighbours.nc','w',format='NETCDF4_CLASSIC')
lat=index_netcdf.createDimension('lat',index_grid.shape[0])
lon=index_netcdf.createDimension('lon',index_grid.shape[1])
index_netcdf.createVariable('lat',np.float64,('lat',))
index_netcdf.createVariable('lon',np.float64,('lon',))

index_netcdf.createVariable('neighbour',np.int16,('lat','lon'))
index_netcdf['neighbour'][:,:]=index_grid

#for debug reasons, and possible future uses, lets store the distance (for now based on straight and diagonal movements) too
index_netcdf.createVariable('distance',np.float64,('lat','lon'))
index_netcdf['distance'][:,:]=distance_grid
index_netcdf['distance'].unit='distance in cell grids'
index_netcdf['lat'][:]=lats
index_netcdf['lon'][:]=lons
index_netcdf.description="Neighbour tells which direction next cell towards nearest land resides.\n\
                          It is a sum of the directions:\n\
                          North 8\n\
                          South 4\n\
                          East  2\n\
                          West  1\n\
                          (such NorthEast=10)\n\
                          Value 0 indicates either land or that this is itself the closest to land"
index_netcdf.close()

if(make_fake_ice_mask):
    fast_ice_mask=bathymetry[:,:].copy()
    fast_ice_mask[fast_ice_mask<15.]=1
    fast_ice_mask[fast_ice_mask>=15.]=0
    ice_netcdf=Dataset(output_dir+'fast_ice_mask.nc','w',format='NETCDF4_CLASSIC')
    lat=ice_netcdf.createDimension('lat',index_grid.shape[0])
    lon=ice_netcdf.createDimension('lon',index_grid.shape[1])
    ice_netcdf.createVariable('lat',np.float64,('lat',))
    ice_netcdf.createVariable('lon',np.float64,('lon',))
    index_netcdf['lat'][:]=lats
    index_netcdf['lon'][:]=lons
    
    ice_netcdf.createVariable('fastice_mask',np.int32,('lat','lon'))
    ice_netcdf['fastice_mask'][:,:]=fast_ice_mask
    ice_netcdf.close()

