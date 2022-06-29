# -*- coding: utf-8 -*-
"""
Created on Tue May 15 13:26:16 2018

@author: siirias
"""

# -*- coding: utf-8 -*-
"""
Created on Mon Apr 11 18:11:08 2016

@author: siirias
"""
import sys

import matplotlib as mp
import matplotlib.pyplot as plt
import numpy as np
from scipy.io import netcdf
from mpl_toolkits.basemap import Basemap
from mpl_toolkits.basemap import Basemap, shiftgrid, cm
from netCDF4 import Dataset

#Full baltic
lon_min=9;lat_min=53.5;lon_max=30.3;lat_max=66; #26 GoB
#Gotlands deep
#lon_min=18;lat_min=55;lon_max=21;lat_max=59;

plot_bathymetry=False
plot_salt=True
plot_legends=True
plot_routes=True
DataDirectory=".\\"
fig=plt.figure(figsize=(13,13))
#fig=plt.figure(figsize=(10,10))


plt.clf()
bmap = Basemap(llcrnrlon=lon_min,llcrnrlat=lat_min,urcrnrlon=lon_max,urcrnrlat=lat_max, \
resolution = 'i',fix_aspect=False)

bmap.drawcoastlines(linewidth=0.2)
bmap.fillcontinents()
bmap.drawparallels(np.arange(50.,69,2.),labels=[1,0,0,0],linewidth=0,dashes=[5,10])
bmap.drawmeridians(np.arange(12.,30,2.),labels=[0,0,0,1],linewidth=0,dashes=[5,10])

salt_min=0
salt_max=12
dcolormap='jet'
levs=[0,1,2,3,4,5,6,7,8,9,10,11,12]

topodata = Dataset(DataDirectory+'bathy_meter.nc')

topo = topodata.variables['Bathymetry'][:,:]
lons = topodata.variables['lon'][:]
lats = topodata.variables['lat'][:]
lons=lons+(lons[1]-lons[0])
lats=lats+(lats[1]-lats[0])
x=np.tile(lons,(lats.shape[0],1))
y=np.tile(lats,(lons.shape[0],1)).T
bmap.pcolor(x,y,topo,cmap=dcolormap)

dim=topo.shape
string=""
for x in range(dim[1]):
    for y in range(dim[0]):
        string+="{} ".format(topo[y,x])
    string+="\n"
topofile=open("nemotopo.txt",'w')
topofile.write(string)
topofile.close()