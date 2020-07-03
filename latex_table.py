# -*- coding: utf-8 -*-
import datetime as dt
import numpy as np
import pandas as pd
import os
import re
import netCDF4
from netCDF4 import Dataset
import xarray as xr
in_dir = '/scratch/project_2001635/siiriasi/smartsea_data/A001/'
files = os.listdir(in_dir)
files = [f for f in files if  bool(re.match('2000.*shark.*nc',f))]
names = []
lats = []
lons = []
depths = []
for f in files: 
    name = re.search("shark(.*)\.nc",f).group(1) # sharkpoint name
    d = xr.open_dataset(in_dir + f)
    lat = float(d.nav_lat)
    lon = float(d.nav_lon)
    depth = np.sum(np.array(d.votemper[0,:,0,0] != 0.0)) # assuming temperature exists in and only in valid depths.
    depth = float(d.deptht_bounds[depth-1,1])
    depths.append(depth)
    names.append(name)
    lats.append(lat)
    lons.append(lon)
    print(name)
#    print("{}  &{:.3f}  &{:.3f}  &{:.0f}  \\\\".format(name, lat, lon, depth)) 
formatters = {\
'Name':lambda x:'{:s}'.format(x),\
'Latitude':lambda x:'{:.2f}'.format(x),\
'Longitude':lambda x:'{:.2f}'.format(x),\
'Model depth':lambda x:'{:.2f}'.format(x),\
}

table = pd.DataFrame({  'Name':names, \
                        'Latitude':lats,\
                        'Longitude':lons,\
                        'Model depth':deths})

table = table.sort_values('Latitude',ascending = False)
table_str = table.to_latex(formatters = formatters, index = False)
print(table_str)
