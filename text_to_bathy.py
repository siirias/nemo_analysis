"""
A simple script to write out nemo bathymetry nc from an xml text-file.
Suitable xml files can be created with Bathymetry Editor BatMan.
Simo Siiria, FMI 2019
"""

import imp
import iris
import gsw
import siri_omen.utility
import numpy as np
import siri_omen.nemo_reader as nrd

data_folder = '/arch/smartsea/analysis/run_configuration/NAMELISTS/'
output_folder = '/arch/smartsea/analysis/derived_data/'
file_name = 'bathy_meter.nc'
out_file_name = 'bathy_meter_nemo_SS.txt'

data = iris.load(data_folder+file_name,'Bathymetry')[0]
bathy_data = data.data
lat_min = np.min(data.coord('latitude').points)
lat_max = np.max(data.coord('latitude').points)
lon_min = np.min(data.coord('longitude').points)
lon_max = np.max(data.coord('longitude').points)

out_file = open(output_folder+out_file_name,'w')
out_file.write(\
"""\
<?xml version="1.0"?>
<depth_negative dn="0" />
<coordinates lat_min="{}" lat_max="{}" lon_min="{}" lon_max="{}" />
<Data>
""".format(lat_min,lat_max,lon_min,lon_max))

point = [0,0]
for lon in range(bathy_data.shape[data.coord_dims('longitude')[0]]):
    line_to_write = ""
    for lat in range(bathy_data.shape[data.coord_dims('latitude')[0]]):
        point[data.coord_dims('longitude')[0]]=lon
        point[data.coord_dims('latitude')[0]]=lat
        line_to_write+="{} ".format(bathy_data[tuple(point)])
    out_file.write(line_to_write[:-1])  # to get rid of the last space
    out_file.write("\n")

out_file.write(\
"""
</Data>
""")
out_file.close()
