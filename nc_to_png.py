"""
Script for saving a png image from a slice of netcdf
Simo Siiria, FMI 2019
"""

import iris
import gsw
import siri_omen.utility
import numpy as np
import siri_omen.nemo_reader as nrd
import PIL
from matplotlib import cm

x_axis = 0 
y_axis = 1 
time_axis = np.nan
#data_folder = '/arch/smartsea/analysis/tmp/'
#file_name = 'iowtopo2_rev03.nc'
#input_field = 7
 
#data_folder = '/arch/smartsea/analysis/run_configuration/NAMELISTS/'
#file_name = 'bathy_meter.nc'
#input_field = 'Bathymetry'


#data_folder = '/arch/smartsea/analysis/A002/'
data_folder = '/arch/smartsea/analysis/tmp/'
#file_name = 'border_slice_potential_temperature_A001_depthlat.nc'
#file_name = 'NORDIC-GOB_1d_20060101_20061231_grid_T.nc'
file_name = 'A002JustSSS.nc'
input_field = 0
x_axis = 1
y_axis = 2
time_axis = 0

output_folder = '/arch/smartsea/analysis/tmp/video_tmp/'
out_file_base = 'SSS_A002'
out_file_format = '.png'
out_file_quality = 100  # va value of 1 - 100, 100 being best
the_colormap = cm.gray
auto_boundaries=False
if (not auto_boundaries):
    d_min = 0.0
    d_max = 8.5  # These valueas are used, unless auto_boundaries
                   # is set, in which case the actual min/max values
                   # of the data will be used.
print('Loading field {}'.format(input_field))
if(type(input_field) == int):  # we want the n:th field found.
    data = iris.load(data_folder+file_name)[input_field]
else:                          # we want a field with specific name
    data = iris.load(data_folder+file_name,input_field)[0]
print('loaded {}'.format(data.name()))
bathy_data = data.data
if( auto_boundaries):
    d_max = bathy_data.max()
    d_min = bathy_data.min()
print("Boundaries: min: {:.2f} MAX: {:.2f}".format(d_min, d_max))
if (not d_max == d_min):  # normalize the data.
    bathy_data = bathy_data - d_min
    bathy_data = bathy_data / (d_max - d_min)
else:
    print("Warning: All values identical!")
    bathy_data[:] = 0.5
if np.isnan(time_axis):
    loops = [0]  # no time axis, just one image
else:
    loops = range(bathy_data.shape[time_axis])
for frame in loops:
    # lets slice the needed part, providing there are any preferences:
    indices = [slice(1)]*len(bathy_data.shape)
    if len(indices)>x_axis:
        indices[x_axis]=slice(None)
    if len(indices)>y_axis:
        indices[y_axis]=slice(None)
    if len(indices)>time_axis:
        indices[time_axis]=slice(frame,frame+1)
    bathy_data2 = np.squeeze(bathy_data[tuple(indices)])
    
    image_data = the_colormap((bathy_data2.data))
    
    #mark masked areas as transparent
    tmp_alpha = image_data[:,:,3]
    tmp_alpha[bathy_data2.mask] = 0.
    image_data[:,:,3] = tmp_alpha
    if(len(loops)>1): #we need numbers
        out_file_name  = out_file_base + "{:04d}".format(frame) + out_file_format
    else:
        out_file_name  = out_file_base + out_file_format

    the_image = PIL.Image.fromarray(np.uint8(image_data*255))
    if(out_file_format in ['.jpg', '.JPG']): #  jpg's can't hande alpha
        the_image = the_image.convert('RGB')
    the_image.save(output_folder + out_file_name, quality = out_file_quality)
    print('Saved frame {}'.format(frame))
