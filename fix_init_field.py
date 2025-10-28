"""
Simo Siiria, FMI 2019
"""

import imp
import iris
import gsw
import matplotlib as mp
import matplotlib.pyplot as plt
import siri_omen.utility
import numpy as np
import siri_omen.nemo_reader as nrd
import iris.plot as iplt
import iris.quickplot as iqplt

data_folder = '/arch/smartsea/analysis/run_configuration/NAMELISTS/'
output_folder = '/arch/smartsea/analysis/derived_data/'
file_name = 'ini_GB_TS.nc'
out_file_name = 'ini_GB_TS_fixed.nc'

def fix_point(point,data_mask,test_data):
    did_change = False
    if(data_mask[point]>0.0):  # not masked, no operation
        return (test_data,False)
    else:
        surroundings=[]
        for i,max_val in zip(point,test_data.shape):
            surroundings.append(slice(i-1,i+2))
        surroundings = tuple(surroundings)
        if(np.sum(data_mask[surroundings])>0.):
            old_val = test_data[point]
            new_val = np.sum(test_data[surroundings]*data_mask[surroundings])/\
                        np.sum(data_mask[surroundings])
            test_data[point] = new_val 
            did_change = True
    return (test_data,did_change)

def fix_mask(test_data,border=1):
    data_mask=np.ones(data_size)
    for i in range(data_size[0]):
        for j in range(data_size[1]):
            if(test_data[i,j]>=fill_number):
                data_mask[i-border:i+border+1,j-border:j+border+1]=0.0
    return data_mask
def fix_layer(test_data):
    border_width=2
    data_mask = fix_mask(test_data,border_width)    
    plt.figure()
    plt.imshow(test_data);plt.colorbar();
    plt.draw()
    print(np.sum(data_mask))
    for iteration in range(10):
        changed = np.zeros(data_size)
        for i in range(data_size[0]):
            for j in range(data_size[1]):
                (test_data,ch) = fix_point((i,j),data_mask,test_data)
                if(ch):
                    changed[i,j] = 1.0
        data_mask += changed 
        print("iteration",iteration,np.sum(data_mask))
    test_data[0,:]=test_data[1,:].copy()
    return test_data
data = iris.load(data_folder+file_name,'vosaline')[0]
data.data[data.data<-1000]=5  # a gludge to get rid of the broken values
for layer in range(data.shape[0]):
    test_data = data.data[layer,:,:]
    test_data[0,:]=test_data[2,:].copy()
    test_data[1,:]=test_data[2,:].copy()
    fill_number = 3.0
    lat_min = np.min(data.coord('latitude').points)
    lat_max = np.max(data.coord('latitude').points)
    lon_min = np.min(data.coord('longitude').points)
    lon_max = np.max(data.coord('longitude').points)
    data_size=test_data.shape
    data_mask=np.ones(data_size)
    test_data = fix_layer(test_data)
    plt.figure()
    plt.imshow(test_data);plt.colorbar();plt.draw()

print("ready!")

plt.show()
