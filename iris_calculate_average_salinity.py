# !/usr/bin/env ipython
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  3 14:11:45 2018

@author: siirias
"""

import datetime
import matplotlib as mp
import matplotlib.pyplot as plt
import numpy as np
from scipy.io import netcdf
from smartseahelper import smh
import os
import cmocean
import pandas as pd
import iris
import iris.plot as iplt
import iris.quickplot as iqplt
import siri_omen
import siri_omen.utility as sou
import siri_omen.nemo_reader as nrd
import cf_units
import iris.util
import gsw  # TEOS-10
import warnings

ss = smh()
ss.grid_type = 'T'
ss.interval = 'm'
#ss.root_data_in = "/lustre/tmp/siirias/o/tmp/"  # gludge as the main disk is not sure enough.
ss.root_data_in = "/arch/smartsea/analysis/"
# folder_start = 'OUTPUT'
#name_markers = ['new_REANALYSIS']
name_markers = ['SSR_6']
#name_markers = ['REANALYSIS_ERA40_new_ice']
#name_markers = ['new_REANALYSIS','REANALYSIS_SMHI','REANALYSIS']
variable_temperature = 'potential_temperature'
variable_salinity = 'salinity'
analysis_type = 'average'  # this is just for naming the output
#collapse_style={'name':'depth','coords':['longitude', 'latitude']}    
collapse_style={'name':'depthlat','coords':['longitude']}    
#collapse_style={'name':'depthlon','coords':['latitude']}    
#collapse_style={'name':'depthlatlon','coords':[]}    
#collapse_style={'name':'total','coords':['longitude','latitude','depth']}    
#domain_limits={'lat_min':59.,'lat_max':60.,'lon_min':-180.0,'lon_max':180.}
domain_limits={'lat_min':49.,'lat_max':80.,'lon_min':-180.0,'lon_max':180.}
def cut_to_limits(cube_in):
    cube_out  = cube_in.extract(iris.Constraint(
                coord_values = \
                {'latitude':lambda l: \
                    domain_limits['lat_min']<l\
                    <domain_limits['lat_max'],\
                 'longitude': lambda l: \
                    domain_limits['lon_min']<l\
                    <domain_limits['lon_max']\
                }))
    return cube_out


for name_marker in name_markers:
    folder_start = 'OUTPUT'
    ss.save_interval = 'year'
    if '1' in name_marker: # the 001 series are hindcasts, all other scenarios
        startdate = datetime.datetime(1975, 1, 1)
        enddate = datetime.datetime(2005, 12, 31)
    elif 'SSR' in name_marker:
        startdate = datetime.datetime(1980, 1, 1)
        enddate = datetime.datetime(2012, 12, 31)
        ss.save_interval = 'year'
        folder_start = ''
        ss.file_name_format = 'NORDIC-GoB_1{}_{}_{}_grid_{}.nc'
        if 'new_' in name_marker:
            ss.file_name_format = 'SS-GOB_1{}_{}_{}_grid_{}.nc'
            ss.save_interval = 'year'
        elif '_SMHI' in name_marker:
            ss.file_name_format = 'NORDIC-GOB_1{}_{}_{}_grid_{}.nc'
            ss.save_interval = 'month'
        elif 'SSR' in name_marker:
            ss.file_name_format = 'NORDIC-GOB_1{}_{}_{}_grid_{}.nc'
            ss.save_interval = 'month'
            
    else:
        startdate = datetime.datetime(2006, 1, 1)
        enddate = datetime.datetime(2058, 12, 31)
    datadir = ss.root_data_out+"/derived_data/" # where everyt output is stored
    
    ss.main_data_folder= ss.root_data_in+"/{}{}/".format(folder_start, name_marker)
    depth_ax = 'deptht'
    ss.grid_type = 'T'
    ss.interval = 'm'
    depth_ax = 'deptht'
    
    
    if '1' in name_marker: # the 001 series are hindcasts, all other scenarios
        series_name = 'Hindcast_{}_{}'.format(name_marker, variable_temperature)
    elif 'REANALYSIS' in name_marker:
        series_name = '{}_{}'.format(name_marker, variable_temperature)
    else:
        series_name = 'Scenario_{}_{}'.format(name_marker, variable_temperature)
    
    filenames = ss.filenames_between(startdate, enddate)
    ok_files = 0
    files_working = []
    for f in filenames:
        if(os.path.isfile(ss.main_data_folder+f)):
            ok_files += 1
            files_working.append(f)
        else:
            print(f)
    print()
    print("ok {} out of {}".format(ok_files, len(filenames)))
    print(ss.main_data_folder, f)
        
    running_number = 0
    
    just_one = False
    if(just_one):
        files_working = [files_working[0]]
    mean_variables = None
    time_axis = None    
    is_first = True
    yearly_variable_data = {}
    iris_list = iris.cube.CubeList([])
    for num, f in enumerate(files_working):
        temperature = None
        salinity = None
        with warnings.catch_warnings():
            # this to not show warnings of wrong units
            warnings.simplefilter("ignore")  
            variables = iris.load(ss.main_data_folder+f,\
                         [variable_temperature, variable_salinity])
            for v in variables:
                if v.name() == variable_temperature:
                    temperature = v
                if v.name() == variable_salinity:
                    salinity = v
            temperature = nrd.remove_null_indices(temperature,fill_value=0.0)
            nrd.fix_cube_coordinates(temperature)
            temperature = cut_to_limits(temperature)

            salinity = nrd.remove_null_indices(salinity,fill_value=0.0)
            nrd.fix_cube_coordinates(salinity)
            salinity = cut_to_limits(salinity)
            salinity = sou.abs_sal_from_pract_sal(salinity) 
            # Now we should have salinity and temperature set.

            v = sou.cube_volumes(temperature)
            depth = v.coord('depth').points
            thicknesses = sou.cube_cell_thicknesses(v,return_dictionary=True)
            total_volume = np.sum(v.data)*1e-9  # 1e-9 as w want km^3

            pressure = sou.cube_pressure(salinity).data
            d = sou.cube_density(salinity, temperature)
            total_mass = np.sum(v.data*d.data)*1e-3  
            salt_content = salinity.copy()
            salt_content.data = (salinity.data*v.data*d.data)*1e-6
            total_salt = np.sum(salt_content.data)
            average_salt = total_salt/total_mass
            print("Total salt {:.1f} GT\tAverage salinity {:.2} g/kg"\
                    .format(total_salt*1e-9,average_salt*1e3))
            d = salinity.collapsed(collapse_style['coords'],
                                        iris.analysis.MEAN,
                                        weights = v.data)
        all_coords=[]
        for coord in d.coords():
            if coord.name() is not 'time' and coord.shape[0]>1:
                all_coords.append(coord.name())
        iris_list.append(d)
        print("Analysing {} ({} of {})".format(f, num+1, len(files_working)))
    iris_salt_content = siri_omen.concatenate_cubes(iris_list)
    iris_salt_content.attributes['name'] = "{}-salinity-{}-{}".format(\
                                                        name_marker,\
                                                        startdate, enddate)
    out_file_name = datadir+'{}_{}_{}_{}.nc'.\
                format(analysis_type,variable_salinity,\
                 name_marker, collapse_style['name'])
    iris.save(iris_salt_content, out_file_name)
    print("{}:Cube saved succesfully", out_file_name)
    print(iris_salt_content)
