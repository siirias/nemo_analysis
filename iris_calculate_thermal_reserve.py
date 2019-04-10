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
name_markers = ['new_REANALYSIS','REANALYSIS_SMHI','REANALYSIS']
variable_temperature = 'potential_temperature'
variable_salinity = 'salinity'
#collapse_style={'name':'depth','coords':['longitude', 'latitude']}    
#collapse_style={'name':'depthlat','coords':['longitude']}    
#collapse_style={'name':'depthlatlon','coords':[]}    
collapse_style={'name':'total','coords':['longitude','latitude','depth']}    
for name_marker in name_markers:
    folder_start = 'OUTPUT'
    ss.save_interval = 'year'
    if '1' in name_marker: # the 001 series are hindcasts, all other scenarios
        startdate = datetime.datetime(1975, 1, 1)
        enddate = datetime.datetime(2005, 12, 31)
    elif 'REANALYSIS' in name_marker:
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

            salinity = nrd.remove_null_indices(salinity,fill_value=0.0)
            nrd.fix_cube_coordinates(salinity)
            salinity = sou.abs_sal_from_pract_sal(salinity) 
            # Now we should have salinity and temperature set.

            heat_content = sou.cube_heat_content(salinity, temperature)
            d = heat_content.collapsed(collapse_style['coords'],
                                        iris.analysis.SUM)
        all_coords=[]
        for coord in d.coords():
            if coord.name() is not 'time' and coord.shape[0]>1:
                all_coords.append(coord.name())
        if len(all_coords)>0:
            max_heat_content = np.max(d.collapsed(all_coords,iris.analysis.SUM).data)
        else:
            max_heat_content = np.max(d.data)
        print("MAX heat Content:{}".format(max_heat_content))
        iris_list.append(d)
        print("Analysing {} ({} of {})".format(f, num+1, len(files_working)))
    iris_heat_content = siri_omen.concatenate_cubes(iris_list)
    iris_heat_content.attributes['name'] = "{}-heat-content-{}-{}".format(\
                                                        name_marker,\
                                                        startdate, enddate)
    out_file_name = datadir+'reserve_{}_{}_{}.nc'.\
                format(variable_temperature, name_marker, collapse_style['name'])
    iris.save(iris_heat_content, out_file_name)
    print("{}:Cube saved succesfully", out_file_name)
    print(iris_heat_content)
