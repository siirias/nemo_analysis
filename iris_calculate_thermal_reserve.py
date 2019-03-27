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
import siri_omen.nemo_reader as nrd
import cf_units
import iris.util

ss = smh()
ss.grid_type = 'T'
ss.interval = 'm'
ss.file_name_format = 'NORDIC-GoB_1{}_{}_{}_grid_{}.nc'
ss.root_data_in = "/lustre/tmp/siirias/o/tmp/"  # gludge as the main disk is not sure enough.
# folder_start = 'OUTPUT'
name_markers = ['new_REANALYSIS']
# name_markers = ['REANALYSIS', 'A001', 'D001', 'C001']
# name_markers = ['A001', 'B001', 'C001']
# variables = ['SSS', 'SST', 'SSH_inst']
variable_temperature = 'potential_temperature'
variable_salinity = 'salinity'
# variables= ['vosaline', 'potential_temperature']
make_climatology = True
# name_marker = 'A001'  # this one tells which dataseries is handled.
climatology_time_slots = 366
if(ss.interval =='m'):
    climatology_time_slots = 12
    
for name_marker in name_markers:
    folder_start = 'OUTPUT'
    ss.save_interval = 'month'
    if '1' in name_marker: # the 001 series are hindcasts, all other scenarios
        startdate = datetime.datetime(1975, 1, 1)
        enddate = datetime.datetime(2005, 12, 31)
    elif 'REANALYSIS' in name_marker:
        startdate = datetime.datetime(1980, 1, 1)
        enddate = datetime.datetime(2012, 12, 31)
        ss.save_interval = 'year'
        folder_start = ''
        if 'new_' in name_marker:
          ss.file_name_format = 'SS-GOB_1{}_{}_{}_grid_{}.nc'
    else:
        startdate = datetime.datetime(2006, 1, 1)
        enddate = datetime.datetime(2058, 12, 31)
    datadir = ss.root_data_out+"/derived_data/" # where everyt output is stored
    
    ss.main_data_folder= ss.root_data_in+"/{}{}/".format(folder_start, name_marker)
    depth_ax = 'deptht'
    if variable_temperature in ['SST', 'SSS', 'SSH_inst']:
        ss.grid_type = 'T'
        ss.interval = 'd'
    if variable_temperature in ['vosaline', 'potential_temperature']:
        ss.grid_type = 'T'
        ss.interval = 'm'
        depth_ax = 'deptht'
    
    
    if '1' in name_marker: # the 001 series are hindcasts, all other scenarios
        series_name = 'Hindcast_{}_{}'.format(name_marker, variable_temperature)
    elif 'REANALYSIS' in name_marker:
        series_name = '{}_{}'.format(name_marker, variable_temperature)
    else:
        series_name = 'Scenario_{}_{}'.format(name_marker, variable_temperature)
    
    enddate = datetime.datetime(1982, 12, 31) # REMOVE!
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
        p_temperature = iris.load(ss.main_data_folder+f, variable_temperature)[0]
        p_temperature = nrd.remove_null_indices(p_temperature)
        nrd.fix_cube_coordinates(p_temperature)
#        p_temperature.convert_units('K')  # Kelvins, to be sure
        
        salinity = iris.load(ss.main_data_folder+f, variable_salinity)[0]
        if salinity.units.name == 'unknown':
            salinity.units='1e-3'
        salinity = nrd.remove_null_indices(salinity)
        nrd.fix_cube_coordinates(salinity)
        # Now we should have salinity and tempereture set.

        volumes = siri_omen.utility.cube_volumes(p_temperature)
        # now we have the volume set.

        # TODO: figure out the actual formula for heat content
        # Should be Volume*density*specific heat*potential_temperature
        # actually: TEOS-10 prefers to use Conservative Temperature
        # instead of potential temperature. gws.CT_from_pt, and after that:
        # "Conservative Temperature accurately represents the Heat Content
        # per unit of mass of seawater"
        # as such equation should be: Volume*Density*CT
        p_temperature.data=gsw.CT_from_pt(salinity.data,p_temperature.data)
        #convert to Conservative Temperature
        d = p_temperature.collapsed(['longitude', 'latitude'], iris.analysis.SUM, weights = volumes)
        iris_list.append(d)
        print("Analysing {} ({} of {})".format(f, num+1, len(files_working)))
    iris_p_temperature = siri_omen.concatenate_cubes(iris_list)
    iris_p_temperature.attributes['name'] = "{}-{}-{}-{}".format(\
                                                        name_marker,\
                                                        variable_temperature,\
                                                        startdate, enddate)
    out_file_name = datadir+'reserve_{}_{}.nc'.format(variable_temperature, name_marker)
    iris.save(iris_p_temperature, out_file_name)
    print("{}:Cube saved succesfully", out_file_name)
    print(iris_p_temperature)
