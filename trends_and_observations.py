# -*- coding: utf-8 -*-
"""
Created on Wed Mar 22 11:00:46 2023

@author: siirias
"""
import datetime as dt
import matplotlib as mp
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
#import seaborn as sns
import os
import re
import netCDF4
import smartseahelper
from netCDF4 import Dataset
import xarray as xr
from pandas.plotting import register_matplotlib_converters
register_matplotlib_converters()

import warnings
warnings.filterwarnings("ignore")
def read_observations(measuring_point, variable, depth, depth_window = 3.0):
    sm = smartseahelper.smh()
    sm.root_data_in = "C:/Data/SmartSeaModeling/"
    sm.root_data_out = "C:/Data/SmartSeaModeling/"
    out_dir = sm.root_data_out+"/tmp/"
    in_dir = sm.root_data_in + "/observations/"
    model_in_dir = sm.root_data_in + "/Extracted_profiles/"
    fig_factor = 1.0#1.5 #0.8  #1.5
    fig_size = (10*fig_factor,5*fig_factor)
    
    obs_dtypes = {
        'Cruise': str,
        'Station': str,
        'Type': str,
        'yyyy-mm-ddThh:mm': str,
        'Latitude [degrees_north]': float,
        'Longitude [degrees_east]': float,
        'Bot. Depth [m]': float,
        'Secchi Depth [m]:METAVAR:FLOAT': float,
        'PRES [db]': float,
        'TEMP [deg C]': float,
        'PSAL [psu]': float,
        'DOXY [ml/l]': float,
        'PHOS [umol/l]': float,
        'TPHS [umol/l]': float,
        'SLCA [umol/l]': float,
        'NTRA [umol/l]': float,
        'NTRI [umol/l]': float,
        'AMON [umol/l]': float,
        'NTOT [umol/l]': float,
        'PHPH []': float,
        'ALKY [meq/l]': float,
        'CPHL [ug/l]': float
    }
    obs_nans = ['Nan']
    
    def handle_specials(d):
        #get rid of non numeric strings in dataobs
        if(d==''):
            return np.nan
        if(d[0]=='<'):
            return 0.0
        return float(d)
    
    convert_dict = {col: handle_specials 
                    for col in obs_dtypes.keys() if obs_dtypes[col] == float}
    
    depth_variable = 'PRES [db]'
    time_variable = 'yyyy-mm-ddThh:mm'
    T_variable = 'TEMP [deg C]'
    S_variable = 'PSAL [psu]'
    
    all_variables = [{'obs':T_variable, 'model':'votemper'},
                     {'obs':S_variable, 'model':'vosaline'}]
    var = list(filter(lambda x: x['model']==variable, all_variables))[0]    
    file = measuring_point+'.csv'
    data = pd.read_csv(in_dir+file, 
                       parse_dates = [time_variable],
                       na_values = obs_nans,
                       converters = convert_dict,
                       dtype = obs_dtypes)
    indices = np.abs(data[depth_variable]-d) <= depth_window
    d_to_plot = data.loc[indices, var['obs']]
    t_to_plot = data.loc[indices, time_variable]
    return t_to_plot, d_to_plot


point_depths = [[0.5,245.0],[0.5,110.0],[0.5,120.0],[0.5,100.0]]
for measuring_point,depths in zip(['F64','SR5','US5B','BO3'],point_depths):
    for var in ['vosaline','votemper']:
        for d in depths:
            t_to_plot,d_to_plot = read_observations(measuring_point, var, d)
            plt.figure()
            plt.plot(t_to_plot, d_to_plot,'k.')
            plt.title(f'{var},{measuring_point}, {d} db')
            plt.xlabel('Time')
            plt.ylabel(var)
