#!/usr/bin/env ipython
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
#from mpl_toolkits.basemap import Basemap
#from mpl_toolkits.basemap import Basemap, shiftgrid, cm
from netCDF4 import Dataset
from smartseahelper import smh
import os
import cmocean
import re
import pandas as pd


name_marker='REANALYSIS'  #this one tells which dataseries is handled.
compare_two=False
name_marker2='A001'
var1 = 'mean_votemper'  
if 'REANALYSIS' in name_marker:
    data_folder="/{}/".format(name_marker)
else:
    data_folder="/OUTPUT{}/".format(name_marker)
ss=smh()
ss.grid_type='T'
ss.interwall='d'
ss.main_data_folder= ss.root_data_in+data_folder
datadir = ss.root_data_out+"/tmp/Images_for_video2/" #where everyt output is stored
derived_folder=ss.root_data_in+'/derived_data/'

def C2K(x):
    return float(x)+273.15


data1=pd.read_csv(derived_folder+'{}_{}.csv'.format(var1,name_marker))
data2=pd.read_csv(derived_folder+'{}_{}.csv'.format(var1,name_marker2))
data1[var1]=data1[var1].apply(C2K)
data2[var1]=data2[var1].apply(C2K)

fig=data1.plot('time',var1,ylim=(0,290))
data2.plot('time',var1,ax=fig)
plt.show()
print('loppui')
