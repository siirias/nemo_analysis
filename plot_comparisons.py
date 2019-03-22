#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 13 16:38:59 2018

@author: siirias
"""
import datetime as dt
import calendar
import matplotlib as mp
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from smartseahelper import smh
ss=smh()
colors=['r','g','b','m']

datadir=ss.root_data_out+"/derived_data/"
sets=['A001','B001','C001']
fig1=plt.figure(figsize=(6,6))
fig2=plt.figure(figsize=(16,6))

for i,the_set in zip(range(len(sets)),sets):
    yearly_data=pd.read_csv(datadir+'yearly_max_ice_{}.csv'.format(the_set))
    full_data=pd.read_csv(datadir+'ice_extent_{}.csv'.format(the_set),parse_dates=['time'])
    
    plt.figure(fig1.number)
    plt.plot(yearly_data['year'],yearly_data['max_ice'],colors[i],label=the_set)
    
    plt.figure(fig2.number)
    plt.plot(full_data['time'],full_data['ice_extent'],colors[i],label=the_set)

plt.legend()
plt.figure(fig1.number)
plt.legend()
