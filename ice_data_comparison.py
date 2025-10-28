# -*- coding: utf-8 -*-
"""
Created on Mon Mar 22 18:57:56 2021

@author: siirias
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

fig_dpi = 300

for data_set in ['A001', 'B001', 'D001']:
    output_dir = "D:\\Data\\Figures\\SmartSea\\new_run\\ice\\"
    new_dat=pd.read_csv("E:\\SmartSea\\derived_data\\ice_extent_{}.csv".format(data_set),parse_dates=['time'])
    old_dat=pd.read_csv("E:\\SmartSea\\derived_data_old\\ice_extent_{}.csv".format(data_set),parse_dates=['time'])
    
    l = len(new_dat.time)
    
    fig,ax1 = plt.subplots(figsize=[20,10])
    plt.plot(old_dat.time[:l], old_dat.ice_extent[:l],'r',label='old')
    plt.plot(new_dat.time[:l], new_dat.ice_extent[:l],'b',label='new')
    plt.title(data_set)
    plt.legend(loc='upper left')
    ax2 = ax1.twinx()
    #plt.figure()
    plt.plot(old_dat.time[:l], old_dat.ice_extent[:l]-new_dat.ice_extent[:l],
             'k',label='old - new', alpha=0.4)
    plt.legend(loc='upper right')
    print("{} difference in averages: {}, ({}%)".format(
    data_set,
    np.mean(old_dat.ice_extent[:l]-new_dat.ice_extent[:l]),
    100.0*(1.0-np.mean(old_dat.ice_extent[:l])/np.mean(new_dat.ice_extent[:l]))
    ))
    
#    print("{} difference in averages: {}, ({}%)".format(
#    data_set,
#    np.mean(np.abs(old_dat.ice_extent[:l]-new_dat.ice_extent[:l])),
#    100.0*(np.mean(np.abs(old_dat.ice_extent[:l]-new_dat.ice_extent[:l]))/np.mean(old_dat.ice_extent[:l]))
#    ))
    
    filename = "{}_ice_concentration.png".format(data_set)
        
    plt.savefig(output_dir+filename,\
                    facecolor='w',dpi=fig_dpi,bbox_inches='tight')
    print("Saved: {}".format(\
          output_dir+filename))
