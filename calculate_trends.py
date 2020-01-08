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
import netCDF4 as nc4
from netCDF4 import Dataset
from smartseahelper import smh
import os
import cmocean
import pandas as pd
import re
plot_stuff=True
in_dir = "/arch/smartsea/analysis/derived_data/test/"
out_dir = "/arch/smartsea/analysis/derived_data/test/output/"

out_file = open(out_dir+"statistics.txt",'w')
#fname_vars=["SSS", "SST", "SBS", "SBT", "vosaline", "votemper" ,"icevolume"]
fname_vars=["SSS", "SST", "SBS", "SBT", "vosaline", "votemper"]
#fname_var = ["SBT"]
for fname_var in fname_vars:

    if(plot_stuff):
        plt.close('all')
        plt.figure(1)
        plt.title("{}, Historical runs".format(fname_var))
        plt.figure(2)
        plt.title("{}, RCP 4.5".format(fname_var))
        plt.figure(3)
        plt.title("{}, RCP 8.5".format(fname_var))
        target_fig={'1':1,'2':2,'5':3}
        target_col={'A':'k','B':'b','C':'c','D':'g'}

    in_files =[ 
            "yearly_{}_B001.csv".format(fname_var),
            "yearly_{}_A001.csv".format(fname_var),
            "yearly_{}_C001.csv".format(fname_var),
            "yearly_{}_D001.csv".format(fname_var),
            "yearly_{}_B002.csv".format(fname_var),
            "yearly_{}_A002.csv".format(fname_var),
            "yearly_{}_C002.csv".format(fname_var),
            "yearly_{}_D002.csv".format(fname_var),
            "yearly_{}_D005.csv".format(fname_var),
            "yearly_{}_A005.csv".format(fname_var),
            "yearly_{}_B005.csv".format(fname_var)]
    results={}
    mean_results={}
    for in_file in in_files:
        if(os.path.isfile(in_dir+in_file)):
            data_set=re.search("[ABCD]\d\d\d",in_file).group()
            clim_mod=data_set[0]
            clim_scen=data_set[-1]
            data=pd.read_csv(in_dir+in_file)
            fitted=np.polyfit(data['year'],data['mean'],1)
            change_dec=fitted[0]*10
            if(plot_stuff):
                y_limits = [None,None]
                if( fname_var in ['SST', 'SBT', 'votemper']):
                    y_limits = [4.0, 10.0]
                if( fname_var in ['SBS', 'vosaline']):
                    y_limits = [3.0, 8.0]
                if( fname_var in ['SSS']):
                    y_limits = [3.0, 6.0]
                plt.figure(target_fig[clim_scen])
                ydat = data['year']
                plt.plot(ydat,data['mean'],target_col[clim_mod],label = "",alpha=0.1)
                plt.plot(ydat,data['mean'],target_col[clim_mod],marker='.',\
                                label = "{}({:.2}/dec)".format(data_set,change_dec),alpha=0.1)
                plt.plot(ydat,ydat*fitted[0]+fitted[1],target_col[clim_mod], label = "")
                plt.xlabel('year')
                plt.ylabel(fname_var)
                plt.gca().set_ylim(y_limits[0],y_limits[1])
            print("{}: Change of '{}' is {:.3} per decade".format(\
                                        data_set,\
                                        fname_var,\
                                        change_dec))
            results[data_set]=change_dec
            if(clim_scen not in mean_results.keys()):
                mean_results[clim_scen]={}
                mean_results[clim_scen]['sum']=change_dec
                mean_results[clim_scen]['num']=1
            else:
                mean_results[clim_scen]['sum']+=change_dec
                mean_results[clim_scen]['num']+=1
        else:
            print("{} not found.".format(in_dir+in_file))
    scen_names={'1':'History','2':'RCP4.5','5':'RCP8.5'}
    for i in mean_results:
        mean_results[i]['mean'] = mean_results[i]['sum']/mean_results[i]['num']
    for i in mean_results:
        print("{}: {} changes {:.3} per decade on average".format(\
                                    scen_names[i],\
                                    fname_var,\
                                    mean_results[i]['mean']))

    for i in results:
        out_file.write("{}\t{}\t{}\n".format(i,fname_var,results[i]))

    for i in mean_results:
        out_file.write("{}\t{}\t{}\n".format(scen_names[i],fname_var,mean_results[i]['mean']))

    if(plot_stuff):
        plt.figure(1)
        plt.legend()
        plt.savefig(out_dir+'history_'+fname_var+'.png')
        plt.figure(2)
        plt.legend()
        plt.savefig(out_dir+'RCP45_'+fname_var+'.png')
        plt.figure(3)
        plt.legend()
        plt.savefig(out_dir+'RCP85_'+fname_var+'.png')
out_file.close()
