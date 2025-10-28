# !/usr/bin/env ipython
# -*- coding: utf-8 -*-
"""
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


data_dir = "/arch/smartsea/analysis/derived_data/"
files = [
         "reserve_potential_temperature_new_REANALYSIS_depthlat.nc",
         "reserve_potential_temperature_REANALYSIS_SMHI_depthlat.nc",
         "reserve_potential_temperature_REANALYSIS_depthlat.nc"
        ]
labels = ["FMI2","SMHI","FMI1"]
styles = ['b','r','b--']
data=[]
cut_depth = 40.0

#set_name =""
#latitude_min = 50.0
#latitude_max = 75.0

set_name ="60to60_2"
latitude_min = 60.0
latitude_max = 60.2

#set_name ="60to61"
#latitude_min = 60.0
#latitude_max = 61.0

#set_name ="63to66"
#latitude_min = 63.0
#latitude_max = 66.0


heat_total=[]
heat_top=[]
heat_bottom=[]
datemin = iris.time.PartialDateTime(year=1980, month=1, day=1, hour=0)
datemax = iris.time.PartialDateTime(year=2008, month=12, day=31, hour=0)
def plot_changes(cubes, offset=0.0, title=""):
    plt.figure(figsize=(15,5))
    for cube,label,style in zip(cubes,labels,styles):
        iplt.plot(cube-offset,style,label=label,linewidth=1)
    plt.title(title)
    plt.legend()
    

def plot_difference(cube1, cube2, title=""):
    plt.figure(figsize=(15,5))
    iplt.plot(cube1*0.0,'g-',linewidth=1,alpha=0.3)
    iplt.plot(cube1-cube2,styles[0],linewidth=1)
    plt.title(title)


def compare_two(first,second):
    #gludge code to avoid too much rewriting
    plot_difference(heat_total[second],heat_total[first],\
                    title="Difference Total {}-{} {}"\
                    .format(labels[second],labels[first],set_name))
    plt.savefig(data_dir+"heat_content_difference_{}_{}_total{}.png".\
                format(labels[first],labels[second],set_name))
    plt.draw()

    plot_difference(heat_top[second],heat_top[first],\
                    title="Difference Top 0-{} m {}-{} {}".\
                    format(cut_depth,labels[second],labels[first],set_name))
    plt.savefig(data_dir+"heat_content_difference_{}_{}_top{}.png".\
                format(labels[first],labels[second],set_name))
    plt.draw()

    plot_difference(heat_bottom[second],heat_bottom[first],\
                    title="Difference Bottom {}+ m {}-{} {}"\
                    .format(cut_depth,labels[second],labels[first],set_name))
    plt.savefig(data_dir+"heat_content_difference_{}_{}_bottom{}.png".\
                format(labels[first],labels[second],set_name))
    plt.draw()
    
for ifile in files:
    data.append(iris.load(data_dir+ifile)[0])
    data[-1] = data[-1].extract(iris.Constraint(time = \
                 lambda t: datemin < t < datemax))
    data[-1] = data[-1].extract(iris.Constraint(latitude =\
                 lambda l: latitude_min < l < latitude_max))
    data[-1] = data[-1].collapsed('latitude',iris.analysis.SUM)
    heat_total.append(data[-1].collapsed('depth',iris.analysis.SUM))
    heat_top.append(data[-1].extract(
                iris.Constraint(coord_values = 
                {'depth': lambda d: d<cut_depth})).collapsed('depth',iris.analysis.SUM))
    heat_bottom.append(data[-1].extract(
                iris.Constraint(coord_values = 
                {'depth': lambda d: d>=cut_depth})).collapsed('depth',iris.analysis.SUM))
plot_changes(heat_total,title="Total {}".format(set_name))
plt.savefig(data_dir+"heat_content_total{}.png".format(set_name))
plt.draw()

plot_changes(heat_top,title="Top 0-{} m {}".format(cut_depth, set_name))
plt.savefig(data_dir+"heat_content_top{}.png".format(set_name))
plt.draw()

plot_changes(heat_bottom,title="Bottom {}+ m {}".format(cut_depth,set_name))
plt.savefig(data_dir+"heat_content_bottom{}.png".format(set_name))
plt.draw()

compare_two(0,1)
compare_two(2,0)



plt.show()
