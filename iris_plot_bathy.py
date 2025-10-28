"""
A simple script to write out nemo bathymetry as an xml text-file.
Created xml file can be opened with Bathymetry Editor BatMan.
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
file_name = 'bathy_meter.nc'

data = iris.load(data_folder+file_name,'Bathymetry')[0]
bathy_data = data.data
#plt.figure(figsize=(5,10))
#iqplt.contourf(data,np.arange(0,300,5))
iqplt.contour(data,[0.],colors='k',linewidth=1.)

plt.show()
