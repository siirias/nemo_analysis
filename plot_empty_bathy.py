# -*- coding: utf-8 -*-
"""
Created on Tue Nov 23 14:26:56 2021

@author: siirias
"""

import matplotlib.pyplot as plt
import xarray as xr
import numpy as np

in_dir = "c:\\Data\\NemoTest\\"
bathy_file = "bathy_meter.nc"

dat = xr.open_dataset(in_dir + bathy_file)
the_mask = np.array(dat.bdy_msk)

the_data = np.random.random(the_mask.shape)
the_data[the_mask == 0.0] = np.nan

plt.imshow(the_data)
plt.gca().invert_yaxis()
plt.gca().set_facecolor('#6a3e25')

