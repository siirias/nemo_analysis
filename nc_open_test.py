# -*- coding: utf-8 -*-
"""
Created on Mon Oct  5 11:11:13 2020

@author: siirias
"""
import xarray as xr
import iris  # if this causes problems, comment this, and the iris test.

in_file = "icevolume_A005_2040.nc"
in_dir = "C:\\Users\\siirias\\Downloads\\anton_testi\\"


#test with xarray, should only need the xarray package
print("******Loading with XArray******")
d=xr.open_dataset(in_dir+in_file)
print(d.data_vars) #show the available variables
print(d.icevolume.data[0,180,170]) #example of extracting one value, t,x,y
d.close()

# test with iris. Requires iris, bit tricier sometimes to get to work
# installing iris within anaconda/miniconda should work, sometimes problems
# with udunits2 there might come error: 
#  [UT_OPEN_ARG] Failed to open UDUNITS-2 XML unit database: "b'No such file or directory'"
# in this case, installation has gone somehow wrong, and the directory
# is not where conda/python thinks it is.
# easiest is to find where udunits2.xml is, and just copy that and other
# xmls foudn in same locations to where system thinks tehy should be.
# for example, in my case they were:
# C:\Users\siirias\AppData\Local\Continuum\anaconda3\pkgs\udunits2-2.2.27.6-h252784a_1001\Library\share\udunits
# and system search them from:
# C:\Users\siirias\AppData\Local\Continuum\anaconda3\share\udunits

# ACTUAL TEST STARTS
variables = iris.load(in_dir+in_file)
print("******Loading with Iris*******")
print(variables) # lists variables found (iris calls them cubes) 
                 # should have just one in the original  example
print(variables[0].name())
print(variables[0].data[0,180,170])