import datetime
import numpy as np
from scipy.io import netcdf
from netCDF4 import Dataset
from smartseahelper import smh
import os
import cmocean
import re
from datetime import datetime as dt
from datetime import timedelta

start_time = dt.strptime('2006-01-01','%Y-%m-%d')
end_time = dt.strptime('2060-01-01','%Y-%m-%d')
out_file_name = "scenario_timestamps_daily.txt"
now_time = start_time
with open(out_file_name,'w') as out_file:
    while(now_time<end_time):
        print(now_time)
        out_file.write(str(now_time)+'\n')
        now_time += timedelta(days=1)

