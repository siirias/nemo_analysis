import os
import re
import subprocess as sp
import time
import process_helper as ph
#all_series=['A001', 'A002', 'A005', 'B001', 'B002', 'B005', 'D001', 'D002', 'D005']
all_series=['A001', 'A002', 'A005', 'B001', 'B002', 'B005','D001', 'D002', 'D005']
#all_series = ['A001']
in_dir_root = '/scratch/project_2001635/siiriasi/smartsea_data/'
#in_dir_root = '/scratch/project_2001635/siiriasi/smartsea_data/forcings/'
#all_variables=['SST', 'SSS', 'SBS', 'icethic', 'icecon', 'iceuvelo', 'icevvelo', 'snowthic']
#all_variables=['sorunoff']
all_variables=['soicecov', 'icevolume']
file_type = '.*_1d_.*_grid_T.nc'   # daily data
year_search = '_1[dmh]_(....)'     # daily data
time_axis = 'time_counter'
#file_type = '.*\.nc'   # river forcing
#year_search = '_y(....)'     # river forcing
#time_axis = 'time'
monthly = False 
one_time_processes = 6 

processes = ph.ProcessCounter(one_time_processes)

for series in all_series:
    for variable in all_variables:
        in_dir = in_dir_root+series+'/'
        out_dir='/scratch/project_2001635/siiriasi/smartsea_data/daily_test/'
        files = os.listdir(in_dir)
        files = [f for f in files if re.match(file_type,f)]
        for in_file in files:
            print(in_file)
            year=re.search(year_search,in_file).groups()[0]
            if monthly:
                    for m in range(12):
                        month = m + 1
                        first_step = int((365.25/12.0)*m)
                        last_step = min(364,int((365.25/12.0)*(m+1)))
                        out_file="{}_{}_{}_{}.nc".format(\
                                variable, series, year, month)
                        the_command = "ncea -O -v {} -d {},{},{} {} {}".format(\
                                    variable,\
                                    time_axis, first_step, last_step,\
                                    in_dir+in_file,\
                                    out_dir+out_file)
                        print(the_command)
                        processes.add(the_command)
            else: #don't split to months
                    out_file="{}_{}_{}.nc".format(\
                            variable, series, year)
                    the_command = "ncea -O -v {} {} {}".format(\
                                variable,\
                                in_dir+in_file,\
                                out_dir+out_file)
                    print(the_command)
                    processes.add(the_command)
 
