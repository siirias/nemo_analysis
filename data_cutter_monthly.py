import os
import re
import subprocess as sp
import time
import process_helper as ph

one_time_processes = 6 
processes = ph.ProcessCounter(one_time_processes)

all_series = ['D007', 'D002', 'D004', 'D005', 'D001']
all_variables = ['kd'] # kd = vertical light attenuation
all_depths = ['0-3', '0-12', '12-80', '80-inf']
file_filter = '.*_1m_.*_ptrc_T.nc'  #regular expression
year_search = '1m_(\d\d\d\d)'
depth_lims = {\
'0-3':'deptht,0,0',\
'0-12':'deptht,0,5',\
'0-9':'deptht,0,2',\
'12-80':'deptht,6,19',\
'80-inf':'deptht,19,35'
}
for series in all_series:
    for variable in all_variables:
        in_dir = '/scratch/project_2001635/siiriasi/smartsea_data/{}'.format(series)
        out_dir = '/scratch/project_2001635/siiriasi/smartsea_data/syke_test/'
        for depth in all_depths:
            depth_lim = depth_lims[depth]
            files = ph.get_file_list(in_dir,file_filter)
            for in_file in files:
                year=re.search(year_search,in_file).groups()[0]
                for m in range(11):
                    month = m + 1
                    out_file = '{}/{}_{}_{}_{}_{}.nc'.format(
                            out_dir,\
                            variable, depth, series,\
                             year, month)
                    the_command = 'ncea -O -d {} -d time_counter,{},{} -v {} {} {}'.format(\
                                    depth_lim, m, m, \
                                    variable,\
                                    '{}/{}'.format(in_dir,in_file),\
                                     out_file)
                    print(the_command)
                    processes.add(the_command)
