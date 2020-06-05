import os
import re
import subprocess as sp
import time
import process_helper as ph

one_time_processes = 6 
processes = ph.ProcessCounter(one_time_processes)

all_series = ['A001', 'A002']
all_variables = ['icevolume'] # kd = vertical light attenuation
#all_depths = ['0-3', '0-12', '12-80', '80-inf']
all_depths = ['']
file_filter = '.*icevolume.*nc'  #regular expression
year_search = '(\d\d\d\d)'
depth_lims = {\
'':'',\
'0-3':'deptht,0,0',\
'0-12':'deptht,0,5',\
'0-9':'deptht,0,2',\
'12-80':'deptht,6,19',\
'80-inf':'deptht,19,35'
}
for series in all_series:
    for variable in all_variables:
        in_dir = '/scratch/project_2001635/siiriasi/smartsea_data/daily_test/'
        out_dir = '/scratch/project_2001635/siiriasi/smartsea_data/daily_test/mean/'
        for depth in all_depths:
            depth_lim = depth_lims[depth]
            if(depth_lim == ''):
                depth_switch = ''
            else:
                depth_switch = "-d {}".format(depth_lim)
                depth = depth + '_'
            files = ph.get_file_list(in_dir,file_filter)
            for in_file in files:
                year=re.search(year_search,in_file).groups()[0]
                for m in range(11):
                    month = m + 1
                    out_file = '{}/{}_{}{}_{}_{}.nc'.format(
                            out_dir,\
                            variable, depth, series,\
                             year, month)
                    the_command = 'cdo timselmean,7 {} {}'.format(\
                                    '{}/{}'.format(in_dir,in_file),\
                                    out_file                                    
                                    ) 
#                    the_command = 'ncea -O  -d time_counter,{},{} -v {} {} {}'.format(\
 #                                    depth_switch, m, m, \
 #                                    variable,\
 #                                    '{}/{}'.format(in_dir,in_file),\
 #                                     out_file)
                    print(the_command)
                    processes.add(the_command)
# in_dir='/scratch/project_2001635/siiriasi/smartsea_data/daily_test/'
# out_dir='/scratch/project_2001635/siiriasi/smartsea_data/daily_test/means/'
# in_files=$(ls $in_dir/*.nc  --color=none)
# for f in $in_files; do
#     f_name_full=$(basename $f)
#     f_extension="${f_name_full##*.}"
#     f_name="${f_name_full%.*}"
#     f_out_name=${f_name}_weekly_mean.${f_extension}
#     cdo timselmean,7 $f $out_dir/$f_out_name 
# done

