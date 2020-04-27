import os
import re
import subprocess as sp
import time
import process_helper as ph

processes = ph.ProcessCounter(16)
stats = ['max','min','mean']
in_dir = '/scratch/project_2001635/siiriasi/smartsea_data/daily_test/'

for stat in stats:
    cdo_oper = {\
        'max':'timmax',
        'min':'timmin',
        'mean':'timmean'}

    out_dir = '/scratch/project_2001635/siiriasi/smartsea_data/daily_test/{}/'.format(stat)
    in_files = ph.get_file_list(in_dir,".*nc")
    for f in in_files:
        f_name_no_ext = re.search("^([^\.]*)",f).groups()[0]
        f_extension = re.search("([^\.]*)$",f).groups()[0]
        f_out_name = "{}_monthly_{}.{}".format(\
                    f_name_no_ext,\
                    stat,\
                    f_extension)
        the_command = "cdo {} {}/{} {}/{}".format(\
                    cdo_oper[stat], \
                    in_dir, f, \
                    out_dir, f_out_name) 
        processes.add(the_command)
print("All done!")
