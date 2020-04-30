import os
import re
import subprocess as sp
import time
import process_helper as ph
import glob
#NOTE: For some reason the zipping halts with this funciton.
# I used it just to print the commands so I could copy paste
# those to prompt....
skip_processed = True
skip_original = False
if( not skip_processed):
    processes = ph.ProcessCounter(3)
    stats = ['max','min','mean']
    variables = ['SSS', 'SST', 'SBS', 'icethic', 'icecon', 'snowthic', 'icevvelo', 'iceuvelo', 'sorunoff']
    in_dir = '/scratch/project_2001635/siiriasi/smartsea_data/daily_test/'
    for var in variables:
        for stat in stats:
    #        print("Handling {}, {}".format(var,stat))
            out_dir = '/scratch/project_2001635/siiriasi/smartsea_data/daily_test/{}/'.format(stat)
    #        the_command = "zip -v {}/{}_monthly_{}.zip {}".format(\
    #                    out_dir, var, \
    #                    stat, \
    #                    " ".join(glob.glob("{}/{}*.nc".format(out_dir,var)))) 
            the_command = "zip -v {}/{}_monthly_{}.zip {}".format(\
                        out_dir, var, \
                        stat, \
                        " {}/{}*.nc".format(out_dir,var)) 
            print(the_command)
    #        processes.add(the_command)
    processes.wait_all()

if(not skip_original):
    stats = ['max','min','mean']
    variables = ['SSS', 'SST', 'SBS', 'icethic', 'icecon', 'snowthic', 'icevvelo', 'iceuvelo', 'sorunoff']
    in_dir = '/scratch/project_2001635/siiriasi/smartsea_data/daily_test/'
    for var in variables:
        for stat in stats:
            out_dir = '/scratch/project_2001635/siiriasi/smartsea_data/daily_test/{}/'.format(stat)
            the_command = "zip -v {}/{}_monthly_{}.zip {}".format(\
                        out_dir, var, \
                        stat, \
                        " {}/{}*.nc".format(out_dir,var)) 
            print(the_command)
    #        processes.add(the_command)
    processes.wait_all()
    

print("All done!")
