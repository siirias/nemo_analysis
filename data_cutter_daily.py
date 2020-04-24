import os
import re
import subprocess as sp
import time
#all_series=['A001', 'A002', 'A005', 'B001', 'B002', 'B005', 'D001', 'D002', 'D005']
#all_series=['B001']
all_series=['A001', 'A002', 'A005', 'B002', 'B005', 'D001', 'D002', 'D005']
all_series=['D005']
#all_variables=['SST', 'SSS']
all_variables=['SBS', 'icethic', 'icecon']
monthly = True
one_time_processes =  16
class ProcessCounter:

    def __init__(self, max_processes):
        self.max_processes = max_processes
        self.process_list = []
        self.time_started = 0.0
        self.time_elapsed = 0.0

    def add(self, new_process):
        # adds process, and if there is too many, waits for 
        # them to finnish, then clears the list.
        if len(self.process_list) == 0:
            self.time_started = time.time()
        self.process_list.append(new_process)
        if(len(self.process_list)>=self.max_processes):
            print("Waiting...",len(self.process_list))
            for p in self.process_list:
                p.wait()
            self.process_list = []
            self.time_elapsed = time.time() - self.time_started
            self.time_started = time.time()
            print("Waited {:0.3f} s, wait per process {:0.3f} s".format(\
                    self.time_elapsed, \
                    self.time_elapsed/self.max_processes))

processes = ProcessCounter(one_time_processes)

for series in all_series:
    for variable in all_variables:
        in_dir='/scratch/project_2001635/siiriasi/smartsea_data/'+series+'/'
        out_dir='/scratch/project_2001635/siiriasi/smartsea_data/daily_test/'
        file_type = ".*_1d_.*_grid_T.nc"
        files = os.listdir(in_dir)
        files = [f for f in files if re.match(file_type,f)]
        for in_file in files:
            print(in_file)
            year=re.search('_1[dmh]_(....)',in_file).groups()[0]
            for m in range(12):
                month = m + 1
                out_file="{}_{}_{}_{}.nc".format(\
                        variable, series, year, month)
                the_command = "ncea -O -v {} {} {}".format(\
                            variable,\
                            in_dir+in_file,\
                            out_dir+out_file)
                sh_output = sp.Popen(the_command,\
                                universal_newlines = True,\
                                stdout =sp.PIPE,\
                                shell = True)
                processes.add(sh_output)
