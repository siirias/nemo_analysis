import os
import re
import subprocess as sp
import time


def get_file_list(in_dir, file_regex ='.*'):
    # get_file_list(in_dir, file_regex)
    # Shortcut function to get file names matching
    # given regexp in given dir
    files = os.listdir(in_dir)
    files = [f for f in files if re.match(file_regex,f)]
    return files


class ProcessCounter:

    def __init__(self, max_processes):
        self.max_processes = max_processes
        self.process_list = []
        self.time_started = 0.0
        self.time_elapsed = 0.0

    def add(self, new_command):
        # adds process, and if there is too many, waits for 
        # them to finnish, then clears the list.
        new_process = sp.Popen(new_command,\
                        universal_newlines = True,\
                        stdout =sp.PIPE,\
                        shell = True)
        if len(self.process_list) == 0:
            self.time_started = time.time()
        self.process_list.append(new_process)
        if(len(self.process_list)>=self.max_processes):
            print("Waiting...",len(self.process_list))
            self.wait_all()

    def wait_all(self):
        for p in self.process_list:
            p.wait()
            self.process_list = []
            self.time_elapsed = time.time() - self.time_started
            self.time_started = time.time()
            print("Waited {:0.3f} s, wait per process {:0.3f} s".format(\
                    self.time_elapsed, \
                    self.time_elapsed/self.max_processes))
