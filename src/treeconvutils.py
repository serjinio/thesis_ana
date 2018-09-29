
import subprocess
import os
import time
import re


class ScattreeProcess:

    command = "src/scattree"
    logs_dir = "scattree_logs"

    def __init__(self, run_no):
        self.run_no = run_no
        self.logfilename = "{}/scattree{}.log".format(
            ScattreeProcess.logs_dir, run_no)
        self._popen = None
        if not os.path.exists(ScattreeProcess.logs_dir):
            os.makedirs(ScattreeProcess.logs_dir)

    def run(self):
        print "Starting conversion for run no. {}".format(self.run_no)
        self.logfile = open(self.logfilename, "w")
        self._popen = subprocess.Popen(
            [ScattreeProcess.command, str(self.run_no)],
            stdout=self.logfile, stderr=self.logfile)

    def check_logfile(self):
        logfile = open(self.logfilename, "r")
        for line in logfile:
            if re.search("error", line, re.IGNORECASE):
              return False
        return True

    def check_result(self):
        if self.logfile:
            self.logfile.close()

        if self._popen and (self._popen.returncode != 0):
            err_msg = ("ERROR: child process for run no. {} "
                       "finished with an error! Return code was {}.".
                       format(self.run_no, self._popen.returncode))
            raise OSError(err_msg)
        elif not self.check_logfile():
            err_msg = ("ERROR: child process for run no. {} "
                       "finished with an error! See its log file: {} "
                       "for details!".format(self.run_no, self.logfilename))
            raise OSError(err_msg)
        else:
            print ("Conversion for run no. {} finished "
                   "successfully.".format(self.run_no))

    def poll(self):
        pollres = self._popen.poll()
        if pollres is not None:
            self.check_result()
        return pollres

    def __str__(self):
        return "scattree for run no. {}".format(self.run_no)


def convert_runs_range(runs_range, max_processes=8):
    processes = set()
    for run_no in runs_range:
        logfile = "scattree{}.log".format(run_no)
        ps = ScattreeProcess(run_no)
        ps.run()
        #ps.communicate()
        processes.add(ps)
        if len(processes) >= max_processes:
            os.wait()
            finished_processes = {p for p in processes if p.poll() is not None}
            processes -= finished_processes
    while(len(processes) > 0):
        os.wait()
        finished_processes = {p for p in processes if p.poll() is not None}
        processes -= finished_processes
