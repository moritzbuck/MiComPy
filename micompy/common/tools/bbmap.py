import subprocess
import os
import shutil
from subprocess import call, check_output, Popen, PIPE
import json
from os.path import join as pjoin

class bbmap(object):
    def __init__(self, executable = "bbmap.sh", jni=True ):
        self.executable = executable
        self.jni = "jni=t" if jni else "jni=f"
        try:
            with open(os.devnull, 'w') as handle:
                subprocess.call([self.executable, "-h"], stdout=handle)
        except OSError as e:
            print("bbmap not found")


    def make_index(self, genome):
            with open(os.devnull, 'w') as handle:
                call([self.executable, "ref="+genome.genome, "path="+genome.path], stdout=handle, stderr=handle)

    def get_ANI(self, genome, reads, target_id = 0.1):
        FNULL = open(os.devnull, 'w')
        cmd = Popen([self.executable, "ref="+genome.genome, "path="+genome.path, "in="+reads, "out=/dev/null", "minid=" +str(target_id), self.jni], stderr = PIPE, stdout=PIPE)
        out_lines = cmd.stderr.readlines()
        FNULL.close()

        coverage = [float(l.split("\t")[1].replace("%","")) for l in out_lines if "mapped:" in l][0]
        ANI = [float(l.split("\t")[3].replace("%","")) for l in out_lines if "Match Rate:" in l][0]

        return {'ANI' : ANI, 'coverage' : coverage}
