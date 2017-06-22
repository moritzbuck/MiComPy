import subprocess
import os
import shutil
from subprocess import call, check_output
import json
from os.path import join as pjoin

class MASH(object):


    def __init__(self, executable = "mash", kmer=21, hashes=10000 ):
        self.executable = executable
        self.kmer = kmer
        self.hashes=hashes

        try:
            with open(os.devnull, 'w') as handle:
                subprocess.call([self.executable, "-h"], stdout=handle)
        except OSError as e:
            print("mash not found")


    def run_mash_sketch(self, g, nproc = 1):
        infolder = g.path
        temp_dir = pjoin(g.path,  "checkm_temp")

        with open(os.devnull, 'w') as FNULL:
            out = check_output(["mash", "sketch", "-k", str(self.kmer), "-s", str(self.hashes), g.genome], stderr = FNULL)
