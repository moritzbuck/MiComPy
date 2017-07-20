import subprocess
import os
import shutil
from subprocess import call, Popen, PIPE, check_output
import json
from os.path import join as pjoin
from micompy.common.tools.tool import Tool

class MASH(Tool):

    def __init__(self, executable = "mash", kmer=21, hashes=1000 ):
        Tool.__init__(self, name = "MASH", executable = executable)
        self.kmer = kmer
        self.hashes=hashes

    def run_mash_sketch(self, g, nproc = 1):
        infolder = g.path
        temp_dir = pjoin(g.path,  "checkm_temp")

        with open(os.devnull, 'w') as FNULL:
            out = check_output(["mash", "sketch", "-k", str(self.kmer), "-s", str(self.hashes), g.genome], stderr = FNULL)

    def mash_compare(self, genome1, genome2):
        FNULL = open(os.devnull, 'w')
        cmd = Popen([self.executable, "dist", genome1.mash, genome2.mash], stderr = FNULL, stdout=PIPE)
        out = cmd.stdout.read().replace(" ","_").split()
        out_dict = {}
        out_dict['dist'] = float(out[2])
        out_dict['pvalue'] = float(out[3])
        out_dict['counts'] = int(out[4].split("/")[0])
        out_dict['nb_kmers'] = int(out[4].split("/")[1])
        return out_dict

    def mash_compare_many(self, genome1, genomes, cpus=8):
        genomes_string = [g.mash for g in genomes]
        genome_map = {g.genome : g for g in genomes}
        FNULL = open(os.devnull, 'w')
        cmd = Popen([self.executable, "dist", "-p", str(cpus) , genome1.mash] + genomes_string, stderr = FNULL, stdout=PIPE)
        out_lines = cmd.stdout.readlines()

        out_dict = {}
        for out in out_lines:
            out = out.replace(" ","_").split()
            key = genome_map[out[1]]
            out_dict[key] = {}
            out_dict[key]['dist'] = float(out[2])
            out_dict[key]['pvalue'] = float(out[3])
            out_dict[key]['counts'] = int(out[4].split("/")[0])
            out_dict[key]['nb_kmers'] = int(out[4].split("/")[1])
        return out_dict
