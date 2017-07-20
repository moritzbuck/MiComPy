import subprocess
import os
import shutil
from subprocess import call, check_output, Popen, PIPE
import json
from os.path import join as pjoin
from micompy.common.tools.tool import Tool
import tempfile
from pandas import read_csv

class HMMer(Tool):
    col_names = ["target_name" , "target_accession" , "tlen" , "query_name" , "query_accession" , "qlen" , "E-value","score" , "bias" , "nb" , "of" , "c-Evalue" , "i-Evalue" , "score" , "bias" , "hmm_from" , "hmm_to" , "ali_from" , "ali_to" , "env_from" , "env_to" , "acc" , "description_of_target"]

    def __init__(self, executable = "hmmsearch"):
        Tool.__init__(self, name = "HMMer", executable = executable)


    def hmmsearch_pfam_presence(self, genome, pfams = "/home/moritz/DataBases/PFAM/Pfam-A.hmm", ncpus = 1, evalue_cutoff = 0.001):
        temp_file = tempfile.NamedTemporaryFile(suffix='.hmmsearch.out', delete = False).name
        print "Running HMMsearch with pfams on", genome
        with open(os.devnull, 'w') as handle:
            call([self.executable, "--domtblout", temp_file, "--cpu", str(ncpus) , pfams, genome.proteom], stdout=handle, stderr=handle)
        output = self.parse_hmmer_tabout(temp_file)
        output = set(output[output['i-Evalue'] < evalue_cutoff]['query_accession'])
        os.remove(temp_file)
        return output

    def parse_hmmer_tabout(self, file):
        table = read_csv(file, delim_whitespace=True, comment="#", names = HMMer.col_names)
        return table
