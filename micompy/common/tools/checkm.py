import subprocess
import os
import shutil
from subprocess import call, check_output
import json
from os.path import join as pjoin
from micompy.common.tools.tool import Tool

class CheckM(Tool):
    checkm_fields = ['Bin Id',
    'Marker lineage',
    '# genomes',
    '# markers',
    '# marker sets',
    '0',
    '1',
    '2',
    '3',
    '4',
    '5+',
    'Completeness',
    'Contamination',
    'Strain heterogeneity']

    checkm_marker_sets = {
    'Bacteria' : "/home/moritz/DataBases/checkm/mksets/bacteria.ms"
    }

    def __init__(self, executable = "checkm" ):
        Tool.__init__(self, name = "CheckM", executable = executable)

    def run_checkm(self, g, marker_set = "Bacteria", nproc = 1):
        infolder = g.path
        temp_dir = pjoin(g.path,  "checkm_temp")

        if os.path.exists(temp_dir):
            shutil.rmtree(temp_dir)

        assert len([f for f in os.listdir(infolder) if f[-4:]==".fna"]) == 1, "More than one (or no) genome file in %s (finishing in .fna)" % (infolder)

        with open(os.devnull, 'w') as FNULL:
            out = check_output(["checkm", "analyze", "-g", "-x", "faa", "-t", str(nproc), self.checkm_marker_sets[marker_set], infolder, temp_dir], stderr = FNULL)
            dat = check_output(["checkm", "qa", self.checkm_marker_sets[marker_set], temp_dir], stderr = FNULL)

        shutil.rmtree(temp_dir)

        dat = {k : v for k,v in zip(self.checkm_fields, dat.split("\n")[3].split())}

        with open(pjoin(infolder, g.name + '.checkm.json'), 'w') as outfile:
            json.dump(dat, outfile)

        return dat
