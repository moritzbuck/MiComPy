import os
from subprocess import call, Popen, PIPE
from Bio import SeqIO
from os.path import join as pjoin
from os import makedirs
import json

class Genome(object):
    def __repr__(self):
        return "a Genome named " + self.name
    
    def __init__(self, name , path, ref=None, manual_metadata=None):
        """
        name : name of genome
        path : where all data gets dumped
        ref : original fasta file
        manual metadata : a dictionary with diverse metadata
        """

        self.name = name
        self.path = path

        if ref:
            self.ref = ref
            
        self.size = None
        self.gc = None
        self.metadata = None
        self.cluster = None
        self.checkm_meta = None        

        self.json = pjoin(self.path, name + ".json")
        if not os.path.exists(self.path):
            makedirs(self.path)

        if os.path.exists(self.json):
            self.load_data()
            
        if not self.metadata:
            self.metadata = {}
        if manual_metadata:
            self.metadata.update(manual_metadata)
        self.proteom = pjoin(self.path, name + ".faa")
        self.genome = pjoin(self.path, name + ".fna")
        self.mash = pjoin(self.path, name + ".msh")
        self.short_proteom = pjoin(self.path, self.metadata['short_name'] + ".faa") if self.metadata['short_name'] == self.metadata['short_name'] else None
    
    def prokka(self, cpus=1, sequence = None):
        FNULL = open(os.devnull, 'w')
        return call(["prokka", "--centre", "X", "--outdir", self.path, "--force",  "--prefix" , self.name, "--locustag", self.name, "--cpus", str(cpus), sequence if sequence else self.ref], stderr = FNULL)

    def compute_mash(self, cpus=1, k = 21, mash_len=10000 ):
        FNULL = open(os.devnull, 'w')
        return call(["mash", "sketch", "-o", self.mash, "-p", str(cpus), "-k", str(k), "-s",str(mash_len), self.genome], stderr = FNULL)

    def mash_compare(self,g, cpus=1, ):
        FNULL = open(os.devnull, 'w')
        cmd = Popen(["mash", "dist", "-p", str(cpus), self.mash, g.mash], stderr = FNULL, stdout=PIPE)
        out = cmd.stdout.read().split()
        out_dict = {}
        out_dict['dist'] = float(out[2])
        out_dict['pvalue'] = float(out[3])
        out_dict['counts'] = int(out[4].split("/")[0])
        out_dict['nb_kmers'] = int(out[4].split("/")[1])
        return out_dict
    
    def is_annotated(self):
        return os.path.exists(self.proteom)

    def compute_size(self):
        with open(self.ref) as handle:
            self.size = sum([len(s) for s in SeqIO.parse(handle, "fasta")])
        self.write_data()

    def compute_gc(self):
        with open(self.ref) as handle:
            self.gc = sum([s.seq.count('G')+s.seq.count('C') for s in SeqIO.parse(handle, "fasta")])
            self.gc = float(self.gc)/self.size
        self.write_data()

        
    def get_meta(self, field, default = None):
        if self.metadata.has_key(field):
            return self.metadata[field] if self.metadata[field] == self.metadata[field] else default
        else :
            return -1

        
    def write_data(self):
        temp = {}
        temp['metadata'] = self.metadata
        temp['size'] = self.size
        temp['gc'] = self.gc
        temp['cluster'] = self.cluster
        temp['checkm'] = self.checkm_meta 
        with open(self.json,"w") as handle:
            json.dump(temp, handle)

    def load_data(self):
        with open(self.json) as handle:
            temp = json.load(handle)
        self.metadata = temp['metadata']
        self.size = temp['size']
#        self.gc = temp['gc']
        self.cluster = temp['cluster']
        self.checkm_meta = temp['checkm']
        
    def completness(self):
        return self.checkm_meta['Completeness']/100.0 if self.checkm_meta.has_key('Completeness') else -1

    def contamination(self):
        return self.checkm_meta['Contamination']/100.0 if self.checkm_meta.has_key('Completeness') else -1

    def is_good(self, min_comp = 0.35, max_cont = 0.05):
        return (self.completness() > min_comp and self.contamination() < max_cont and self.metadata['phylum'] != "Chloroflexi" and self.metadata['phylum'] != "Elusimicrobia") or self.metadata['type'] == "SAG" 

    def conv_name(self):
        return  self.name if self.metadata['short_name'] != self.metadata['short_name'] else self.metadata['short_name']
            
