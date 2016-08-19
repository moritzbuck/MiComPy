import os
from subprocess import call
from Bio import SeqIO
from os.path import join as pjoin
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
        self.metadata = None
        self.cluster = None
        self.checkm_meta = None        

        self.json = pjoin(self.path, name + ".json")

        if os.path.exists(self.json):
            self.load_data()
            
        if not self.metadata:
            self.metadata = {}
        if manual_metadata:
            self.metadata.update(manual_metadata)
        self.proteom = pjoin(self.path, name + ".faa")
        self.genome = pjoin(self.path, name + ".fna")
        self.short_proteom = pjoin(self.path, self.metadata['short_name'] + ".faa") if self.metadata['short_name'] == self.metadata['short_name'] else None
    
    def prokka(self, cpus=1, sequence = None):
        return call(["prokka", "--centre", "X", "--compliant", "--outdir", self.path, "--force",  "--prefix" , self.name, "--locustag", self.name, "--cpus", str(cpus), sequence if sequence else self.ref])

    def is_annotated(self):
        return os.path.exists(self.proteom)

    def compute_size(self):
        with open(self.ref) as handle:
            self.size = sum([len(s) for s in SeqIO.parse(handle, "fasta")])
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
        temp['cluster'] = self.cluster
        temp['checkm'] = self.checkm_meta 
        with open(self.json,"w") as handle:
            json.dump(temp, handle)

    def load_data(self):
        with open(self.json) as handle:
            temp = json.load(handle)
        self.metadata = temp['metadata']
        self.size = temp['size']
        self.cluster = temp['cluster']
        self.checkm_meta = temp['checkm']
        
    def completness(self):
        return self.checkm_meta['Completeness']/100.0 if self.checkm_meta.has_key('Completeness') else -1

    def contamination(self):
        return self.checkm_meta['Contamination']/100.0 if self.checkm_meta.has_key('Completeness') else -1

    def is_good(self, min_comp = 0.38, max_cont = 0.05):
        return (self.completness() > min_comp and self.contamination() < max_cont and self.metadata['phylum'] != "Chloroflexi" and self.metadata['phylum'] != "Elusimicrobia") or self.metadata['type'] == "SAG" 

    def conv_name(self):
        return  self.name if self.metadata['short_name'] != self.metadata['short_name'] else self.metadata['short_name']
            
