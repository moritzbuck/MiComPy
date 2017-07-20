from ete2 import NCBITaxa
from joblib import Parallel, delayed
from os.path import join as pjoin
import os
from tqdm import tqdm


def mash(info):
    genome = info
    genome.compute_mash_hash()

def checkm(info):
    genome = info
    genome.compute_checkm()

def prokka(info):
    genome = info
    genome.prokka()



class Database(object) :
    def __getitem__(self, key):
        if type(key) == int :
            return self.genomes[key]
        else :
            li = [g for g in genomes if g.name == key]
            assert len(li) == 1, "There are more than one genomes with the name " + key
            return li[0]


    def __init__(self,  data_path, workbench = None, genomes = [], taxDb = None):
        self.data_path = data_path
        self.workbench = workbench
        self.metadata_path = pjoin(self.data_path, "metadata")
        if not os.path.exists(self.metadata_path):
            os.makedirs(self.metadata_path)

        self.metadata_file = pjoin(self.metadata_path, "metadata.csv")

        if taxDb:
            self.taxDb = taxDb
        else :
            self.taxDb = NCBITaxa()

        self.genomes = genomes


    def process(self, num_cores = 10):

        print "Running prokka for protein annotation (excedpt if faas already provided)"
        to_prokka = [g for g in self.genomes if not os.path.exists(g.proteom)]
        prokka_stuff = Parallel(n_jobs=num_cores)(delayed(prokka)(i) for i in tqdm(to_prokka))

        to_mash = [g for g in self.genomes if not os.path.exists(g.genome + ".msh")]

        print "running mash hashing"
        mashstuff= Parallel(n_jobs=num_cores)(delayed(mash)(i) for i in tqdm(to_mash))

        print "running CheckM"
        to_check = [g for g in self.genomes if not os.path.exists(g.genome.replace(".fna",".checkm.json")) or not g.checkm_meta]
        checkmstuff= Parallel(n_jobs=num_cores)(delayed(checkm)(i) for i in tqdm(to_check))

        print "computing genome sizes"
        for g in tqdm(self.genomes):
            if not g.size:
                g.compute_size()

        print "computing gc contents"
        for g in tqdm(self.genomes):
            if not g.size:
                g.compute_gc()

        print "making fake reads"
        for g in tqdm(self.genomes):
            if not os.path.exists(g.fakereads):
                g.make_fake_reads(read_len=150)
