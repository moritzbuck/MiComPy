import os
from subprocess import call, Popen, PIPE
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from os.path import join as pjoin
from os import makedirs
import json
from numpy import random as nrandom
from numpy import nan
from random import random
from tqdm import tqdm

class Genome(object):
    def __repr__(self):
        return "a Genome named " + self.name

    def __init__(self, name , path, ref=None, manual_metadata=None, taxDb = None, workbench = None):
        """
        name : name of genome
        path : where all data gets dumped
        ref : original fasta file
        manual metadata : a dictionary with diverse metadata
        """

        self.name = name
        self.path = path
        self.taxDb = taxDb
        self.workbench = workbench

        if ref:
            self.ref = ref

        self.size = None
        self.gc = None
        self.metadata = None
        self.cluster = None
        self.checkm_meta = None
        self.pfams = None
        self.fakereads = pjoin(self.path, self.name + ".fakereads.fasta")

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
        self.mash = pjoin(self.path, name + ".fna.msh")
        self.short_proteom = pjoin(self.path, self.metadata['short_name'] + ".faa") if self.metadata.has_key('short_name') and self.metadata['short_name'] == self.metadata['short_name'] else None
        self.taxo = None


    def clean(self):
        def clean_rm(f) :
            if os.path.exists(f):
                os.remove(f)
        clean_rm(self.fakereads)
        clean_rm(self.json)
        clean_rm(self.mash)
        clean_rm(pjoin(self.path,"ref"))

    def prokka(self, cpus=1, sequence = None):
        ############### PUT IN TOOOOOOOOOOOOOOOOOOOLS 
        FNULL = open(os.devnull, 'w')
        return call(["prokka", "--centre", "X", "--outdir", self.path, "--force",  "--prefix" , self.name, "--locustag", self.name, "--cpus", str(cpus), sequence if sequence else self.ref], stderr = FNULL)
        close(FNULL)

    def compute_checkm(self, cpus=1):
        checker = self.workbench.tools['CheckM']
        self.checkm_meta = checker.run_checkm(self)
        self.write_data()

    def compute_mash_hash(self):
        masher = self.workbench.tools['MASH']
        masher.run_mash_sketch(self)

    def get_pfams(self, evalue_cutoff = 0.001):
        if not self.pfams:
            hmmer = self.workbench['HMMer']
            self.pfams = hmmer.hmmsearch_pfam_presence(self, evalue_cutoff = evalue_cutoff)
            self.write_data()
        return self.pfams

    def mash_compare(self,g, cpus=1, ):
        masher = self.workbench.tools['MASH']
        out_dict = masher.mash_compare(self,g)
        return out_dict

    def mash_compare_many(self,g, cpus=8, ):
        masher = self.workbench.tools['MASH']
        out_dict = masher.mash_compare_many(self,g, cpus)
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

    def clean_proteom(self):
        if not os.path.exists(pjoin(self.path, "original_files", "protein_name_map.csv")):
            name_map = {}

            with open(self.proteom) as handle:
                all_prots = [p for p in SeqIO.parse(handle,"fasta")]


            for i,p in enumerate(all_prots):
                name_map[p.id] = self.name + "_" + str(i).zfill(5)
                p.id =  name_map[p.id]
                p.description = ""

            with open(pjoin(self.path, "original_files", "protein_name_map.csv"), "w") as handle:
                handle.writelines([k + "," + v + "\n" for k, v in name_map.items()])

            with open(self.proteom, "w") as handle:
                SeqIO.write(all_prots, handle, "fasta")

    def write_data(self):
        temp = {}
        temp['metadata'] = self.metadata
        temp['size'] = self.size
        temp['gc'] = self.gc
        temp['cluster'] = self.cluster
        temp['checkm'] = self.checkm_meta
        temp['pfams'] = list(self.pfams)

        with open(self.json,"w") as handle:
            json.dump(temp, handle)

    def load_data(self):
        with open(self.json) as handle:
            temp = json.load(handle)
        self.metadata = temp.get('metadata')
        self.size = temp.get('size')
        self.gc = temp.get('gc')
        self.cluster = temp.get('cluster')
        self.checkm_meta = temp.get('checkm')
        self.pfams = temp.get('pfams')
        if self.pfams:
            self.pfams = set(self.pfams)

    def completness(self):
        return float(self.checkm_meta['Completeness'])/100.0 if self.checkm_meta.has_key('Completeness') else -1

    def contamination(self):
        return float(self.checkm_meta['Contamination'])/100.0 if self.checkm_meta.has_key('Completeness') else -1

    def is_good(self, min_comp = 0.35, max_cont = 0.05):
        return self.completness() > min_comp and self.contamination() < max_cont

    def conv_name(self):
        return  self.name if not self.metadata.has_key('short_name') or self.metadata['short_name'] != self.metadata['short_name'] else self.metadata['short_name']

    def get_sequence(self):
        with open(self.genome) as handle:
            seqs = [s for s in SeqIO.parse(handle, "fasta")]
        return seqs

    def get_ANI(self,genome):
        mapper = self.workbench.tools['BBMap']
        if not os.path.exists(genome.fakereads):
            genome.make_fake_reads()

        if not os.path.exists(pjoin(self.path,"ref")):
            mapper.make_index(self)

        return mapper.get_ANI(self,genome.fakereads)

    def make_bbmap_index(self):
        mapper = bbmap()
        mapper.make_index(self)


    def make_fake_reads(self, read_len = 150):
        seqs = self.get_sequence()
        split_seq = lambda seq, x : [seq[(i*x):(i*x+x)]for i in  range(len(seq)/x)]

        reads = sum([ split_seq(str(s.seq),read_len) for s in seqs if len(s) > read_len], [])
        zz = len(str(len(reads)))

        with open(self.fakereads, "w") as handle:
            handle.writelines([">" + self.name + "_fakeread_" + str(i).zfill(zz) + "\n" + r + "\n" for i,r in enumerate(reads)])

    def get_taxo(self):
        if self.taxo:
            return self.taxo
        elif self.metadata.get('species_taxid'):
            taxid = self.metadata.get('species_taxid')
            rank = self.taxDb.get_rank(self.taxDb.get_lineage(taxid))
            taxa = self.taxDb.get_taxid_translator(self.taxDb.get_lineage(taxid))
            self.taxo = {rank[k] : taxa[k] for k in rank}
            return self.taxo
        else :
            return None

    def make_simu_mapping(self, rate, type = "full"):
# type full : fully random position of murations, rate is the mutation rate per base
# type uniform : mutation at regular interval, rate is the spacer between each mutation
# type core : like full but a proportion of the genome is left unchanged, rate[0] is the mut rate, rate[1] the proportion left unchanged

        nucs = {'A', 'T', 'C', 'G'}


        seq = "".join([str(s.seq) for s in self.get_sequence()])
        GC = float(seq.count('G')+seq.count('C'))/len(seq)
        base_rates = {'A' : (1-GC)/2, 'T' : (1-GC)/2, 'G' : GC/2, 'C' : GC/2}
        sub_rates = {b : {bb : base_rates[bb]/(1-base_rates[b]) for bb in nucs.difference(b)} for b in nucs}
        sub_rates = {b : (rates.keys() , rates.values()) for b, rates in sub_rates.items()}
        if type == "full":
            mutated = [s if random() < rate or not s in nucs else nrandom.choice(sub_rates[s][0], p=sub_rates[s][1]) for s in tqdm(seq)]
        elif type == "uniform":
            mutated = [s if i % rate != 0  or not s in nucs else nrandom.choice(sub_rates[s][0], p=sub_rates[s][1]) for i,s in tqdm(enumerate(seq))]
        elif type == "core" :
            core_len = int(len(seq) * rate[1])
            core = seq[:core_len]
            aux = seq[core_len:]
            mutated = [s if random() < rate[0]*(1-rate[1]) or not s in nucs else nrandom.choice(sub_rates[s][0], p=sub_rates[s][1]) for s in tqdm(aux)]
            mutated = [s for s in core] + mutated

        with open("/tmp/fake.fna", "w") as handle :
            SeqIO.write(SeqRecord(Seq("".join(mutated)), id = "fake"), handle, "fasta")

        fake = Genome("fake", "/tmp/", "fake.fna", manual_metadata = {'short_name' : nan})

        mapping = self.get_ANI(fake)
        fake.clean()

        return mapping
