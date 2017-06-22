from ftplib import FTP
import StringIO
from os.path import join as pjoin
import os
from tqdm import tqdm
from pandas import DataFrame
from joblib import Parallel, delayed
import multiprocessing
from ete2 import NCBITaxa
import shutil
from subprocess import call, check_output
import json
from micompy.common.genome import Genome
from numpy import nan
import signal

ncbi = "ftp.ncbi.nlm.nih.gov"
def download_ftp(genome_path, fhead, dir):
    def handler(signum, frame):
        print "FTP download timeout"
        raise Exception("FTP Downlaod timeout")

    signal.signal(signal.SIGALRM, handler)
    signal.alarm(60)

    try :
        ftp = FTP(ncbi)
        ftp.login()
        ftp.cwd(dir)
        ftp.retrbinary("RETR " + fhead + "_genomic.fna.gz", open(pjoin(genome_path, fhead +".fna.gz"),"w").write)
        ftp.retrbinary("RETR " + fhead + "_protein.faa.gz", open(pjoin(genome_path, fhead +".faa.gz"),"w").write)
        ftp.close()
        if not os.path.exists(pjoin(genome_path, "original_files")) :
            os.makedirs(pjoin(genome_path, "original_files"))
        fnas = [f for f in os.listdir(genome_path) if f.endswith(".fna.gz")]
        faas = [f for f in os.listdir(genome_path) if f.endswith(".faa.gz")]
        with open(genome_file, "w") as handle:
            call(["zcat", pjoin(genome_path, *fnas)], stdout = handle)
        with open(proteom_file, "w") as handle:
            call(["zcat", pjoin(genome_path, *faas)], stdout = handle)
        for f in faas + fnas:
            shutil.move(pjoin(genome_path,f), pjoin(genome_path, "original_files") )
    except Exception, exc:
        print exc

    signal.alarm(0)

def download(info):

    ftp_head, genome_path, genome_file = info
    proteom_file = genome_file.replace(".fna", ".faa")
    dpat = ftp_head.split("/")
    dir = "/".join(dpat [3:])

    fhead = dpat[-1]

    download_ftp(genome_path,fhead, dir)




class Refseq(object) :

    def __init__(self, data_path = "/home/moritz/DataBases/genomes/RefSeq/", check = True, marker_set = 'Bacteria', mash_hashes = 10000, mash_kmer = 21):
        self.data_path = data_path
        self.marker_set = marker_set
        self.mash_hashes = mash_hashes
        self.mash_kmer = mash_kmer

        self.taxDb = NCBITaxa()

        ftp =  FTP(ncbi)

        print("getting metadata from ncbi")

        FNULL = open(os.devnull, 'w')
        ftp.login()
        ftp.cwd('genomes/refseq/bacteria/')
        info = StringIO.StringIO()
        ftp.retrbinary("RETR " + "assembly_summary.txt", info.write)
        info.seek(0)
        self.metadata = DataFrame.from_csv(info, sep="\t", header=1)
        ftp.close()
        self.metadata['short_name'] = nan
        self.metadata = self.metadata.transpose().to_dict()

        self.genomes = []

        for k,v in tqdm(self.metadata.items()):
            genome_path = pjoin(self.data_path, v['assembly_level'], k)
            genome_file = pjoin(genome_path, k + ".fna")
            self.genomes += [Genome(k, genome_path, ref=genome_file, manual_metadata = v)]

    def process(self, num_cores = 10):
        to_dl = []
        for k,v in tqdm(self.metadata.items()):
            genome_path = pjoin(self.data_path, v['assembly_level'], k)
            genome_file = pjoin(genome_path, k + ".fna")
            if not os.path.exists(genome_path) :
                os.makedirs(genome_path)
            if not os.path.exists(genome_file):
                to_dl += [ (v['ftp_path'],genome_path, genome_file)]

        dlstuff= Parallel(n_jobs=num_cores)(delayed(download)(i) for i in tqdm(to_dl))

        self.mash()

    def mash(self, clean = False):
        for g in tqdm([gg for gg in self.genomes if os.path.exists(gg.genome) and (clean or os.path.exists(gg.genome)) ]):
            g.compute_mash_hash()

    def checkm(self):
        return None
