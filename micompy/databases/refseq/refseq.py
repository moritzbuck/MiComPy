from ftplib import FTP
import StringIO
from os.path import join as pjoin
import os
from tqdm import tqdm
from pandas import DataFrame
from joblib import Parallel, delayed
import multiprocessing
import shutil
from subprocess import call, check_output
import json
from micompy.common.genome import Genome
from micompy.databases.database import Database
from micompy.pipes.analyses import *
from numpy import nan
import signal

ncbi = "ftp.ncbi.nlm.nih.gov"
def download_ftp(genome_path, genome_file, fhead, dir):
    def handler(signum, frame):
        print "FTP download timeout"
        raise Exception("FTP Downlaod timeout")

    proteom_file = genome_file.replace(".fna", ".faa")
    signal.signal(signal.SIGALRM, handler)
    signal.alarm(60)

    try :
        ftp = FTP(ncbi)
        ftp.login()
        ftp.cwd(dir)
        ftp.retrbinary("RETR " + fhead + "_genomic.fna.gz", open(pjoin(genome_path, fhead +".fna.gz"),"w").write)
        ftp.retrbinary("RETR " + fhead + "_protein.faa.gz", open(pjoin(genome_path, fhead +".faa.gz"),"w").write)
#        ftp.retrbinary("RETR " + fhead + "_protein.faa.gz", open(pjoin(genome_path, fhead +".gff.gz"),"w").write)
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

    dpat = ftp_head.split("/")
    dir = "/".join(dpat [3:])

    fhead = dpat[-1]

    download_ftp(genome_path,genome_file, fhead, dir)

class Refseq(Database) :

    def __init__(self, workbench, data_path = "/home/moritz/DataBases/genomes/RefSeq/", clean = False):
        Database.__init__(self,workbench = workbench, data_path = data_path)

        if not os.path.exists(self.metadata_file) or clean:
            ftp =  FTP(ncbi)

            print "Getting metadata from ncbi"

            FNULL = open(os.devnull, 'w')
            ftp.login()
            ftp.cwd('genomes/refseq/bacteria/')
            info = StringIO.StringIO()
            ftp.retrbinary("RETR " + "assembly_summary.txt", info.write)
            info.seek(0)
            self.metadata = DataFrame.from_csv(info, sep="\t", header=1)
            ftp.close()
            self.metadata['assembly_level'] = self.metadata['assembly_level'].apply(lambda x: x.replace(" ","_"))
            self.metadata = self.metadata.transpose().to_dict()

            DataFrame.from_dict(self.metadata).to_csv(self.metadata_file)

        else :
            print "Loading metadata"
            self.metadata = DataFrame.from_csv(self.metadata_file).to_dict()

        print "Loading genomes"
        for k,v in tqdm(self.metadata.items()):
            genome_path = pjoin(self.data_path, v['assembly_level'].replace(" ","_"), k)
            genome_file = pjoin(genome_path, k + ".fna")
            self.genomes += [Genome(k, genome_path, ref=genome_file, manual_metadata = v, taxDb = self.taxDb, workbench = self.workbench)]

    def process(self, num_cores = 10):
        to_dl = []
        for k,v in tqdm(self.metadata.items()):
            genome_path = pjoin(self.data_path, v['assembly_level'].replace(" ","_"), k)
            genome_file = pjoin(genome_path, k + ".fna")
            if not os.path.exists(genome_path) :
                os.makedirs(genome_path)
            if not os.path.exists(genome_file) or (not os.path.exists(genome_file.replace(".fna",".faa"))):
                to_dl += [ (v['ftp_path'],genome_path, genome_file)]
        dlstuff= Parallel(n_jobs=num_cores)(delayed(download)(i) for i in tqdm(to_dl))

        Database.process(self)

        self.taxos = DataFrame.from_dict({g.name : g.get_taxo(self.taxDb) for g in self.genomes if g.metadata['assembly_level'] == "Complete_Genome"}).transpose().to_csv(pjoin(self.metadata_path, "complete_genomes_pretty_taxo.csv"))
