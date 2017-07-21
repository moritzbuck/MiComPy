from os.path import join as pjoin
import os
from tqdm import tqdm
from pandas import DataFrame, Index
import json
from micompy.common.genome import Genome
from micompy.databases.database import Database
from openpyxl import load_workbook
import shutil
from Bio import SeqIO

class SubImg(Database) :

    
    def __init__(self, workbench, data_path, name = "from_img_" , manual_taxo = None, img_xls = None, img_fasta = None, clean = False):
        Database.__init__(self,workbench = workbench, data_path = data_path)



        if img_xls and img_fasta:
            print "Extracting genomes from img files" 
            self.metadata = DataFrame.from_csv(img_xls,sep="\t",header=0,index_col=0)
            self.metadata['Genome ID'] = self.metadata['Genome ID'].apply(str)
            contigs = DataFrame.from_csv(img_xls,sep="\t",header=0,index_col=0)
            self.metadata = {name + str(g) : { 'IMG_ID' : g , 'name' : name + str(g), 'species_taxid' : manual_taxo, 'long_name' : contigs.loc[contigs['Genome ID'] == g]['Genome'].iloc[0] } for g in set(contigs['Genome ID'])}

            seq_dict = {k : [] for k in self.metadata}
            
            with open(img_fasta,"r") as file:
                for i,c in enumerate(SeqIO.parse(file,"fasta")):
                    seq_dict[name + str(contigs.iloc[i]['Genome ID']) ]+= [c]

            for k, seqs in seq_dict:
                genome_path = pjoin(self.data_path, v['region'], k)
                genome_file = pjoin(genome_path, k + ".fna")
                if not os.path.exists(genome_file):
                    os.makedirs(pjoin(genome_path, 'original_files'))
                with open(dir +  clean_g + ".fasta","w") as file:
                    SeqIO.write(seqs, pjoin(genome_path, 'original_files', k + ".fna"), "fasta")
                self.genomes += [Genome(k, genome_path, ref=pjoin(genome_path, 'original_files', k + ".fna"), manual_metadata = v, taxDb = self.taxDb, workbench = self.workbench)]
            self.metadata.to_csv(self.metadata_file
        print "Loading genomes"


        for k,v in tqdm(self.metadata.items()):
            genome_path = pjoin(self.data_path, v['region'], k)
            genome_file = pjoin(genome_path, k + ".fna")
            with open(dir +  clean_g + ".fasta","w") as file:
                SeqIO.write(seq_sets[g],file,"fasta")

            if not os.path.exists(genome_file):
                os.makedirs(pjoin(genome_path, 'original_files'))
                shutil.move(self.data_path + k + ".fna", pjoin(genome_path, 'original_files'))
            self.genomes += [Genome(k, genome_path, ref=pjoin(genome_path, 'original_files', k + ".fna"), manual_metadata = v, taxDb = self.taxDb, workbench = self.workbench)]

    def process(self, num_cores = 10):
        Database.process(self)




        

