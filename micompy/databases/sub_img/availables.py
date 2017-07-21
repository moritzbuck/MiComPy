from pandas import DataFrame
from Bio import SeqIO
from pandas import Index
from ete2 import NCBITaxa

data_path = "/home/moritz/people/MoreData/genomes/img_od1s"
img_fasta = "/home/moritz/people/MoreData/raw_imgs/od1s.fasta"
img_xls = "/home/moritz/people/MoreData/raw_imgs/od1s.xls"
name = "parcu_from_img_"
taxDb = NCBITaxa()
     
contigs = DataFrame.from_csv(img_xls,sep="\t",header=0,index_col=0)
manual_taxo = taxDb.get_name_translator(['Candidatus Parcubacteria']).values()[0][0]
metadata = {name + str(g) : { 'IMG_ID' : g , 'name' : name + str(g), 'species_taxid' : manual_taxo, 'long_name' : contigs.loc[contigs['Genome ID'] == g]['Genome'].iloc[0] } for g in set(contigs['Genome ID'])}

seq_dict = {k : [] for k in metadata}

with open(img_fasta,"r") as file:
    for i,c in enumerate(SeqIO.parse(file,"fasta")):
        seq_dict[name + str(contigs.iloc[i]['Genome ID']) ]+= [c]


    


