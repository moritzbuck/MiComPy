import sys, os
from subprocess import check_output, call
from os.path import join as pjoin

core = False

def printerr(l):
	print('\033[93m [ERR] ' + l + '\033[0m')

if __name__ == '__main__':  
    class Unbuffered(object):
        def __init__(self, stream):
            self.stream = stream
        def write(self, data):
            self.stream.write(data)
            self.stream.flush()
        def __getattr__(self, attr):
            return getattr(self.stream, attr)
    sys.stdout = Unbuffered(sys.stdout)

devnull = open(os.devnull, 'w')


print "Checking if all python packages are there"
try : 
	from pandas import DataFrame
	from tqdm import tqdm
	from Bio import SeqIO
	print "Required packages loaded"
except :
	printerr("Some python packages are wrong")


print "Checking if all internal stuff is running"
try :
	from micompy.common.genome import Genome
	from micompy.pipes.analyses import *
	from micompy.common.utils.renaming_tree import renaming_tree
	from micompy.common.utils.intrasimilarity import NIC_similarity
	from micompy.gene_clusterings.orthomcl.orthoMCL import orthoMCL
	from micompy.gene_clusterings.orthomcl.clustering import Clustering
except :
	printerr("Internal imports are bad")


print "Checking if all exes there"
#try:
exes = ["find", "blastn", "makeblastdb", "prokka", "hmmsearch", "cat"]
if not all([call(["which",cmd], stdout = devnull) == 0 for cmd in exes]) :
	printerr("Missing executables : ")
	printerr(",".join([cmd for cmd in exes if call(["which",cmd], stdout = devnull) == 1]))

print "Loading  test genomes"

root = "test_data/"
data_root = pjoin(root, "data/")
analyses_root = pjoin(root, "analyses/")
google_data = pjoin(root, "metadata.csv")
manual_metadata = DataFrame.from_csv(google_data).transpose().to_dict()
cpus = 1
all_genomes = []


all_files = check_output(["find" , data_root]).split()
for g in manual_metadata.keys():
    fasta = [f for f in all_files if ".fasta" in f and g in f and ".fasta." not in g]
    assert len(fasta) == 1
    all_genomes += [Genome(g, os.path.dirname(fasta[0]), fasta[0], manual_metadata[g])]


if not core:
	try :
		print "testing genome clustering with NICsimilarity"
#		cluster_genomes(all_genomes,pjoin(analyses_root,"rifle_clusters.tsv"),cutoff=0.95)
		print "testing annotation"
		annotation(all_genomes, cpus = 8 , )
	except : 
		printerr("non-core functions broken")
	print "non-core functions work"




print "All went well!"