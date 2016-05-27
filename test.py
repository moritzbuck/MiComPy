import sys


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



try : 
	from pandas import DataFrame
	from tqdm import tqdm
	from Bio import SeqIO
	print "Required packages loaded"
except :
	print "Some python packages are wrong"
	pass


try :
	from micompy.common.genome import Genome
	from micompy.pipes.analyses import *
	from micompy.common.utils.renaming_tree import renaming_tree
	from micompy.common.utils.intrasimilarity import NIC_similarity
	from micompy.gene_clusterings.orthomcl.orthoMCL import orthoMCL
	from micompy.gene_clusterings.orthomcl.clustering import Clustering
except :
	print "Internal imports are bad"
	pass


print "loading genomes"

root = "test_data/"
data_root = pjoin(root, "data/")
analyses_root = pjoin(root, "analyses/")
google_data = pjoin(root, "metadata.csv")
manual_metadata = DataFrame.from_csv(google_data).transpose().to_dict()
cpus = 1
all_genomes = []

for dir in os.listdir(data_root):
    dir = pjoin(data_root,dir)
    for g in os.listdir(dir):
        g_dir = pjoin(dir,g)
        fasta = [f for f in os.listdir(g_dir) if ".fasta" in f and not ".fasta." in f]
        assert len(fasta) == 1
        all_genomes += [Genome(g, g_dir, pjoin(g_dir,fasta[0]), manual_metadata[g])]


print "All went well!"