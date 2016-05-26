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


print "All went well!"