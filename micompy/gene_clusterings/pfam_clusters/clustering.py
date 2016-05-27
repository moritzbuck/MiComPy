from micompy.gene_clusterings.clustering import Clustering
from subprocess import call
import shutil
from micompy.common.utils.params import cpus, pfama_loc, temp_loc, pjoin

class PfamClustering(Clustering):

	def __init__(self,proteoms, out_path, name, gff = None, seq_type="proteins", checkm = None, name_map = None):

		Clustering.__init__(self,proteoms, out_path, name, gff = None, seq_type="proteins", checkm = None, name_map = None)

    def run(self):
		call(" ".join(["cat"] + self.proteoms  + [ " > " + pjoin(temp_loc, "full_proteom.faa")]) 
		shutil.copy(pfama_loc, temp_loc)   	
    	call(["hmmsearch --cpu 16 -o hmmer_good_prots.raw Pfam-A.hmm" , pjoin(temp_loc, "full_proteom.faa") ])
