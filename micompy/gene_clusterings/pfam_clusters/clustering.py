from micompy.gene_clusterings.clustering import Clustering
from micompy.gene_clusterings.gene_cluster import GeneCluster
from subprocess import call
import shutil
from micompy.common.utils.params import cpus, pfama_loc, temp_loc, pjoin
from tqdm import tqdm 

class PfamClustering(Clustering):
	def __init__(self,proteoms, out_path, name, gff = None, seq_type="proteins", checkm = None, name_map = None):
		Clustering.__init__(self,proteoms, out_path, name, gff = None, seq_type="proteins", checkm = None, name_map = None)

	def run(self):
		call(" ".join(["cat"] + self.proteoms + [ " > " + pjoin(temp_loc, "full_proteom.faa")]))
		shutil.copy(pfama_loc, temp_loc)
		call(["hmmsearch --cpu 16 -o", pjoin(temp_loc,"hmmer_good_prots.raw"), pjoin(temp_loc,"Pfam-A.hmm"), pjoin(temp_loc, "full_proteom.faa") ])
		self.hmm_dict = self.parse_hmmer_results(pjoin(temp_loc,"hmmer_good_prots.raw"))

	def parse_hmmer_results(self, file):

		entries = []
		t = []
		post = False
		with open(file) as handle:
			print "Parsing hmmsearch output!"
			for l in tqdm(handle.readlines()):
				if l[0:17] == "Domain annotation" and not post :
					post = True
				if l[0] != "#" and len(l[:-1]) > 0 and not post :
					t.append(l[:-1])
				if l == "//\n":
					entries.append(t)
					t = []
					post = False 

		hmm_dict = {}
		print "Making a nice usable data structure!"
		for e in tqdm(entries):
			acc = e[1].split()[1].split(".")[0]
			hmm_dict[acc] = {}
			hmm_dict[acc]['version'] = "NA"
			if "KO:" in e[1].split()[1]:
				if len(e[2].split()) > 1:
					hmm_dict[acc]['version'] = " ".join(e[2].split()[2:])
			if "PF" in e[1].split()[1]:
				hmm_dict[acc]['version'] = e[1].split()[1].split(".")[1]
				
			hmm_dict[acc]['short_name'] = e[0].split()[1]
			hmm_dict[acc]['description'] = " ".join(e[2].split()[1:])
			hmm_dict[acc]['cdss'] = []
			for l in e[7:]:
				if "inclusion" in l or "[No hits detected that satisfy reporting thresholds]" in l :
					break 
				hmm_dict[acc]['cdss'].append(l.split()[8])
		return hmm_dict

	def post_process(self):
		print "Post processing PFAM cluster:"
		self.clusters=[GeneCluster(self, name=k, genes =   h['cdss']) for k,h in tqdm(self.hmm_dict.iteritems())]
		print "Post processing single genes:"
		non_singletons = set(sum([c.genes for c in self.clusters],[]))
		for i in tqdm(self.id2name_map.keys()):
			if i not in non_singletons:
				self.clusters += [GeneCluster(self, name = i,  genes =  [self.gene2genome[i] + "|" + i ])]
			
		with open(self.processed_clusters, 'w') as outfile:
			json.dump([c.to_dict() for c in self.clusters], outfile,  indent=4, sort_keys=True)

