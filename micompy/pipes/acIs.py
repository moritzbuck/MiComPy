#!/home/moritz/.pyenv/shims/python
#SBATCH -D /home/moritz/repos/MiComPy/
#SBATCH -J acIs
#SBATCH -o /home/moritz/acI_clustering.out
#SBATCH -e /home/moritz/acI_clustering.err
#SBATCH -A b2014036
#SBATCH -t 5-00:00:00
#SBATCH -n 16
#SBATCH -p core
#SBATCH --mail-user murumbii@gmail.com
#SBATCH --mail-type=ALL

import os
from pandas import DataFrame
from os.path import join as pjoin
from subprocess import call
import sys
from tqdm import tqdm
import sh
from micompy.common.genome import Genome
from micompy.pipes.analyses import *
from micompy.common.utils.renaming_tree import renaming_tree
from micompy.common.utils.intrasimilarity import NIC_similarity
from micompy.gene_clusterings.orthomcl.orthoMCL import orthoMCL
from micompy.gene_clusterings.orthomcl.clustering import Clustering as MCLClustering
from micompy.gene_clusterings.clustering import Clustering
from micompy.gene_clusterings.pfam_clusters.clustering import PfamClustering
from itertools import groupby
from pylab import *
from micompy.common.utils.iotl_annotations import *



root = "/home/moritz/people/sarahi/all_enrichmentss/"
data_root = pjoin(root, "all_AGs/")
analyses_root = pjoin(root, "")
google_data = pjoin(root, "ag_metadata.csv")
manual_metadata = DataFrame.from_csv(google_data).transpose().to_dict()
cpus = 16


all_genomes = [ Genome(g, pjoin(data_root,g), pjoin(data_root, m['genomes']), manual_metadata[g]) for g,m in manual_metadata.iteritems()]

for g in tqdm(all_genomes):
    if not g.size:
        g.compute_size()

all_genomes.sort(key=lambda x: x.size, reverse = True)

annotation(all_genomes, cpus)
sh.cat(*[g.proteom.replace(".faa",".gff") for g in all_genomes],_out ="temp.gff")
sh.grep("CDS","temp.gff",_out = pjoin(analyses_root,"all_gff.gff"))

#checkm(all_genomes, pjoin(analyses_root,"checkm"), cpus)

mcl = orthoMCL(pjoin(analyses_root, "orthoMCL/"), all_genomes, "big_clustering")
#mcl.start_server()
#mcl.full_pipe()
#mcl.stop_server()
#mcl.post_blast()
clusters = Clustering(all_genomes, pjoin(analyses_root, "clustering/"),"acIs", checkm = pjoin(analyses_root,"checkm"), gff = pjoin(analyses_root,"all_gff.gff"))

#cooc = clusters.cooccurence_matrix()
#cooc.to_csv(pjoin(clusters.path,"coocurence.txt"))
#bmft = clusters.make_cluster_bmft()
#bmft.to_csv(pjoin(clusters.path,"bmft.txt"))
for c in clusters:                                                                                                               
    c.genomes = set(c.genomes)
    
#intersects = {c1.name : {c2 : 2.0*float(len(c1.genomes.intersection(c2.genomes)))/(len(c1.genomes)+ len(c2.genomes)) for c2 in clusters} for c1 in clusters}
intersects = {c1.name : {c2 : float(len(c1.genomes.intersection(c2.genomes)))/min(len(c1.genomes),len(c2.genomes)) for c2 in clusters if len(c2.genomes) > 1} for c1 in clusters if len(c1.genomes) > 1}
DataFrame.from_dict(intersects).to_csv(pjoin(clusters.path,"cluster_coocs.txt"))

with open(pjoin(analyses_root,"soft_core.txt")) as handle:
    softcore =  [c[:-1] for c in handle.readlines()]

with open(pjoin(analyses_root,"hard_core.txt")) as handle:
    hardcore =  [c[:-1] for c in handle.readlines()]
    
with open(pjoin(analyses_root, "cluster_denses.txt"),"w") as handle:
    handle.writelines([str(sum([len(c.genomes) == i for c in clusters])) +"\n" for i in range(1,43) ])

with open(pjoin(analyses_root, "hypo_prob.txt"),"w") as handle:
    handle.writelines([str(float(sum([c.annotation == "hypothetical protein" for c in clusters if len(c.genomes) == i]))/sum([len(c.genomes) == i for c in clusters])) + "\n" for i in range(1,43)])

    
