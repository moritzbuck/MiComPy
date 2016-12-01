#!/home/moritz/.pyenv/shims/python
#SBATCH -D /home/moritz/repos/MiComPy/
#SBATCH -J OD1s
#SBATCH -o /home/moritz/people/moritz/CDs/run3.out
#SBATCH -e /home/moritz/people/moritz/CDs/run3.err
#SBATCH -A b2014036
#SBATCH -t 10-00:00:00
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
from dendropy import *
import csv

root = "/home/moritz/people/moritz/CDs/"
data_root = pjoin(root, "data/")
analyses_root = pjoin(root, "analyses")
google_data = pjoin(root, "OD1s and more - Sheet1.csv")
manual_metadata = DataFrame.from_csv(google_data).transpose().to_dict()
cpus = 16

all_raws = set(os.listdir(pjoin(data_root, "raw")))
assert all([k + ".fasta" in all_raws for k in manual_metadata.keys()])


all_genomes = [ Genome(g, pjoin(data_root,"genomes",g), pjoin(data_root,'raw', g + ".fasta"), manual_metadata[g]) for g,m in manual_metadata.iteritems()]

for g in tqdm(all_genomes):
    if not g.size:
        g.compute_size()

all_genomes.sort(key=lambda x: x.completness(), reverse = True)

annotation(all_genomes, cpus)

all_genomes = [g for g in all_genomes if g.is_good()]

derep_clusters = cluster_genomes_mash(genomes = all_genomes,simi_matrix = DataFrame.from_csv(pjoin(analyses_root, "mash_matrix.csv")), output = pjoin(analyses_root, "mash_clusters.txt"), cutoff= 3500)

all_genomes = [g for g in all_genomes if g.name == g.cluster]

#mcl = orthoMCL(pjoin(analyses_root, "dereped_orthoMCL/"), all_genomes, "derep_clustering")
#mcl.start_server()
#mcl.full_pipe()
#mcl.stop_server()

# clusters = MCLClustering(all_genomes, pjoin(analyses_root, "clustering/"),"dereped_parcus", mcl, checkm = pjoin(analyses_root,"all_genomes_quals"), gff = pjoin(analyses_root,"all_gff.gff"))
# clusters.post_process()
clusters = Clustering(all_genomes, pjoin(analyses_root, "clustering/"),"dereped_parcus", checkm = pjoin(analyses_root,"all_genomes_quals"), gff = pjoin(analyses_root,"all_gff.gff"))

with open(pjoin(analyses_root, "genome2phylum.csv"),"w") as handle:
    handle.writelines([g.short_proteom.split("/")[-1].split(".")[0] + "," + (g.metadata['phylum'] if g.metadata['phylum'] == g.metadata['phylum'] else "OD1")+"\n"  for g in all_genomes])

bmft = clusters.make_cluster_bmft()

with open(pjoin(clusters.path,"core.txt")) as handle:
    core_clusters = set([l[:-1] for l in handle.readlines()])
core_clusters = [c for c in clusters if c.name in core_clusters]
used_core = [ c for c in core_clusters if float(len(c.genes))/len(c.genomes) < 1.03]

#cls = cluster_genomes(all_genomes,pjoin(analyses_root,"genome_clusters.txt"), cutoff=0.95)
#mat = mash_matrix(all_genomes, pjoin(analyses_root,"mash_matrix.csv"))
#checkm(all_genomes, pjoin(analyses_root, "all_genomes_quals"), cpus=16)
#sh.cat(*[g.proteom.replace(".faa",".gff") for g in all_genomes],_out ="temp.gff")
#sh.grep("CDS","temp.gff",_out = pjoin(analyses_root,"all_gff.gff"))

for g in all_genomes :                                                                       
    g.new_completeness = float(sum([g.get_meta('short_name') in c.genomes for c in core_clusters]))/len(core_clusters)


clust_data("GC-content", {g.get_meta('short_name'): '%6.2f' % (g.gc) for g in all_genomes},  pjoin(analyses_root, "trees", "core_tree", "gc_content.txt"))

clust_data("genome size", {g.get_meta('short_name'): '%6.2f' % (g.size/g.new_completeness/1000000) for g in all_genomes if g.metadata['type'] != 'SAG'},  pjoin(analyses_root, "trees", "core_tree", "genome_len.txt"))
