#!/home/moritz/.pyenv/shims/python
#SBATCH -D /home/moritz/people/moritz/CDs/
#SBATCH -J OD1s
#SBATCH -o /home/moritz/people/moritz/CDs/run2.out
#SBATCH -e /home/moritz/people/moritz/CDs/run2.err
#SBATCH -A b2014204
#SBATCH -t 3-00:00:00
#SBATCH -n 4
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
from micompy.gene_clusterings.orthomcl.clustering import Clustering



root = "/home/moritz/people/moritz/CDs/"
data_root = pjoin(root, "data/")
analyses_root = pjoin(root, "analyses/")
google_data = pjoin(root, "OD1s and more - Sheet1.csv")
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

all_genomes.sort(key=lambda x: x.size, reverse = True)
all_clusters = [g for g in all_genomes if g.name == g.cluster]
short_proteoms = [g.short_proteom if g.short_proteom else g.proteom for g in all_genomes if g.is_good()]
proteoms = [g.proteom for g in all_genomes if g.is_good()]
name_map = {g.name  : g.conv_name() for g in all_genomes if g.name  !=  g.conv_name() }
rev_name_map = {v:k for k,v in name_map.iteritems()}
name_map.update({g.conv_name : g.conv_name for g in all_genomes })
rev_name_map.update({g.name : g.name for g in all_genomes })

mcl = orthoMCL(pjoin(analyses_root, "orthoMCL/"), short_proteoms, "big_clustering")
clusters = Clustering(proteoms, pjoin(analyses_root, "clustering/"),"candidate_divs", mcl, checkm = pjoin(analyses_root,"checkm"),  name_map = name_map, rev_name_map = rev_name_map)
taxo_map  = {g.name : g.name + ("--" + str(g.metadata['taxonomy_external']) if g.metadata['taxonomy_external'] == g.metadata['taxonomy_external'] else "") for g in all_genomes}
g2clusters = {}

for g in all_genomes:
    if g.is_good():
        g2clusters[g.conv_name()] = []

for c in clusters.single_copy_clusters():
    if len(c.genomes) > 10:
        for cc in c.genome_2_gene_map:
            g2clusters[name_map[cc]] += [ c.name ]

for c in tqdm(clusters.single_copy_clusters()):
    if len(c.genomes) > 8:
        c.coreness=c.compute_coreness()

pfams = PfamClustering(proteoms, pjoin(analyses_root, "pfam_clustering/"), "pfam_cdivs", checkm=pjoin(analyses_root,"checkm"))

pfams.hmm_dict = pfams.parse_hmmer_results("/home/moritz/people/moritz/CDs/hmmer_good_prots.raw")

sc_pfam_clusts = [c for c in pfams if c.name in sc_pfams]

        
#if __name__ == '__main__':
#    
#    class Unbuffered(object):
#        def __init__(self, stream):
#            self.stream = stream
#        def write(self, data):
#            self.stream.write(data)
#            self.stream.flush()
#        def __getattr__(self, attr):
#            return getattr(self.stream, attr)
#
#    sys.stdout = Unbuffered(sys.stdout)
#    
##   annotation(all_genomes)
##    cluster_genomes(all_genomes,pjoin(analyses_root,"rifle_clusters.tsv"),cutoff=0.95)
##    phylophlan(all_clusters, pjoin(analyses_root), proj_name = "all_od1s_plus", cpus=16)
##    for g in all_genomes:
##        if g.metadata['short_name'] == g.metadata['short_name']:
##            os.system("sed 's/%s/%s/g' %s > %s" %(g.name, g.metadata['short_name'], g.proteom, pjoin(g.path, g.metadata['short_name'] + ".faa")))
##    checkm(all_genomes,pjoin(analyses_root,"checkm"), cpus=16)
##    for g in all_genomes:
##        parse_checkm_results(g, pjoin(analyses_root, "checkm"))
##    mcl.start_server()
##    mcl.full_pipe()
##    mcl.stop_server()
##    clusters.post_process()
    all_trees = []
    ccs = []
#    threash = len([g for g in all_genomes if g.is_good()])*0.05
#    short_names = {g.name : g.metadata['short_name'] for g in all_genomes if g.is_good()}
    for c in tqdm(clusters.single_copy_clusters()):
        if len(c.genomes) > 3 and c.coreness < 0 :
            ccs += [c.name]
#            all_trees += [pjoin(clusters.path, "clusters/", c.name, "tree",  c.name + ".tree")] 
##            if not os.path.exists(pjoin(clusters.path, "clusters/", c.name, "align")):
##                os.makedirs(pjoin(clusters.path, "clusters/", c.name, "align"))
##            if not os.path.exists(pjoin(clusters.path, "clusters/", c.name, "tree")):
##                os.makedirs(pjoin(clusters.path, "clusters/", c.name, "tree"))
##            c.align(pjoin(clusters.path, "clusters/", c.name, "align", c.name + "_aligned.faa"), block=False, genome_names = True)
##            if os.path.exists(pjoin(clusters.path, "clusters/", c.name, "align", c.name + "_aligned_blocked.faa")):
##                c.tree_construction(pjoin(clusters.path, "clusters/", c.name, "align", c.name + "_aligned_blocked.faa"), pjoin(clusters.path, "clusters/", c.name, "tree", c.name + ".tree"))
##                renaming_tree(pjoin(clusters.path, "clusters/", c.name, "tree", c.name + ".tree"), pjoin(clusters.path, "clusters/", c.name, "tree", "short_" + c.name + ".tree"), short_names)

#    sh.cat(*all_trees, _out = pjoin(clusters.path , "tree_set.tree"))
#    concat_core_tree(clusters, pjoin(analyses_root, "tree_trunk"))
#    
#"""
#prset aamodelpr = mixed
#mcmc nchains = 1 ngen = 300000
#
#"""
