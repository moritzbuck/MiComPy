#!/home/moritz/.pyenv/shims/python
#SBATCH -D /home/moritz/repos/MiComPy/
#SBATCH -J OD1s
#SBATCH -o /home/moritz/people/moritz/CDs/run3.out
#SBATCH -e /home/moritz/people/moritz/CDs/run3.err
#SBATCH -A b2014036
#SBATCH -t 3-00:00:00
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
from micompy.gene_clusterings.orthomcl.clustering import Clustering
from micompy.gene_clusterings.pfam_clusters.clustering import PfamClustering
from itertools import groupby
from pylab import *
from micompy.common.utils.iotl_annotations import *

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

name_map = {g.name  : g.conv_name() for g in all_genomes if g.name  !=  g.conv_name() }
rev_name_map = {v:k for k,v in name_map.iteritems()}
name_map.update({g.conv_name : g.conv_name for g in all_genomes })
rev_name_map.update({g.name : g.name for g in all_genomes })

mcl = orthoMCL(pjoin(analyses_root, "orthoMCL/"), all_genomes, "big_clustering")
clusters = Clustering(proteoms, pjoin(analyses_root, "clustering/"),"candidate_divs", mcl, checkm = pjoin(analyses_root,"checkm"),  name_map = name_map, rev_name_map = rev_name_map)
taxo_map  = {g.name :  (str(g.metadata['taxonomy_external'] +  "--" ) if g.metadata['taxonomy_external'] == g.metadata['taxonomy_external'] else "") + (str(g.metadata['phylum'] + "--") if g.metadata['phylum'] == g.metadata['phylum'] else "") + g.name for g in all_genomes}
g2clusters = {}
with open(pjoin(root, "big_gff.gff")) as handle:
    temp_dat = [l[:-1].split("\t")[8] for l in handle.readlines()]
big_gff_dict = {l.split(";")[0][3:] : {v.split("=")[0] : v.split("=")[1] for v in l.split(";")[1:]} for l in temp_dat}

with open("marker_hmms.txt") as handle:
    temp_markers = [l[:-1].split("\t") for l in handle.readlines()]
marker_hmms = { p : [] for p  in list(set([l[0] for l in temp_markers]))}
for l in tqdm(temp_markers):
    marker_hmms[l[0]] += [l[1] ]
del marker_hmms['pathway']
for v in big_gff_dict.itervalues():
    t=v['inference'].split(",")
    tt=t[0].split(":")
    v['inference_type'] = tt[0]
    v['inference_tool'] = tt[1]
    v['inference_quality'] = tt[2]
    tt=t[1].split(":")
    v['annotation_type'] = tt[0]
    v['annotation_tool'] = tt[1]
    v['annotation_id'] = tt[2]

with open("micompy/common/utils/pfam2go.txt") as handle:
    temp_pfam2go = [l[:-1].split() for l in handle.readlines()if l[0] != "!"]

pfam2go = { p.split(":")[1] : [] for p  in list(set([l[0] for l in temp_pfam2go]))}
for l in tqdm(temp_pfam2go):
    pfam2go[l[0].split(":")[1]] += [" ".join(l[3:]) ]
    
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

#pfams = PfamClustering(proteoms, pjoin(analyses_root, "pfam_clustering/"), "pfam_cdivs", checkm=pjoin(analyses_root,"checkm"), rev_name_map = rev_name_map)
#pfams.hmm_dict = pfams.parse_hmmer_results(pjoin(pfams.path,"hmmer_good_prots.raw"))
#pfams.post_process()

#tigrfams = PfamClustering(proteoms, pjoin(analyses_root, "tigrfam_clustering/"), "tigrfam_cdivs", checkm=pjoin(analyses_root,"checkm"), rev_name_map = rev_name_map)

#tigrfams.hmm_dict = tigrfams.parse_hmmer_results(pjoin(tigrfams.path,"tigrfam_all_prots.raw"))
#tigrfams.post_process()

for c in clusters:
    c.genes = set(c.genes)
for c in clusters:
    c.genomes = set(c.genomes)

#for c in pfams:
#    c.genes = set(c.genes)
#for c in pfams:
#    c.genomes = set(c.genomes)
        
#if __name__ == '__main__':
    
#    class Unbuffered(object):
#        def __init__(self, stream):
#            self.stream = stream
#        def write(self, data):
#            self.stream.write(data)
#            self.stream.flush()
#        def __getattr__(self, attr):
#            return getattr(self.stream, attr)

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
#    all_trees = []
#    ccs = []
#    threash = len([g for g in all_genomes if g.is_good()])*0.05
#    short_names = {g.name : g.metadata['short_name'] for g in all_genomes if g.is_good()}
#    for c in tqdm(clusters.single_copy_clusters()):
#        if len(c.genomes) > 3 and c.coreness < 0 :
#            ccs += [c.name]
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
#concat_sc_pfam_tree(pfams, pjoin(analyses_root, "tree_trunk","pfam_trunk"))

#
#cluster_overlaps = {(c1.name, c2.name) : float(len(c2.genomes.intersection(c1.genomes)))/min(len(c2.genomes), len(c1.genomes))  for  c2 in tqdm(clusters) for c1 in clusters if c1!=c2 and len(c1.genes) >20 and len(c2.genes)>20}
#overlaps = { v : {} for v in set([ k[0] for k in cluster_overlaps.keys() ])}
#for k,v in tqdm(cluster_overlaps.iteritems()):
#    overlaps[k[0]][k[1]] = v
#cluster_simi = DataFrame.from_dict(overlaps)
#cluster_simi = cluster_simi.fillna(0)

#cluster_overlaps_scg = {(c1.name, c2.name) : float(len(c2.genomes.intersection(c1.genomes)))/max(len(c2.genomes), len(c1.genomes))  for  c2 in tqdm(clusters.single_copy_clusters()) for c1 in clusters.single_copy_clusters() if c1!=c2# and len(c1.genes) >20 and len(c2.genes)>20}
#overlaps = { v : {} for v in set([ k[0] for k in cluster_overlaps_scg.keys() ])}
#for k,v in tqdm(cluster_overlaps_scg.iteritems()):
#    overlaps[k[0]][k[1]] = v
#cluster_simi_scg = DataFrame.from_dict(overlaps)
#cluster_simi_scg = cluster_simi_scg.fillna(0)

def print_metadata():
    SAG_data(all_genomes,pjoin(analyses_root, "tree_trunk","pfam_trunk","sag_metadat.txt"))
    
    env_data(all_genomes,pjoin(analyses_root, "tree_trunk","pfam_trunk","env_metadat.txt"))

"""
library(pheatmap)
library(RColorBrewer) 
dd = read.csv("~/repos/MiComPy/composition_matrix.csv",h=T, row.names=1)
mat = apply(dd > 0, 2, function(x) apply(dd > 0 , 2, function(y)  sum((x + y) == 2)/min(sum(x),sum(y))))

ddis = dist(mat)
tree = hclust(ddis)
clss = data.frame(class=cutree(tree,k=32))
clss$class = as.factor(clss$class)
write.table(clss, file="big_clusters.csv", quote=FALSE, col.names=FALSE, sep=",")
 pheatmap(mat, cluster_rows=tree, cluster_cols=tree, annotation_row=clss, show_rownames=FALSE, show_colnames=FALSE, legend=FALSE,annotation_legend=FALSE, border_color=NA, color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdGy")))(100), breaks = seq(-1,1,0.02))


pheatmap(mat, cluster_rows=tree, cluster_cols=tree, annotation_row=clss)
"""
from __builtin__ import sum

composition_matrix = {c.name : {g: len(v) for g,v in c.genome_2_gene_map.iteritems() } for c in sum([[c for c in l] for l in [clusters]],[])  if len(c.genomes) > 20}
composition_matrix = DataFrame.from_dict(composition_matrix)
composition_matrix = composition_matrix.fillna(0)
composition_matrix.to_csv("composition_matrix.csv")


big_clusts = DataFrame.from_csv("big_clusters.csv",header=None)
big_clusts={"cluster_" + str(cls) : list(big_clusts.loc[big_clusts[1] == cls].index) for cls in set(list(big_clusts[1]))}
big_clusts_compos = {k : (composition_matrix[v]>0).sum(axis=1)/len(v) for k,v in big_clusts.iteritems()}



float_to_rgb = lambda (r,g,b) : '#%02x%02x%02x' % (int(r*255), int(g*255), int(b*255))
cm = get_cmap('Paired')
color_map = {k : float_to_rgb(cm(1.*i/len(big_clusts))[0:3]) for i,k in enumerate(big_clusts_compos)}
denoised = {k: l[l > 0.2] for k,l in big_clusts_compos.iteritems()}

for k, l in denoised.iteritems(): 
    clust_data(k + "",l.to_dict(),pjoin(analyses_root, "tree_trunk","pfam_trunk",k + "_dat.txt"), color_map[k])




#l1= ['PF00687']
#l25 = ['PF01386']
#tt = set(sum(marker_hmms.values(),[])+ l1 +l25)
#hmms_to_expand = [ccc for ccc in pfams if ccc.name in tt] + [ccc for ccc in tigrfams if ccc.name in tt]
#ext_hmms = {c.name : [cc.name for cc in clusters if any([g in set(c.genes) for g in cc.genes]) ]   for c in tqdm(hmms_to_expand)   }
#ext_hmms['dummy'] = []

#marker_compos_ext = {k : (composition_matrix[sum([ext_hmms[c] for c in v + ['dummy']],[])]>0).sum(axis=1)/len(v) for k,v in temp_markers.iteritems() if len(v) > 0}

#temp_markers = {k :  [v for v in l if v in composition_matrix.columns] for k,l in marker_hmms.iteritems()}
#marker_compos = {k : (composition_matrix[v]>0).sum(axis=1)/len(v) for k,v in temp_markers.iteritems() if len(v) > 0}


#for k, l in marker_compos_ext.iteritems(): 
#    clust_data(k + "_ext",l.to_dict(),pjoin(analyses_root, "tree_trunk","pfam_trunk",k + "_ext_dat.txt"), color_map[k])

for k,l in big_clusts.iteritems():
    with open(pjoin(analyses_root, "tree_trunk","pfam_trunk",k + "_genes.txt"),"w") as handle:
        handle.writelines( k + "\n" + "\n".join(list(set([ big_gff_dict[g]['annotation_id'].split("_")[0]  for g in sum([list(c.genes) for c in clusters if c.name in  set(l)], []) if big_gff_dict.has_key(g) and big_gff_dict[g]['annotation_tool'] == "UniProtKB" ]))))

core = [[cc for cc in clusters if c == cc.name][0] for c in big_clusts['cluster_1'] ]
core = [c for c in core if float(len(c.genes))/len(c.genomes) < 1.05  ]
