
import os
import sh
import json
from Bio import SeqIO
import Bio
from collections import Counter
from tqdm import tqdm
import sys
import shutil
from pandas import DataFrame
from numpy import nan_to_num, prod
from orthmcl_tools.orthoMCL import orthoMCL;
import dendropy
import re
import operator
import pandas
from numpy import power
from numpy import log
import shutil
    
class Clustering(object):

    def __repr__(self): return '<%s object "%s with %i clusters">' % (self.__class__.__name__, self.name, len(self))
    def __iter__(self): return iter(self.clusters)
    def __len__(self): return len(self.clusters)
    def __getitem__(self, key): return self.clusters[key]

        
    def __init__(self,proteoms,  out_path, name, mcl, gff = None, seq_type="proteins", checkm = None, name_map = None):

        self.genomes = proteoms
        self.seq_type = seq_type
        self.path = out_path
        self.oMCL_path = mcl.out_dir
        if not os.path.exists(out_path):
            os.makedirs(out_path)

        if checkm:
            self.checkm = pandas.read_table(checkm, skiprows = 3, index_col = 0, sep = r"\s*", names = ["genome", "lineage","nb_gen", "nb_markers", "nb_sets", "0", "1","2","3","4","5+", "completness", "contamination", "heterogeneity"], comment = '-')
            self.completnesses = (self.checkm['completness']/100.0).to_dict()
#            self.completnesses = {unicode(k): v if v < 0.95 else 0.95 for k,v in self.completnesses.iteritems()}
        self.seed = 23
        self.name = name
        self.base = self.path + self.name
        self.scg_tree =  self.base +"_nodes_labeled.tree"
        self.name_map = name_map
        
        self.processed_clusters = self.base + ".json"
        self.align_path = self.base + "_align/"
        self.scc_align_path = self.base + "_scc_align/"
        self.db = self.oMCL_path + "goodProteins.fasta"
        if os.path.exists(self.processed_clusters):
            with open(self.processed_clusters) as file:
                self.clusters= [GeneCluster(self,c) for c in json.load(file)]
        else :
            self.clusters = []

        if gff:
            with open(gff,"r") as handle:
                temp = [ {v.split("=")[0] : v.split("=")[1] for v in l.split("\t")[-1].split(";")} for l in handle.readlines() ]
                self.id2name_map = {cds['locus_tag'] : cds['product'] for cds in temp if cds.has_key('product')}
        else :
            self.id2name_map = {}
            self.gene2genome = {}
            self.genome2len = {}
            for g in self.genomes:
                with open(g) as handle:
                    temp = {s.id : " ".join(s.description.split()[1:]) for s in SeqIO.parse(handle,"fasta")}
                    self.genome2len[g.split("/")[-1].split(".")[0]] = len(temp.keys())
                    self.gene2genome.update({gene : ".".join(g.split(".")[:-1]).split("/")[-1] for gene in temp.keys()})
                    self.id2name_map.update(temp)

                
            
        self.orthoMCL = mcl #orthoMCL(self.oMCL_path, self.genomes, self.name)
        self.raw_clusters = self.orthoMCL.out_mcl
        self.anc_rec_path = out_path + "AncRec/"
    
    def single_copy_clusters(self):
        return [c for c in self.clusters if len(c.genes) == len(c.genomes) and len(c.genomes) > 1]

    def almost_single_copy_clusters(self):
        return [c for c in self.clusters if len(c.genes) == (len(self.assemblies)-1) and (len(c.genomes)==len(self.assemblies)-1)]
    
    def c_name2name(self, line):
        oline = line
        for name , conv_name in self.name_map.iteritems():
            oline = oline.replace(conv_name, name)
        return oline
    
    def post_process(self):
        print "Post processing cluster:"

        with open(self.raw_clusters) as c_file:
            if self.name_map:
                self.clusters=[GeneCluster(self, name=l[:-1].split(": ")[0], genes = self.c_name2name(l[:-1].split(": ")[1]).split()) for l in tqdm(c_file.readlines())]
            else :
                self.clusters=[GeneCluster(self, name=l[:-1].split(": ")[0], genes = l[:-1].split(": ")[1].split()) for l in tqdm(c_file.readlines())]
            
        print "Post processing single genes:"
        non_singletons = set(sum([c.genes for c in self.clusters],[]))
        for i in tqdm(self.id2name_map.keys()):
            if i not in non_singletons:
                self.clusters += [GeneCluster(self, name = i,  genes =  [self.gene2genome[i] + "|" + i ])]
            
        with open(self.processed_clusters, 'w') as outfile:
            json.dump([c.to_dict() for c in self.clusters], outfile,  indent=4, sort_keys=True)

    def run(self):
#        self.orthoMCL.full_pipe()
        print "post-process"
        self.post_process()
        
    def align_all(self):
        print "aligning the hell out of it"
        if os.path.exists(self.align_path):
             shutil.rmtree(self.align_path)
        os.makedirs(self.align_path)

        for i,clust in tqdm(enumerate(self.single_copy_clusters())):
            clust.align(self.align_path + str(i) + ".fasta")

        
    def cat_align(self):
        print "CONCATENATE!!!11"
        cat_seqs = None
        for f in os.listdir(self.align_path):
            if "fasta" in f:
                with open(self.align_path + f,"r") as file:
                    seqs = [s for s in SeqIO.parse(file,"fasta")]
                    if not cat_seqs:
                        cat_seqs = seqs
                        order = [s.id.split("|")[0] for s in cat_seqs]
                        for i,s in enumerate(cat_seqs):
                            s.id = order[i]
                    else :
                        seqs_sorted = [[z for z in seqs if o in z.id][0] for o in order]
                        cat_seqs = [ s + seqs_sorted[i] for i,s in enumerate(cat_seqs) ]

        for i,s in enumerate(cat_seqs):
            s.id = order[i]
            s.description = "composite of " + str(len([f for f in os.listdir(self.align_path) if "fasta" in f])) + " single copy genes-clusters"
        with open(self.base + "_cat_align.fasta","w") as outp:
            SeqIO.write(cat_seqs,outp,"fasta")
        return cat_seqs
        
    def tree_construction(self,root = None, sccs = False):
        threads = 16 
        print "build a tree"
        if os.path.exists(self.base + "RAxML/" ):
            sh.rm("-r", self.base + "RAxML/")
        os.makedirs(self.base + "RAxML/")

        if self.seq_type == "proteins" :
            model = "PROTGAMMALG"
        else:
            model = "GTRGAMMA"

        alignment = self.base + "_scc_cat_align.fasta" if sccs else self.base + "_cat_align.fasta"
        
        sh.raxmlHPC_PTHREADS_AVX("-w", self.base + "RAxML/", "-T", threads-2, "-m", model, "-p", self.seed, "-#", 20, "-s", alignment, "-n", "T13", "-o", root) 
        print "boostrap dat tree"
        sh.raxmlHPC_PTHREADS_AVX("-w", self.base + "RAxML/", "-T", threads-2, "-m", model, "-p", self.seed, "-b", self.seed, "-#", 100, "-s", alignment, "-n", "T14", "-o", root)
        print "combine"
        sh.raxmlHPC_AVX("-m", "GTRCAT", "-w", self.base + "RAxML/", "-p", self.seed, "-f", "b", "-t", self.base + "RAxML/"+"RAxML_bestTree.T13", "-z",self.base + "RAxML/"+ "RAxML_bootstrap.T14", "-n", "T15", "-o", root)
        print "clean up"
        if os.path.exists(self.base + "_branches_labeled.tree"):
            os.remove(self.base + "_branches_labeled.tree")
            os.remove(self.base + "_nodes_labeled.tree")
        sh.ln("-s",  self.base + "RAxML/RAxML_bipartitionsBranchLabels.T15", self.base +"_branches_labeled.tree")
        sh.ln("-s",  self.base + "RAxML/RAxML_bipartitions.T15", self.scg_tree)

        
    def rm_genome(self, name):
        self.assemblies = [a for a in self.assemblies if name not in a.name]
        for c in self:                                                             
            c.genomes = [g for g in  c.genomes if name not in g]
        for c in self:                                                             
            c.genes = [g for g in  c.genes if name not in g]


    def keep_genomes(self, genomes):
        self.assemblies = [a for a in self.assemblies if a.name in genomes]
        for c in self:                                                             
            c.genomes = [g for g in  c.genomes if g in genomes ]
        for c in self:                                                             
            c.genes = [g for g in  c.genes if len([o for o in genomes if g.count(o) ==1]) > 0 ]


    def cooccurence_matrix(self):
        matrix = DataFrame(data=0, index = [a.name for a in self.assemblies], columns=[a.name for a in self.assemblies])
        for a in self.assemblies:
            for b in self.assemblies:
                count = 0
                for c in self:
                    if a.name in c.genomes and b.name in c.genomes:
                        count += 1
                matrix[a.name][b.name] = count
                
        return matrix

    def make_cluster_bmft(self):
        cluster_table = DataFrame.from_dict({i : {k : len(v) for k,v in c.to_dict()['genes'].iteritems()} for i,c in enumerate(self)}, orient='index')
        cluster_table = cluster_table.apply(nan_to_num)
        cluster_table['annotations'] = [c.annotation for c in self]
        cluster_table['qual_annot'] = [c.annot_fraction for c in self]
        cluster_table['genes'] = [";".join(c.genes) for c in self]
        return cluster_table
        

    def ancestral_reconstr(self, outgroup = None):
        print "Make trait File"
        if not os.path.exists(self.anc_rec_path):
            os.makedirs(self.anc_rec_path)
        json2nexus(self.processed_clusters, self.anc_rec_path + "Trait_genome.fasta", self.anc_rec_path + "Traits.nex"  )
        with open(self.anc_rec_path + "Trait_genome.fasta") as handle:
            seqs=[s for s in SeqIO.parse(handle,"fasta") ]
        for s in seqs:
            s.seq.alphabet = Bio.Alphabet.DNAAlphabet()
        with open(self.anc_rec_path + "Matrix.nex","w") as handle:
            SeqIO.write(seqs,handle,"nexus")
        os.system("sed -i 's/format datatype=dna missing=? gap=- interleave;/format interleave datatype=standard   gap=- symbols=\"0123456789ABCDEF\";/' " + self.anc_rec_path + "Matrix.nex")

        print "Make tree file"
        new = dendropy.Tree.get(path=self.scg_tree, schema="newick")
#        new.write(path = self.anc_rec_path + "Tree.nex", schema="nexus")
        tree_string = re.sub("[0-9]{2,}","", re.sub(r':[0123456789.]+','', str(new)))
        taxons = {x : i for i,x in enumerate([t.taxon.label for t in new.nodes() if t.is_leaf()])}
        for k in taxons.keys():
            tree_string=re.sub(k+"([),])",str(taxons[k]+1)+r"\1", tree_string)

        taxons_t = sorted(taxons.items(), key=operator.itemgetter(1))

        t_table = "\n".join([("\t\t" + str(p[1]+1) + "\t" + p[0].replace(" ","_") +",")  for p in taxons_t])
        paup_tree = """#NEXUS
Set AllowPunct=Yes;
BEGIN TREES;
\tTRANSLATE
%s\b
;
TREE Bioperl_1 = [&U] %s ;
END; 
""" % (t_table , tree_string)
            
        paup_script = """
#NEXUS
BEGIN PAUP;
      log file=%sPaup.log replace=yes start=yes;
      exe %sMatrix.nex;
      exe %sAssumptions.nex;
      exe %sTree.nex;
      ctype GeneCopy:All;
      pset Opt=AccTran;

      taxset out= %d;
      outgroup %s;
      set root=outgroup outroot=monophyl;
      root root=outgroup;      

      describetrees 1/Xout=both ApoList=Yes ChgList=yes; [brLens=yes;]
      savetrees from=1 to=1 format=nexus root=yes file=%sOutTree.rtree;
      q;
END;
""" % (self.anc_rec_path, self.anc_rec_path, self.anc_rec_path, self.anc_rec_path, len(seqs), outgroup if outgroup else "out",  self.anc_rec_path)

        assumptions="""
#NEXUS

BEGIN ASSUMPTIONS;
USERTYPE GeneCopy (STEPMATRIX) = 16
      0   1   2   3   4   5   6   7   8   9   A   B   C   D   E   F
[0]   .   10.0 11.0 11.2 11.4 11.6 11.8 12.0 12.2 12.4 12.6 12.8 13.0 13.2 13.4 13.6
[1]   5.0 .   1.0 1.2 1.4 1.6 1.8 2.0 2.2 2.4 2.6 2.8 3.0 3.2 3.4 3.6
[2]   5.2 0.2 .   0.2 0.4 0.6 0.8 1.0 1.2 1.4 1.6 1.8 2.0 2.2 2.4 2.6
[3]   5.4 0.4 0.2 .   0.2 0.4 0.6 0.8 1.0 1.2 1.4 1.6 1.8 2.0 2.2 2.4
[4]   5.6 0.6 0.4 0.2 .   0.2 0.4 0.6 0.8 1.0 1.2 1.4 1.6 1.8 2.0 2.2
[5]   5.8 0.8 0.6 0.4 0.2 .   0.2 0.4 0.6 0.8 1.0 1.2 1.4 1.6 1.8 2.0
[6]   6.0 1.0 0.8 0.6 0.4 0.2 .   0.2 0.4 0.6 0.8 1.0 1.2 1.4 1.6 1.8
[7]   6.2 1.2 1.0 0.8 0.6 0.4 0.2 .   0.2 0.4 0.6 0.8 1.0 1.2 1.4 1.6
[8]   6.4 1.4 1.2 1.0 0.8 0.6 0.4 0.2 .   0.2 0.4 0.6 0.8 1.0 1.2 1.4
[9]   6.6 1.6 1.4 1.2 1.0 0.8 0.6 0.4 0.2 .   0.2 0.4 0.6 0.8 1.0 1.2
[A]   6.8 1.8 1.6 1.4 1.2 1.0 0.8 0.6 0.4 0.2 .   0.2 0.4 0.6 0.8 1.0
[B]   7.0 2.0 1.8 1.6 1.4 1.2 1.0 0.8 0.6 0.4 0.2 .   0.2 0.4 0.6 0.8
[C]   7.2 2.2 2.0 1.8 1.6 1.4 1.2 1.0 0.8 0.6 0.4 0.2 .   0.2 0.4 0.6
[D]   7.4 2.4 2.2 2.0 1.8 1.6 1.4 1.2 1.0 0.8 0.6 0.4 0.2 .   0.2 0.4
[E]   7.6 2.6 2.4 2.2 2.0 1.8 1.6 1.4 1.2 1.0 0.8 0.6 0.4 0.2 .   0.2
[F]   7.8 2.8 2.6 2.4 2.2 2.0 1.8 1.6 1.4 1.2 1.0 0.8 0.6 0.4 0.2 .
;
charset All =1-.;
END;
"""
        print "Write nexus files"
        with open(self.anc_rec_path + "Script.nex", "w") as handle:
            handle.writelines(paup_script)

        with open(self.anc_rec_path + "Assumptions.nex", "w") as handle:
            handle.writelines(assumptions)

        with open(self.anc_rec_path + "Tree.nex", "w") as handle:
            handle.writelines(paup_tree)

        print "Run and parse paup"
        os.system("paup4a146_centos64 -f -n " + self.anc_rec_path + "Script.nex > " +  self.anc_rec_path + "OutPut.nex")
        sh.perl("/pica/h1/moritz/repos/Pyscratches/20150610_orthomcl_tools/orthmcl_tools/parsePaupLog.pl", self.anc_rec_path + "OutPut.nex")
        os.remove(self.anc_rec_path + "OutPut.nex.nodes.johan")
        print "Merge data with clusters"
        merged_json = mergeJSONwPAUP(self.processed_clusters, self.anc_rec_path + "OutPut.nex.changes")
        with open(self.anc_rec_path + "clusters.json","w") as handle:
            json.dump(merged_json, handle,  indent=4, sort_keys=True)

        return merged_json
    
        
    def sccs_cat_align(self):
        all_sccs = {}
        for a in self.assemblies:
            with open(a.sccs_genes) as fasta:
                all_sccs[a.name] = [s for s in SeqIO.parse(fasta,"fasta")]
        cogs = list(set(sum([[s.id for s in seqs] for seqs in all_sccs.values()],[])))
        cogs.sort()

        a2id = {a.name : self.assembly_id(a.name + a.name)  for a in self.assemblies}
        id2a = {a2id[a] : a for a in a2id}
        
        for c in cogs:
            c_seqs = []
            for a,seqs in all_sccs.iteritems():
                seq = [s for s in seqs if s.id == c]
                if len(seq) == 1:
                    seq = seq[0]
                else:
                    continue
                seq.id = a2id[a]
                seq.description = ""
                c_seqs.append(seq)
            if len(c_seqs) == len(self.assemblies):
                with open("temp.faa","w") as fasta:
                    SeqIO.write(c_seqs,fasta,"fasta")
                print "aligning genes for",c
                sh.muscle("-in", "temp.faa","-out", "temp_aligned.faa")
                os.remove("temp.faa")
                try:
                    sh.Gblocks("temp_aligned.faa", "-tD")
                except:
                    pass
        
                os.remove("temp_aligned.faa")
                if os.path.exists("temp_aligned.faa-gb.htm"):
                    os.remove("temp_aligned.faa-gb.htm")
                    if not os.path.exists(self.scc_align_path):
                        os.makedirs(self.scc_align_path)
                    shutil.move("temp_aligned.faa-gb", self.scc_align_path + c + ".ffn")


        print "CONCATENATE!!!11"
        cat_seqs = None
        for f in os.listdir(self.scc_align_path):
            if "ffn" in f:
                with open(self.scc_align_path + f,"r") as file:
                    seqs = [s for s in SeqIO.parse(file,"fasta")]
                    if not cat_seqs:
                        cat_seqs = seqs
                        order =  [s.id for s in cat_seqs]
                        for i,s in enumerate(cat_seqs):
                            s.id = order[i]
                    else :
                        seqs_sorted = [[z for z in seqs if o in z.id][0] for o in order]
                        cat_seqs = [ s + seqs_sorted[i] for i,s in enumerate(cat_seqs) ]
            for i,s in enumerate(cat_seqs):
                s.id = id2a[order[i]]
                s.description = "composite of " + str(len([f for f in os.listdir(self.scc_align_path) if "ffn" in f])) + " single copy cogs"
            with open(self.base + "_scc_cat_align.fasta","w") as outp:
                SeqIO.write(cat_seqs,outp,"fasta")
                
        return all_sccs



def trash(data):
    lines = []
    data_sub = [d for d in data if any([g in genomes_of_interest for g in d['genes'].keys()] ) ]
    head = ["ID", "hypothetical.name", "name.confidence"]+  genomes_of_interest+ [ "ancestor" , ancestry , "genomes.of.interest", "other.genomes", "genes"]
    for d in data_sub:
        name = d['name'] 
        hname = d['annotation']
        c_name = d['annot_fraction']
        if d.has_key('ancestral') and d['ancestral'].has_key("node_55_to_node_54"):
            ancestor = d['ancestral']["node_55_to_node_54"]
            frmo = ancestor['from']
            diff = ancestor['to']-ancestor['from']
        else :
            frmo = "NA"
            diff = 0
        if d.has_key('ancestral') and d['ancestral'].has_key("node_54_to_node_53"):
            ancestor = d['ancestral']["node_54_to_node_53"]
            frmo2 = ancestor['from']
            diff2 = ancestor['to']-ancestor['from']
        else :
            frmo2 = "NA"
            diff2 = 0
    
        counts = [len(d['genes'][g]) if d['genes'].has_key(g) else 0 for g in genomes_of_interest]
        icounts = len([g for g in genomes if g in genomes_of_interest and d['genes'].has_key(g)])
        ocounts = len([g for g in genomes if g not in genomes_of_interest and d['genes'].has_key(g)])
        genes = ";".join(d['mapping'].keys())
        line = [ name, hname, c_name ] + counts + [ frmo, diff, frmo2, diff2, icounts, ocounts, genes ] 
        lines +=  [[str(l) for l in line]]
