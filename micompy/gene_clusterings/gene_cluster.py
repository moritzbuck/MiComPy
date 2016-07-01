from collections import Counter
from numpy import power
from numpy import log
from numpy import nan_to_num, prod



class GeneCluster(object):

    def __repr__(self): return '<%s object %s, annotated as  %s with %i genes from %i genomes>' % (self.__class__.__name__, self.name, self.annotation, len(self.genes), len(self.genomes))

    def __init__(self, clustering , genes, name = None, annotation = None ):

        self.clustering = clustering
        self.coreness = None
        if type(genes) == dict :
#            self.from_dict(genes)
            print "Don't forget to repair"
        else:
            self.name = name
            self.from_list(genes, annotation)
             
    def to_dict(self):
        return {u'name': self.name, u'annot_fraction': self.annot_fraction, u'annotation': self.annotation,  u'genes': self.genome_2_gene_map,  u'mapping': self.mapping, "coreness" : self.coreness}

    def from_dict(self,imp):
        self.name = imp['name']
        self.genes = imp['mapping'].keys()
        self.genomes = imp['genes'].keys()
        self.genome_2_gene_map =  imp['genes']
        
        self.annotation = imp['annotation']
        self.annot_fraction = imp['annot_fraction']
        self.mapping = imp['mapping']
        self.coreness = imp['coreness'] if imp.has_key('coreness') else None
        
    def from_list(self, genes, annotation):
        self.genomes = list(set([g.split("|")[0] for g in genes]))
        self.genes = [g.split("|")[1] for g in genes]
        self.genome_2_gene_map = {go : [ge.split("|")[1] for ge in genes if go == ge.split("|")[0]] for go in self.genomes}
       
        if annotation:
            self.annotation = annotation
            self.annot_fraction = None
            self.mapping = None

        else :
            sub_dict = {g : self.clustering.id2name_map[g] for g in self.genes}
            name_counts = Counter(sub_dict.values())
            total = sum([name_counts[z] for z in name_counts])
            annot_frac = float(name_counts.most_common()[0][1])/float(total)
            self.annotation = name_counts.most_common(1)[0][0]
            self.annot_fraction = annot_frac
            self.mapping = sub_dict
#        self.coreness = self.compute_coreness()

    def to_sequences(self, short=False, genome_name = False, subset = None):
        if not subset:
            subset = self.genomes
        seqs = []
        for g in self.genomes:
            if g in subset:
                genome = [f for f in self.clustering.genomes if f.split("/")[-1].split(".")[0] == g ][0]
                with open(genome, "r") as handle:
                    t_seqs = [s for s in SeqIO.parse(handle, "fasta") if s.id in self.genes]
                    if genome_name:
                        for s in t_seqs:
                            s.id = g
                            s.name = g
                    if short:
                        for s in t_seqs:
                            s.description = ""
                    seqs += t_seqs
        return seqs

    def calc_checksum(self, s):
        return str(sum(ord(c)*i for i,c in enumerate(s)))

    def align(self,output_file, block = True, genome_names = True, subset = None):
        if not subset:
            subset = self.genomes
        with open("temp.faa","w") as unalign:
            temp_seqs = self.to_sequences(short=True, genome_name = genome_names, subset = subset)
            SeqIO.write(temp_seqs, unalign, "fasta")
        sh.muscle("-in", "temp.faa","-out", "temp_aligned.faa")
        os.remove("temp.faa")
        if block:
            try:
                    sh.Gblocks("temp_aligned.faa", "-t=p", "-b5=h", "-b4=2", "-b2=0", "-b3=2000", "-b1=0.3")
            except:
                pass

            if os.path.exists("temp_aligned.faa-gb.htm"):
                shutil.move("temp_aligned.faa-gb", output_file)
                os.remove("temp_aligned.faa-gb.htm")
                sh.seqret("-sequence", "temp_aligned.faa", "-outseq", "nexus:" + ".".join(output_file.split(".")[:-1]) + ".nex")
                os.remove("temp_aligned.faa")
                return 1
            else :
                return 0
        else:
            sh.seqret("-sequence", "temp_aligned.faa", "-outseq", "nexus:" + ".".join(output_file.split(".")[:-1]) + ".nex")  
            shutil.move("temp_aligned.faa", output_file)
            return 1

    def tree_construction(self,alignment, outputtree):
        sh.FastTree("-out", outputtree,  alignment)

        
    def core_probability(self):
        present = prod([self.clustering.completnesses[self.clustering.rev_name_map[g]] for g in self.genomes])
        abscent = prod([1-v for k,v in self.clustering.completnesses.iteritems() if k not in self.genomes and k in self.clustering.genome2len.keys()])
        return present*abscent

    def non_core_probability(self):
        prob_of_random_pres = lambda g: 1.0 - power(float(len(self.clustering)-1)/len(self.clustering) , self.clustering.genome2len[self.clustering.rev_name_map[g]])
        prob_of_random_absc = lambda g: power(float(len(self.clustering)-1)/len(self.clustering) , self.clustering.genome2len[self.clustering.rev_name_map[g]])
        present = prod([prob_of_random_pres(g) for g in self.genomes])
        abscent = prod([prob_of_random_absc(k) for k in self.clustering.genome2len.keys() if k not in self.genomes])
        return present*abscent

    def non_core_probability_plural(self):
        prob_of_random_pres = lambda g: 1.0 - power(float(sum(self.clustering.genome2len.values())-len(self.genes))/sum(self.clustering.genome2len.values()) , self.clustering.genome2len[g])
        prob_of_random_absc = lambda g: power(float(sum(self.clustering.genome2len.values())-len(self.genes))/sum(self.clustering.genome2len.values()) , self.clustering.genome2len[g])
        present = prod([prob_of_random_pres(g) for g in self.genomes])
        abscent = prod([prob_of_random_absc(k) for k in self.clustering.genome2len.keys() if k not in self.genomes])
        return present*abscent

    def compute_coreness(self) :
        return log(self.core_probability()/self.non_core_probability())
