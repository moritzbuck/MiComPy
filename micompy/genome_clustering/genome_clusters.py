class MashCluster(object):
    """docstring for ."""
    def __repr__(self): return '<%s object %s, with %i genomes, of probable family %s>' % (self.__class__.__name__, self.name, len(self.genomes), self.consensus_family )



    def __init__(self, mash_matrix, genomes, taxo, name = None):
        self.genomes = genomes
        c = [g.name  for g in self.genomes]

        if type(mash_matrix) != type(None) :
            self.mash_matrix = mash_matrix[c].loc[c]
            tt = self.mash_matrix.sum()
            self.representative_genome = [g for g in self.genomes if g.name == min(zip(tt.index,tt.values),key = lambda x : x[1])[0]][0]

        if not name:
            self.name = self.representative_genome.name

        self.taxoDb = taxo
        self.consensus_family = self.get_taxo()
        self.consensus_family = max(self.consensus_family, key = lambda x: self.consensus_family[x])

    def get_taxo(self, level = "family"):
        all_taxa = [g.get_taxo(self.taxoDb).get(level) for g in self.genomes]
        set_taxa = set(all_taxa)
        return {t : all_taxa.count(t) for t in list(set_taxa)}

    def set_consensus_family(self):
        self.consensus_family = self.get_taxo()
        self.consensus_family = max(self.consensus_family, key = lambda x: self.consensus_family[x])

    def add_genome(self,g):
        self.genomes += [g]
        self.set_consensus_family()
