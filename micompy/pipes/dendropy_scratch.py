import dendropy
import itertools
from fractions import Fraction


forest = dendropy.TreeList.get_from_path(pjoin(clusters.path , "tree_set.tree"), "newick" )
trunk = dendropy.Tree.get_from_path(pjoin(analyses_root, "tree_trunk",  "full_alignment.tree") , "newick", taxon_namespace = forest[0].taxon_namespace )

distance = lambda x,y : dendropy.calculate.treecompare.robinson_foulds_distance(x,y)
tree_intersect = lambda x,y : set([l.taxon.label for l in x.leaf_iter()]).intersection(set([l.taxon.label for l in y.leaf_iter()]))
edit_dist_bitstrings = lambda x,y: sum([a != b for a,b in zip(list(x),list(y))])

overlaps = {i : {j : 2*float(len(tree_intersect(t1,t2)))/(len(t1.leaf_nodes())+len(t2.leaf_nodes())) for j,t2 in enumerate(forest)}  for i,t1 in tqdm(enumerate(forest))}

test = {t1 : {t2 : distance(t1,t2) for t2 in forest}  for t1 in tqdm(forest)}


def test(min_new_completness, min_count = None):
    if min_count:
       core = [c for c in clusters.single_copy_clusters() if len(c.genomes) > min_count]
    else :
       core = [c for c in clusters.single_copy_clusters() if c.coreness > min_coreness]
    temp = sum([c.genomes for c in core],[])
    gInCore = set(temp)
    main_trunk = [g for g in gInCore if float(temp.count(g))/len(core) > min_new_completness]
    use_core = [ c for c in core if len([g for g in main_trunk if g in c.genomes]) > len(main_trunk)/2]
    return len(use_core), len(main_trunk), len(use_core) * len(main_trunk)

[(x,test(float(x)/100)) for x in [t for t in xrange(0,50,1)]]

def test2(min_count):
    core = [c for c in clusters.single_copy_clusters() if len(c.genomes) > min_count]
    temp = sum([c.genomes for c in core],[])
    gInCore = set(temp)
    main_trunk = [g for g in gInCore if float(temp.count(g))/len(core) > min_new_completness]
    use_core = [ c for c in core if len([g for g in main_trunk if g in c.genomes]) > len(main_trunk)/2]
    return len(use_core), len(main_trunk), len(use_core) * len(main_trunk)


def closest_partition(tree, bipart):
    bitlist = [b == "1" for b in list(bipart.leafset_as_bitstring()[::-1])]
    tree_taxas = set([l.taxon for l in tree.leaf_iter()])
    trunk_taxas = set([l.taxon for l in trunk.leaf_iter()])
    common_taxas = trunk_taxas.intersection(tree_taxas)
    if len(common_taxas) < 4:
        return (0,0,0)
    subbitlist = [b if taxa in common_taxas else False for b, taxa in zip(bitlist, trunk.taxon_namespace) ]
    if (2*sum(subbitlist)) > len(subbitlist):
        subbitstring="".join(['0' if b else '1' for b in subbitlist])[::-1]
    else:
        subbitstring="".join(['1' if b else '0' for b in subbitlist])[::-1]

    subtree = tree.clone()
    subtree.prune_taxa(tree_taxas.difference(common_taxas))
    min_diff = min([edit_dist_bitstrings(f.leafset_as_bitstring(), subbitstring) for f in subtree.encode_bipartitions()])
    return float(min_diff),len(common_taxas), subbitstring.count("1")


dists = [-1]*len(trunk.encode_bipartitions())
for i,bipart in enumerate(trunk.encode_bipartitions()):
    if bipart.leafset_as_bitstring().count("1") > 1:
        print "Running bipart" , i
        tt = [closest_partition(tree,bipart) for tree in tqdm(forest) ]
        dists[i] [t[0]/t[2] for t in tt if t[2] > 1])
        

def supporting_partitions(tree, tree_trunk, shuffle = False):
    tree_taxas = set([l.taxon for l in tree.leaf_iter()])
    trunk_taxas = set([l.taxon for l in tree_trunk.leaf_iter()])
    common_taxas = trunk_taxas.intersection(tree_taxas)
    if len(common_taxas) < 4:
        return [(-1,-1)]*len(trunk.bipartition_encoding)

    subtree = tree.clone()
    if shuffle:
        subtree.shuffle_taxa()
    subtree.prune_taxa(tree_taxas.difference(common_taxas))
    
    out = []
    for i,bipart in enumerate(tree_trunk.encode_bipartitions()):
        bitlist = [b == "1" for b in list(bipart.leafset_as_bitstring()[::-1])]
        subbitlist = [b if taxa in common_taxas else False for b, taxa in zip(bitlist, tree_trunk.taxon_namespace) ]
        if (2*sum(subbitlist)) > len(subbitlist):
            subbitstring="".join(['0' if b else '1' for b in subbitlist])[::-1]
        else:
            subbitstring="".join(['1' if b else '0' for b in subbitlist])[::-1]
        if subbitstring.count("1") < 2 :
            out += [(-1, 0)]
        else :
            out += [(min([edit_dist_bitstrings(f.leafset_as_bitstring(), subbitstring) for f in subtree.encode_bipartitions()]), subbitstring.count("1"))]
    return out


support = [supporting_partitions(tree, trunk) for tree in tqdm(list(forest))]
vec = [[] for i in xrange(len(support[0]))]
counts = [0]*len(support[0])

#for l in support:
#    for i,v in enumerate(l):
#        if v != None :
#            counts[i] += 1
#            vec[i] += [ v ]

#scores = [ float(a)/b if b != 0 else None for a,b in zip(vec,counts)]
tt = {c : {"edge_" + str(i) : v for i,v in  enumerate(l)} for l,c in zip(support, ccs)}
support_mat = DataFrame.from_dict(tt)
support_mat=support_mat.loc[support_mat.apply(lambda l: not all([v[1]<1 for v in l]), axis=1)]
support_mat=support_mat[support_mat.columns[support_mat.apply(lambda l: not all([v[1]<1 for v in l]), axis=0)]]

bstraps = []
for z in tqdm(xrange(100)):
    shuffled = trunk.clone()
    sup = [supporting_partitions(tree, shuffled, shuffle = True) for tree in list(forest)]
    perfects = [0]*len(support[0])
    cou = [0]*len(support[0])
    for l in sup:
        for i,v in enumerate(l):
            if v != None:
                perfects[i] += v
                cou[i] += 1
    bstraps += [ [ float(a)/b if b != 0 else None for a,b in zip(perfects,cou)] ]

support3 = [sum([s < bs[i] for bs in bstraps]) for i,s in enumerate(scores)]
for i,bipart in enumerate([b for b in trunk.encode_bipartitions() if b.leafset_as_bitstring().count("1") > 1]):
            for e in trunk.edges():
                    if e.bipartition == bipart:
                            e.label = str(scores[i]/bipart.leafset_as_bitstring().count("1"))

def nCk(n,k):
    mul=lambda x,y:x*y 
    return int( reduce(mul, (Fraction(n-i, i+1) for i in range(k)), 1) )

def supporting_quartets(tree, tree_trunk, shuffle = False, i=None):
    tree_taxas = set([l.taxon for l in tree.leaf_iter()])
    trunk_taxas = set([l.taxon for l in tree_trunk.leaf_iter()])
    common_taxas = trunk_taxas.intersection(tree_taxas)
    if len(common_taxas) < 4:
        return "To few common taxa"

    if shuffle:
        subtree.shuffle_taxa()
    print "making tree", i, "preparing", nCk(len(common_taxas), 4), "quartets from ", len(common_taxas), "taxa"
    quartets = [s for s in tqdm(itertools.combinations(common_taxas, 4))]
    counts = 0
    print "and evaluating them"
    for q in tqdm(quartets):
        subtree = tree.clone()
        subtrunk = tree_trunk.clone()
        subtree.prune_taxa(tree_taxas.difference(q))
        subtrunk.prune_taxa(trunk_taxas.difference(q))
        
        strip_tree(subtrunk)
        strip_tree(subtree)
        subtrunk.reroot_at_edge([e for e in subtrunk.edges() if e.head_node.taxon == q[0]][0])
        subtree.reroot_at_edge([e for e in subtree.edges() if e.head_node.taxon == q[0]][0])
        if str(subtree) == str(subtrunk):
            counts += 1

    return (counts, len(quartets), len(common_taxas))

def strip_tree(tree):
    for e in tree.edges():
        e.length = None
    for n in tree.internal_nodes():
        n.label = None
