import dendropy

forest = dendropy.TreeList.get_from_path(pjoin(clusters.path , "tree_set.tree"), "newick" )
trunk = dendropy.Tree.get_from_path(pjoin(analyses_root, "tree_trunk",  "full_alignment.tree") , "newick", taxon_namespace = forest[0].taxon_namespace )

distance = lambda x,y : dendropy.calculate.treecompare.symmetric_difference(x,y)
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
        

def supporting_partitions(tree, tree_trunk):
    tree_taxas = set([l.taxon for l in tree.leaf_iter()])
    trunk_taxas = set([l.taxon for l in tree_trunk.leaf_iter()])
    common_taxas = trunk_taxas.intersection(tree_taxas)
    if len(common_taxas) < 4:
        return (0,0,0)

    subtree = tree.clone()
    subtree.prune_taxa(tree_taxas.difference(common_taxas))
    
    out = []
    for i,bipart in enumerate(tree_trunk.encode_bipartitions()):
        if bipart.leafset_as_bitstring().count("1") > 1:
            bitlist = [b == "1" for b in list(bipart.leafset_as_bitstring()[::-1])]
            subbitlist = [b if taxa in common_taxas else False for b, taxa in zip(bitlist, tree_trunk.taxon_namespace) ]
            if (2*sum(subbitlist)) > len(subbitlist):
                subbitstring="".join(['0' if b else '1' for b in subbitlist])[::-1]
            else:
                subbitstring="".join(['1' if b else '0' for b in subbitlist])[::-1]
            if subbitstring.count("1") < 2 :
                out += [None]
            else :
                out += [min([edit_dist_bitstrings(f.leafset_as_bitstring(), subbitstring) for f in subtree.encode_bipartitions()]) == 0]
    return out

bstraps = []
for z in tqdm(xrange(100)):
    shuffled = trunk.clone()
    shuffled.shuffle_taxa()
    support = [supporting_partitions(tree, shuffled) for tree in forest]
    perfects = [0]*len(support[0])
    counts = [0]*len(support[0])
    for l in support:
        for i,v in enumerate(l):
            if v:
                perfects[i] += 1
                if v != None :
                    counts[i] += 1
    bstraps += [ [ float(a)/b if b != 0 else None for a,b in zip(perfects,counts)] ]

