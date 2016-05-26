import dendropy

trunk = dendropy.Tree.get_from_path(pjoin(analyses_root, "tree_trunk",  "full_alignment.tree") , "newick" )
forest = dendropy.TreeList.get_from_path(pjoin(clusters.path , "tree_set.tree"), "newick" )

distance = lambda x,y : dendropy.calculate.treecompare.symmetric_difference(x,y )

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
