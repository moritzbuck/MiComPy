from micompy.common.utils.intrasimilarity import NIC_similarity
import os
from os.path import join as pjoin
import json
from pandas import read_table
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from subprocess import call
from micompy.gene_clusterings.pfam_clusters.clustering import sc_pfams
from micompy.gene_clusterings.pfam_clusters.clustering import pfam_equivs
from tqdm import tqdm
import sh

def annotation(genomes, cpus = 1, clean = False):
    for g in genomes:
        if clean or not g.is_annotated():
            g.prokka(cpus = cpus)

def concat_core_tree(clusters, path, min_clust_size = 10, min_coreness = 0, min_new_completness = 0.3 ):
    
    if not os.path.exists(path):
        os.makedirs(path)
        
    for c in tqdm(clusters.single_copy_clusters()):
        if len(c.genomes) > min_clust_size:
            c.coreness = c.compute_coreness()
            
    core = [c for c in clusters.single_copy_clusters() if c.coreness > min_coreness]
    temp = sum([c.genomes for c in core ],[])
    gInCore = set(temp)
    main_trunk = [g for g in gInCore if float(temp.count(g))/len(core) > min_new_completness]
    use_core = [ c for c in core if len([g for g in main_trunk if g in c.genomes]) > len(main_trunk)/2]
 
    align_path = pjoin(path, "alignments")
    
    if not os.path.exists(align_path):
        os.makedirs(align_path)
    
    for c in tqdm(use_core):
        c.align(pjoin(align_path, c.name + ".faa"), block=True, subset = main_trunk)

    concat_seqs = {}
    for g in main_trunk:
        concat_seqs[g] = ""

    for f in [pjoin(align_path, t) for t in os.listdir(align_path) if "blocked.faa" in t]:
        with open(f) as handle:
            entries = {s.id : str(s.seq) for s in SeqIO.parse(handle,"fasta")}
        length = len(entries.values()[0])
        assert all([len(s) == length for s in entries.values()])

        for g in concat_seqs:
            if entries.has_key(pfam_clusters.rev_name_map[g]):
                concat_seqs[g] += entries[pfam_clusters.rev_name_map[g]]
            else:
                concat_seqs[g] += "X"*length
    alignment = pjoin(path, "full_alignment.faa")
    with open(alignment,"w") as handle:
        SeqIO.write([SeqRecord(Seq(s), id=g) for g,s in concat_seqs.iteritems()], handle, "fasta")
            
    #    sh.FastTree("-out", pjoin(path, "full_alignment.tree"), pjoin(path, "full_alignment.faa") )

    print "build a tree"
    rax_path = pjoin(path, "RAxML/")
    if os.path.exists(rax_path ):
            sh.rm(rax_path)
    os.makedirs(rax_path)

    model = "PROTGAMMALGF"
    seed = 42 

    sh.raxmlHPC_AVX('-m', "PROTGAMMALGF", "-T", threads-2 , '-p', sede, '-s', alignment, '-n', 'tree', '-w', rax_path, '-f', 'a', '-x', 1, '-N', 'autoMR')

    
    taxo = {g.name : g.name +( "|" + g.metadata['taxonomy_external'] if g.metadata['taxonomy_external'] == g.metadata['taxonomy_external'] else "" ) + ("|" + g.metadata['phylum'] if g.metadata['phylum'] == g.metadata['phylum'] else "" )  for g in all_genomes if g.is_good()}
    renaming_tree(pjoin(path, "full_alignment.tree"),pjoin(path, "full_alignment.taxo.tree"), taxo)

def concat_sc_pfam_tree(pfam_clusters, path, min_new_completness = 0.3 ):
    
    if not os.path.exists(path):
        os.makedirs(path)

    sc_pfam_clusts = [c for c in pfam_clusters if c.name in sc_pfams]

    remove = []
    for ll in pfam_equivs:
        tt = [c for c in sc_pfam_clusts if c.name in ll]
        keep = max(tt, key=lambda x: len(x.genomes))
        remove += [t for t in tt if t != keep]

    remove = set(remove)
    sc_pfam_clusts = [c for c in sc_pfam_clusts if c not in remove]
        
    for c in sc_pfam_clusts :
        bads = {g : l for g,l in c.genome_2_gene_map.iteritems() if len(l) >1}
        c.black_list = []
        for g in bads.iterkeys():
            es = [v for v in zip(pfam_clusters.hmm_dict[c.name]['cdss'],pfam_clusters.hmm_dict[c.name]['es']) if v[0].split("|")[0]==g]
            es.sort(key=lambda x : x[1])
            c.black_list += [e[0] for e in es][1:]
            
    temp = sum([list(c.genomes) for c in sc_pfam_clusts ],[])
    gInCore = set(temp)
    main_trunk = [pfam_clusters.rev_name_map[g] for g in gInCore if float(temp.count(g))/len(sc_pfam_clusts) > min_new_completness]
    use_core = [ c for c in sc_pfam_clusts if len([g for g in main_trunk if g in c.genomes]) > len(main_trunk)/2]
 
    align_path = pjoin(path, "alignments")
    
    if not os.path.exists(align_path):
        os.makedirs(align_path)
    
    for c in tqdm(use_core):
        c.align(pjoin(align_path, c.name + ".faa"), block=True, subset = main_trunk)

    concat_seqs = {}
    for g in main_trunk:
        concat_seqs[pfam_clusters.rev_name_map[g]] = ""

    for f in [pjoin(align_path, t) for t in os.listdir(align_path) if ".faa" in t]:
        with open(f) as handle:
            entries = {s.id : str(s.seq) for s in SeqIO.parse(handle,"fasta")}
        length = len(entries.values()[0])
        assert all([len(s) == length for s in entries.values()])

        for g in concat_seqs:
            if entries.has_key(pfam_clusters.rev_name_map[g]):
                concat_seqs[g] += entries[pfam_clusters.rev_name_map[g]]
            else:
                concat_seqs[g] += "X"*length
    alignment = pjoin(path, "full_alignment.faa")

    concat_seqs = {k : v for k,v in concat_seqs.iteritems() if v.count("-") + v.count("X") != len(v)}
    with open(alignment,"w") as handle:
        SeqIO.write([SeqRecord(Seq(s), id=g) for g,s in concat_seqs.iteritems()], handle, "fasta")
            
    #    sh.FastTree("-out", pjoin(path, "full_alignment.tree"), pjoin(path, "full_alignment.faa") )

    print "build a tree"
    rax_path = pjoin(path, "RAxML/")
    if os.path.exists(rax_path ):
            sh.rm(rax_path)
    os.makedirs(rax_path)

    model = "PROTGAMMALG"
    seed = 42 

    os.system(" ".join(['raxmlHPC_AVX', '-m', 'PROTGAMMALGF', '-T', str(14) , '-p', str(23), '-s', alignment, '-n', 'tree', '-w', rax_path, '-f', 'a', '-x', str(1), '-N', 'autoMR']))

    
    taxo = {g.name : g.name +( "|" + g.metadata['taxonomy_external'] if g.metadata['taxonomy_external'] == g.metadata['taxonomy_external'] else "" ) + ("|" + g.metadata['phylum'] if g.metadata['phylum'] == g.metadata['phylum'] else "" )  for g in all_genomes if g.is_good()}
    renaming_tree(pjoin(path, "full_alignment.tree"),pjoin(path, "full_alignment.taxo.tree"), taxo)



    
def cluster_genomes(genomes,output, cutoff=0.95):
    for g in genomes:
        if not g.size:
            g.compute_size()
    genomes.sort(key = lambda x : x.size, reverse=True)

    left = set(genomes)
    clusters = {}
    for g in genomes:
        if g in left:
            #make blast db of g (larger genome)
            blast_db_files = [g.ref + ".nhr", g.ref + ".nin",  g.ref + ".nsq"]
            blast_db_cmd = ["makeblastdb" ,"-in", g.ref, "-dbtype", "nucl", "-out", g.ref]
            with open("/dev/null") as null:
                blastdb_return = call(blast_db_cmd, stdout=null)
            #find the similar genomes
            sim = [bitch for bitch in left if NIC_similarity(bitch.ref,g.ref,blast_db = False, chunk_size = 250, length_threshold=0.8, identity_threshold=0.75
                                                             ) > cutoff]
            #remove them for the list of left genomes
            for s in sim:
                left.remove(s)
            for f in blast_db_files:
                os.remove(f)
            print g.name+"\t"+":".join([ s.name for s in sim])
            clusters[g.name] = sim
            for k in sim:
                k.cluster = g.name
                k.write_data()
            with open(output,"w") as handle:
                json.dump(output, handle)
    return clusters

def pfam_clusters(genomes,output):
    os.system("hmmsearch --cpu 16 -o hmmer_good_prots.raw ~/glob/data/pfam/Pfam-A.hmm  analyses/orthoMCL/goodProteins.fasta")
    return clusters



def checkm(genomes, output, cpus = 1):
    try: 
        call(["which", "hmmscan"])
    except: 
        print "No hmmscan"
    genome_dir = output + "_temp"
    checkm_dir = output + "_checkm"
    if os.path.exists(checkm_dir):
        os.removedirs(checkm_dir)
    if os.path.exists(genome_dir):
        os.removedirs(genome_dir)
    os.makedirs(genome_dir)
    for g in genomes:
        call(["cp", g.genome, genome_dir])

    with open(output,"w") as handle: 
        call(["checkm", "lineage_wf", "-t", str(cpus), "-x", "fna", genome_dir, checkm_dir], stdout=handle)
        
def phylophlan(genomes, output, cpus = 1, phylophlan_folder = "/home/moritz/repos/phylophlan/", clean=False, default_genomes = [], full = False, proj_name = "temp"):
    cwd = os.getcwd()
    os.chdir(phylophlan_folder)
    in_dir = pjoin(phylophlan_folder, "input", proj_name)
#    os.makedirs(in_dir)
#    for g in genomes:
#        call(["cp", g.proteom, in_dir])
#    for g in default_genomes:
#        call(["cp", g, in_dir])
#    for l in os.listdir(in_dir):
#        os.system("sed -i 's/*//g' " + pjoin(in_dir,l))
    if full:
        call(["python", "phylophlan.py", "--nproc", str(cpus), "-i", "-t", proj_name])
    else :
        call(["python", "phylophlan.py", "--nproc", str(cpus), "-u", proj_name])
    os.chdir(cwd)
    call(["cp",pjoin(phylophlan_folder,"output",proj_name, proj_name + ".tree.nwk"), output])
    if clean:
        os.removedirs(in_dir)
        os.removedirs(pjoin(phylophlan_folder,"output", proj_name))
        os.removedirs(pjoin(phylophlan_folder,"data", proj_name))
        
def parse_checkm_results(genome, checkm_out):
    df = read_table(checkm_out, sep=r"\s{2,}", skipinitialspace=True, comment="-", index_col=0)
    line = df.loc[genome.name].to_dict()
    genome.checkm_meta = line 
    genome.write_data()
