def json2nexus(filename, out_name, trait_file, updated_json = None):
    with open(filename) as handle:
        data=json.loads(handle.read())

    genomes = list(set(sum([d['genes'].keys() for d in data],[])))

    trait_vectors = {g: [len(c['genes'][g]) if( c['genes'].has_key(g) ) else 0 for c in data] for g in genomes }
    
    formated_trait_vector = {g : "".join([hex(15 if t >15 else t)[-1].capitalize() for t in ts]) for g,ts in trait_vectors.iteritems()}

    with open(out_name, "w") as handle:
        for g, ts in formated_trait_vector.iteritems():
            handle.write(">" + g + "\n")
            handle.write(ts + "\n")

    zf = len(str(len(data)))
    with open(trait_file,"w") as handle:
        handle.write(
"""#NEXUS
BEGIN CHARACTERS;
DIMENSIONS NCHAR=%d;
CHARSTATELABELS""" %(len(data)))
        for i,d in enumerate(data):
            handle.write("\n%d %s," %(i, "orthoMCL_"+str(i).zfill(zf)))
            d['name'] = "orthoMCL_"+str(i+1).zfill(zf)
        handle.write("\b;\nEND;\n")

    if updated_json:
        with open(updated_json,"w") as handle:
            json.dump(data, handle,  indent=4, sort_keys=True)

def mergeJSONwPAUP(json_file, paup_file):
    with open(json_file) as handle:
        data=json.loads(handle.read())
    changes = DataFrame.from_csv(paup_file, sep="\t", index_col=None, header=0)
    for d in tqdm(data):
        if len(d['genes']) > 1:
            sub = changes[changes['orthoMCL'] == (int(d['name'].split("_")[1])+1)]
            d['ancestral'] = {r[1] : {'from' : r[3], 'to': r[4]} for r in sub.itertuples()}
    return data

def make_bmft(data, genomes_of_interest, ancestry, clustering):
    data_sub = [d for d in data if any([g in genomes_of_interest for g in d['genes'].keys()] ) ]
    head = ["ID", "hypothetical.name", "name.confidence"]+  genomes_of_interest+ [ "ancestor" , ancestry , "genomes.of.interest", "other.genomes", "genes"]
    lines = []
    genomes = set(clustering.gene2genome.values())
    for d in data_sub:
        name = d['name'] 
        hname = d['annotation']
        c_name = d['annot_fraction']
        if d.has_key('ancestral') and d['ancestral'].has_key(ancestry):
            ancestor = d['ancestral'][ancestry]
            frmo = ancestor['from']
            diff = ancestor['to']-ancestor['from']
        else :
            frmo = "NA"
            diff = 0
        counts = [len(d['genes'][g]) if d['genes'].has_key(g) else 0 for g in genomes_of_interest]
        icounts = len([g for g in genomes if g in genomes_of_interest and d['genes'].has_key(g)])
        ocounts = len([g for g in genomes if g not in genomes_of_interest and d['genes'].has_key(g)])
        genes = ";".join(d['mapping'].keys())
        line = [ name, hname, c_name ] + counts + [ frmo, diff , icounts, ocounts, genes ] 
        lines +=  [[str(l) for l in line]]
    bmft = DataFrame.from_records(lines, columns=head, index="ID")
    return bmft
