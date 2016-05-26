def renaming_tree(tree, outputtree, mapping, keep = False):
    with open(tree) as handle:
        t = handle.readlines()
        if len(t) > 0:
            t=t[0][:-1]
        else :
            print tree, "is empty for some reason"
            t = ""
            
    if keep:
        mapping = { k: k + "|" + v for k,v in mapping.iteritems()}
        
    for k,v in mapping.iteritems():
        t = t.replace(k + ":", mapping[k] + ":")
        t = t.replace(k + ",", mapping[k] + ",")
    with open(outputtree,"w") as handle:
        handle.writelines(t + "\n")


