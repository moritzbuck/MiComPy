from pylab import *

shape_star = str(3)
red = "#ff0000"

template_pres_abs = lambda label, color, shape, f_labels, data :"""
DATASET_BINARY

SEPARATOR COMMA

#label is used in the legend table (can be changed later)
DATASET_LABEL,%s

#dataset color (can be changed later)
COLOR,%s
FIELD_SHAPES,%s
FIELD_LABELS,%s

DATA
%s
""" % (label, color, shape, f_labels, data)



template_strips = lambda  label, cm, data :"""
DATASET_COLORSTRIP

SEPARATOR COMMA
DATASET_LABEL,%s
COLOR,#000000
COLOR_BRANCHES,0

LEGEND_TITLE,%s
LEGEND_SHAPES,%s 
LEGEND_COLORS,%s
LEGEND_LABELS,%s


DATA
# taxa_name,"clade",color,label
# eg.: rfile_blablabla,clade,#ff0000,blablablabacteria
%s
""" % (label,label, ",".join(["1"]*len(cm.values())), ",".join([cm[k]for k in sort(cm.keys())]), ",".join(sort(cm.keys())), data  )

template_colours = lambda  data :"""
TREE_COLORS

SEPARATOR COMMA

DATA
# taxa_name,"clade",color,label
# eg.: rfile_blablabla,clade,#ff0000,blablablabacteria
%s
""" % ( data)

template_gradient = lambda label, data, color :"""
DATASET_GRADIENT
SEPARATOR COMMA

DATASET_LABEL,%s
COLOR,%s

COLOR_MIN,#ffffff
COLOR_MAX,%s

DATA
# taxa_name,value
# eg.: rfile_blablabla,200
%s
""" % (label,color,color, data)





float_to_rgb = lambda (r,g,b) : '#%02x%02x%02x' % (int(r*255), int(g*255), int(b*255))

def SAG_data(genomes, file):
    dat = {g.get_meta('short_name') : g.get_meta('type') for g in genomes}
    dat = [ k + "," + str(1 if v == "SAG" else -1) + "\n" for k,v in dat.iteritems()]
    with open(file,"w") as handle:
        handle.writelines(template_pres_abs("SAGs",red,shape_star,"", "".join(dat)))

def taxa_data(genomes, file):
    dat = {g.get_meta('short_name') : g.get_meta('phylum', default = "") + g.get_meta('taxonomy_external', default = "") for g in genomes if g.is_good()}
    factors = list(set(dat.values()))
    cm = get_cmap('Set3')
    color_map = {k : float_to_rgb(cm(1.*i/len(factors))[0:3]) for i,k in enumerate(factors)}
    dat = [ ",".join([k, "range", color_map[v],v])  + "\n" for k,v in dat.iteritems() if v != ""]
    with open(file,"w") as handle:
        handle.writelines(template_colours("".join(dat)))



def env_data(genomes, file):
    dat = {g.metadata['short_name'] : g.get_meta('environment', default = "") for g in genomes if g.is_good()}
    factors = list(set(dat.values()))
    cm = get_cmap('Set1')
    color_map = {k : float_to_rgb(cm(1.*i/len(factors))[0:3]) for i,k in enumerate(factors) if k != "" and k != "unknown"}
    dat = [ ",".join([k, color_map[v],v])  + "\n" for k,v in dat.iteritems() if v != "" and v != "unknown"]
    with open(file,"w") as handle:
        handle.writelines(template_strips("environment",color_map,"".join(dat)))

def clust_data(name, clust_data, file, color="#ff0000"):
    dat = [ k + "," + str(float(v)+0.01) + "\n" for k,v in clust_data.iteritems()]
    with open(file,"w") as handle:
        handle.writelines(template_gradient(name,"".join(dat), color))


def class_data(genomes, class_dict, file , name = "derep_clusters"):
    factors = list(set(class_dict.values()))
    cm = get_cmap('Set1')
    color_map = {k : float_to_rgb(cm(1.*i/len(factors))[0:3]) for i,k in enumerate(factors) if k != "" and k != "unknown"}
    dat = [ ",".join([k, color_map[v],v])  + "\n" for k,v in class_dict.iteritems() if v != "" and v != "unknown"]
    with open(file,"w") as handle:
        handle.writelines(template_strips("environment",color_map,"".join(dat)))

