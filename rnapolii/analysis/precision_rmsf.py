#!/usr/bin/env python

'''precision_rmsf.py
Calculates the within- and between-cluster RMSD (=precision)
It uses the IMP.pmi.analysis.Precision class
Also calculates within-cluter residue mean square fluctuation (RMSF)
'''

from __future__ import print_function
import IMP
import IMP.pmi
import IMP.pmi.analysis
import IMP.pmi.output
import IMP.atom
import glob
from copy import deepcopy
try:
    from itertools import combinations_with_replacement
except ImportError:
    # Provide implementation (from Python's own itertools docs)
    # for Python < 2.7
    def combinations_with_replacement(iterable, r):
        # combinations_with_replacement('ABC', 2) --> AA AB AC BB BC CC
        pool = tuple(iterable)
        n = len(pool)
        if not n and r:
            return
        indices = [0] * r
        yield tuple(pool[i] for i in indices)
        while True:
            for i in reversed(range(r)):
                if indices[i] != n - 1:
                    break
            else:
                return
            indices[i:] = [indices[i] + 1] * (r - i)
            yield tuple(pool[i] for i in indices)

# common settings
test_mode = False                           # runs on every 20 frames
# specify the directory (of clusters) to be analysed
root_cluster_directory = glob.glob('kmeans_*_*')[0]


# choose components for the precision calculation
# key is the named precision item
# value is a list of selection tuples [either "domain_name"
# or (start,stop,"domain_name") ]
selections = {"Rpb4": ["Rpb4"],
              "Rpb7": ["Rpb7"],
              "Rpb4_Rpb7": ["Rpb4", "Rpb7"]}


##############################
# don't change anything below
##############################

# setup Precision calculator
orig_selections = deepcopy(selections)
model = IMP.Model()
pr = IMP.pmi.analysis.Precision(model, resolution=1,
                                selection_dictionary=selections)
pr.set_precision_style('pairwise_rmsd')

# gather the RMF filenames for each cluster
rmf_list = []
frame_list = []
cluster_dirs = glob.glob(root_cluster_directory+'/cluster.*/')
if test_mode:
    # runs on the first 10 structures to test if it runs smoothly
    for d in cluster_dirs:
        rmf_list.append(glob.glob(d+'/*.rmf3')[0::20])
        frame_list.append([0]*len(rmf_list[-1]))
else:
    for d in cluster_dirs:
        rmf_list.append(glob.glob(d+'/*.rmf3'))
        frame_list.append([0]*len(rmf_list[-1]))

print(frame_list)

# add them to the Precision object
for rmfs, frames, cdir in zip(rmf_list, frame_list, cluster_dirs):
    pr.add_structures(zip(rmfs, frames), cdir)


# calculate intra-cluster and inter-cluster precision
print("calculating precision")
for clus1, clus2 in combinations_with_replacement(range(len(rmf_list)), 2):
    pr.get_precision(
        cluster_dirs[clus1],
        cluster_dirs[clus2],
        root_cluster_directory+"/precision."+str(clus1)+"."+str(clus2)+".out")

# compute residue mean-square fluctuation (RMSF)
print("calculating RMSF")
for d in cluster_dirs:
    pr.get_rmsf(structure_set_name=d, outdir=d)
    pr.selection_dictionary = deepcopy(orig_selections)
print("done")
