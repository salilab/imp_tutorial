#!/usr/bin/env python

'''accuracy.py
This script calculates the average accuracy of the structures in a cluster
It uses the IMP.pmi.analysis.Precision class
Requires a native structure for comparison
'''

from __future__ import print_function
import IMP
import IMP.pmi
import IMP.pmi.analysis
import IMP.pmi.output
import IMP.atom
import glob
import itertools

# common settings
reference_rmf = "../data/native.rmf3"
test_mode = False                             # run on every 10 rmf files
rmfs = glob.glob('kmeans_*_*/cluster.0/*.rmf3') # list of the RMFS to calculate on

# choose components for the precision calculation
# key is the named precision item
# value is a list of selection tuples [either "domain_name" or (start,stop,"domain_name") ]
selections = {"Rpb4":["Rpb4"],
              "Rpb7":["Rpb7"],
              "Rpb4_Rpb7":["Rpb4","Rpb7"]}


##############################
# don't change anything below
##############################

# setup Precision calculator
model=IMP.Model()
frames=[0]*len(rmfs)
pr=IMP.pmi.analysis.Precision(model,selection_dictionary=selections)
pr.set_precision_style('pairwise_rmsd')
pr.add_structures(zip(rmfs,frames),"ALL")

# calculate average distance to the reference file
pr.set_reference_structure(reference_rmf,0)
print(pr.get_average_distance_wrt_reference_structure("ALL"))
