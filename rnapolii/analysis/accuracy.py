import IMP
import IMP.pmi
import IMP.pmi.analysis
import IMP.pmi.output
import IMP.atom
import glob
import itertools

rmfs=glob.glob('kmeans_*_1/cluster.0/*.rmf3')[0::10]


selection_dictionary={
                  "Rpb4":["Rpb4"],
                  "Rpb7":["Rpb7"],
                  "Rpb4_Rpb7":["Rpb4","Rpb7"]}

model=IMP.Model()

frames=[0]*len(rmfs)

model=IMP.Model()
pr=IMP.pmi.analysis.Precision(model,selection_dictionary=selection_dictionary)
pr.set_precision_style('pairwise_rmsd')

pr.add_structures(zip(rmfs,frames),"ALL")

refrmf=''
pr.set_reference_structure("../data/native.rmf3",0)

print pr.get_average_distance_wrt_reference_structure("ALL")
