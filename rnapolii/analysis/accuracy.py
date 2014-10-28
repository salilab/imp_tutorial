import IMP
import IMP.pmi
import IMP.pmi.analysis
import IMP.pmi.output
import IMP.atom
import glob
import itertools

rmfs=glob.glob('kmeans_500_1/cluster.0/*.rmf3')[0::10]


selection_dictionary={"Rpb1":[(1,1176,"Rpb1")],
                  "Rpb2":["Rpb2"],
                  "Rpb3":["Rpb3"],
                  "Rpb4":["Rpb4"],
                  "Rpb5":["Rpb5"],
                  "Rpb6":["Rpb6"],
                  "Rpb7":["Rpb7"],
                  "Rpb10":["Rpb10"],
                  "Rpb11":["Rpb11"],
                  "Rpb12":["Rpb12"]}

# accuracy of every pair
#for c in itertools.combinations(selection_dictionary.keys(),2):
#    selection_dictionary.update({c[0]+"_"+c[1]:selection_dictionary[c[0]]+selection_dictionary[c[1]]})

# accuracy of all complex
selection_dictionary.update({"All_solid":[(1,1176,"Rpb1"),"Rpb2","Rpb3","Rpb4","Rpb5","Rpb6","Rpb7","Rpb10","Rpb11","Rpb12"]})




model=IMP.Model()


frames=[0]*len(rmfs)

model=IMP.Model()
pr=IMP.pmi.analysis.Precision(model,'one',selection_dictionary=selection_dictionary)
pr.set_precision_style('pairwise_rmsd')

pr.add_structures(zip(rmfs,frames),"ALL",is_mpi=False)



refrmf=''
pr.set_reference_structure("../create_native_rmf/native.rmf3",0)

print pr.get_average_distance_wrt_reference_structure("ALL")
