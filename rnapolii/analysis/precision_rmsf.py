import IMP
import IMP.pmi
import IMP.pmi.analysis
import IMP.pmi.output
import IMP.atom
import glob
import itertools

is_mpi=False
test_mode=True  # runs on the first 10 structures to test if it runs smoothly

# specify the cluster directory to be analysed

root_cluster_directory='kmeans_50_1'

# choose whatever selection for the precision calculation

selection_dictionary={"Rpb4":["Rpb4"],
                  "Rpb7":["Rpb7"],
                  "Rpb4_Rpb7":["Rpb4","Rpb7"]}


#####################################################################
# don't change anything below
rmfs=[]
frames=[]
clusterdirectories=glob.glob(root_cluster_directory+'/cluster.*/')

if test_mode:
  # runs on the first 10 structures to test if it runs smoothly
  for clusterdirectory in clusterdirectories:
      rmfs.append(glob.glob(clusterdirectory+'/*.rmf3')[0::10])
      frames.append([0]*len(rmfs[-1]))
else:
  for clusterdirectory in clusterdirectories:
      rmfs.append(glob.glob(clusterdirectory+'/*.rmf3')[0::1])
      frames.append([0]*len(rmfs[-1]))
 
model=IMP.Model()
pr=IMP.pmi.analysis.Precision(model,'one',selection_dictionary=selection_dictionary)
pr.set_precision_style('pairwise_rmsd')

for n in range(len(rmfs)):
    pr.add_structures(zip(rmfs[n],frames[n]),clusterdirectories[n],is_mpi=is_mpi)


for pair in itertools.product(range(len(rmfs)), repeat=2):
    clus1=pair[0]
    clus2=pair[1]
    pr.get_precision(root_cluster_directory+"/precision."+str(clus1)+"."+str(clus2)+".out",
                     clusterdirectories[clus1],
                     clusterdirectories[clus2],
                     is_mpi=is_mpi,skip=1)

for n in range(len(rmfs)):
    pr.get_rmsf(clusterdirectories[n],clusterdirectories[n]+"/",is_mpi=is_mpi,skip=1)
