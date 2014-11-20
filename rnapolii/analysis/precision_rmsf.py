import IMP
import IMP.pmi
import IMP.pmi.analysis
import IMP.pmi.output
import IMP.atom
import glob
import itertools

# common settings
test_mode = False                        # runs on every 10 frames
root_cluster_directory = 'kmeans_5_1'  # specify the directory (of clusters) to be analysed


# choose components for the precision calculation
# key is the named precision item
# value is a list of components
selections={"Rpb4":["Rpb4"],
            "Rpb7":["Rpb7"],
            "Rpb4_Rpb7":["Rpb4","Rpb7"]}


##############################
# don't change anything below
##############################

# setup Precision calculator
model = IMP.Model()
pr = IMP.pmi.analysis.Precision(model,'one',selection_dictionary=selections)
pr.set_precision_style('pairwise_rmsd')


# gather the RMF filenames for each cluster
rmfs=[]
frames=[]
cluster_dirs=glob.glob(root_cluster_directory+'/cluster.*/')
if test_mode:
  # runs on the first 10 structures to test if it runs smoothly
  for d in cluster_dirs:
      rmfs.append(glob.glob(d+'/*.rmf3')[0::10])
      frames.append([0]*len(rmfs[-1]))
else:
  for d in cluster_dirs:
      rmfs.append(glob.glob(d+'/*.rmf3')[0::1])
      frames.append([0]*len(rmfs[-1]))

# add them to the Precision object
for n in range(len(rmfs)):
    pr.add_structures(zip(rmfs[n],frames[n]),cluster_dirs[n])

# calculate intra-cluster and inter-cluster precision
print "calculating precision"
for clus1,clus2 in itertools.combinations_with_replacement(range(len(rmfs)),2):
    pr.get_precision(root_cluster_directory+"/precision."+str(clus1)+"."+str(clus2)+".out",
                     cluster_dirs[clus1],
                     cluster_dirs[clus2],
                     skip=1)

# compute resiude mean-square fluctuation (RMSF)
print "calculating RMSF"
for n in range(len(rmfs)):
    pr.get_rmsf(cluster_dirs[n],cluster_dirs[n]+"/",skip=1)

print "done"
