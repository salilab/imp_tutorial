from __future__ import print_function
import IMP
import IMP.pmi
import IMP.pmi.macros
import matplotlib.pyplot as plt
plot = plt.plot
figure = plt.figure
boxplot = plt.boxplot
savefig= plt.savefig
m=IMP.Model()
import warnings
warnings.filterwarnings('ignore')
# ## Stage 4 - Analysis <a name="Analysis_3"></a>
# 
# In the analysis stage we cluster (group by similarity) the sampled models to determine high-probability configurations. Comparing clusters may indicate that there are multiple acceptable configurations given the data. 
# 
# In this stage we perform several analysis.  Here, we will perform calculations for:
# 
# * **Clustering**: Grouping the structure together using similarity via RMSD
# * **Cluster Uncertainty**: Determining the within-group precision and between-group similarity via RMSD
# * **Cluster Accuracy**: Fit of the calculated clusters to the true (known) solution
# * **Sampling Exhaustiveness**: Qualitative and quantitative measurement of sampling completeness
# 
# ### Precomputed results <a name="Precomputed_Results_3"></a>
# 
# A long modeling run was precomputed and analyzed. You can [download](ftp://salilab.org/tutorials/imp/rnapolii/results.tar.gz) it from our website, and you can [download](ftp://salilab.org/tutorials/imp/rnapolii/analysis.tar.gz) the corresponding analysis.
# 
# ### Clustering top models using `analysis.py` <a name="Clustering_3"></a>
# The `long_analysis.py` script, found in the `modeling` directory, calls the [AnalysisReplicaExchange](https://integrativemodeling.org/nightly/doc/ref/classIMP_1_1pmi_1_1macros_1_1AnalysisReplicaExchange.html) macro, which finds top-scoring models, extracts coordinates, runs clustering, and does basic cluster analysis including creating localization densities for each subunit. The script generates RMF, MRC files which should be viewable in Chimera.
# 
# We can choose the number of clusters by changing the distance threshold, the subunits we want to use to calculate the RMSD, and the number of good-scoring solutions to include.

# If we perform sampling multiple times separately, they can all be analyzed at the same time by appending to list of stat files. The `best_models` parameter set the number of best scoring models to be analyzed. Note that we use `alignment=False`. This is needed in case there is no absolute reference frame (like an EM map).

# In[ ]:


are=IMP.pmi.macros.AnalysisReplicaExchange(m,
                 ["./output/rmfs/0.rmf3"],
                 best_models=100,
                 alignment=False)

print(are)


# Then, we start the clustering. 
# We specify the components used in calculating the RMSD between models. 
# Then we cluster using a rmsd threshold of 30 Angstroms.

# In[ ]:


are.set_rmsd_selection(molecules=["Rpb4","Rpb7"])

are.cluster(30.0)


# For each cluster, we can print its information

# In[ ]:


print(are)

for cluster in are:
    print(cluster)


# We can get a given cluster by using the square bracket, as in lists, for instance `are[1]` is the cluster with index 1. Using the list properties of the object `are`, we get the best scoring cluster.

# In[ ]:


from operator import attrgetter
    
best_cluster=min(are,key=attrgetter('average_score'))

print(best_cluster)


# We can iterate on the members of the cluster to display the infos. Afterwords, we save the coordinates of the cluster in a rmf file.

# In[ ]:


for member in best_cluster:
    print(member)
    
are.save_coordinates(best_cluster)


# Next we can examine the distances between all cluster members. A plot is output to a single file in the clustering directory. The first plot is the distance matrix of the models after being grouped into clusters. 
# 
# The second plot is a dendrogram, basically showing the distance matrix in a hierarchical way. Each vertical line from the bottom is a model, and the horizontal lines show the RMSD agreement between models. Sometimes the dendrogram can indicate a natural number of clusters, which can help determine the correct threshold to use. 
# 
# <img src="files/images/rnapolii_dist_matrix.png" alt="Distance matrix and dendrogram" width="600px" />
# 

# In[ ]:


# slow!
are.plot_rmsd_matrix("rmsd_matrix.pdf")


# ### Structural uncertainty of the solutions <a name="uncertainty_3"></a>
# 
# The cluster center can be computed as the median structure. After that one can compute the precision of the cluster, as well as the average distance between two clusters.

# In[ ]:


are.compute_cluster_center(cluster=best_cluster)
print(are.precision(cluster=best_cluster))
print(are.bipartite_precision(cluster1=best_cluster,cluster2=are[0]))


# We can plot the root mean square fluctuation (rmsf) of a given molecule in a given cluster.

# In[ ]:


rmsf1=are.rmsf(cluster=best_cluster,molecule='Rpb4');
plot(rmsf1.keys(),rmsf1.values())
figure()

rmsf2=are.rmsf(cluster=best_cluster,molecule='Rpb7');
plot(rmsf2.keys(),rmsf2.values())


# And compute the the rmsf for all molecules and map the value on the structure. Finally we save the colored cooridnates in a rmf file.
# 
# <img src="files/images/rnapolii_rmsf.all.png" alt="Structural uncertainty" width="600px" />

# In[ ]:


for mol in ['Rpb1','Rpb2','Rpb3','Rpb4','Rpb5','Rpb6','Rpb7','Rpb8','Rpb9','Rpb10','Rpb11','Rpb12']: 
    are.rmsf(cluster=best_cluster,molecule=mol);
ch1=IMP.pmi.tools.ColorHierarchy(are.stath1)
ch1.color_by_uncertainty()
are.save_coordinates(best_cluster)


# We can save the localization densities of a given cluster, for given groups of molecules.
# Now we specify the subunits (or groups or fractions of subunits) for which we want to create density localization maps. 
# `density_names` is a dictionary, where the keys are convenient names like "Rpb4" and the values are a list of selections. 
# The selection items can either be a domain name like "Rpb1" or a list like (200,300,"Rpb1") 
# which means residues 200-300 of component Rpb1. This enables the user to combine multiple selections 
# for a single density calculation.
# 
# The localization densities can give a qualitative idea of the precision of a cluster. Below we show results from `cluster.1` in the provided results: the native structure without Rpb4/7 (in blue), the target density map (in mesh), and the localization densities (Rpb4 in cyan, Rpb7 in purple). The localizations are quite narrow and close to the native solution:
# 
# <img src="files/images/rnapolii_localization.png" alt="Localization densities" width="600px" />

# In[ ]:


density_names={"REST":["Rpb1","Rpb2","Rpb3","Rpb5","Rpb6","Rpb8","Rpb9","Rpb10","Rpb11","Rpb12"],
     "Rpb4":["Rpb4"],
     "Rpb7":["Rpb7"],
     "Rpb1-200-300":[(200,300,"Rpb1")]}

are.save_densities(cluster=best_cluster,density_custom_ranges=density_names,prefix="BestCluster")
are.save_densities(cluster=are[0],density_custom_ranges=density_names,prefix="Cluster-0")


# We can compute the global contact map of the whole complex for the second cluster.

# In[ ]:


# it is slow
are.contact_map(cluster=are[1]);


# ### Accuracy evaluation <a name="Accuracy_3"></a>
# We can evaluate the accuracy of a cluster against a native configuration. This is useful for benchmarking (but obviously is of no use when we don't know the 'real' structure).
# First we add the native structure as an independent cluster

# In[ ]:


are.add_cluster(["../rnapolii/data/native.rmf3"])

print(are)


# Then we compute the `bipartite_precision` between a given cluster and the native structure.

# In[ ]:


print(are[-1].members)
are.bipartite_precision(cluster1=best_cluster,cluster2=are[-1])


# ### Sampling Exhaustiveness <a name="Sampling_Exhaustiveness_3"></a>
# We can also determine sampling exhaustiveness by dividing the models into multiple sets, performing clustering on each set separately, and comparing the clusters. This step is left as an exercise to the reader. Some things you can try:
# * cluster two subsets of the data
# * qualitative analysis: look at the localization densities - they should be similar for the two subsets
# * quantitative analysis: compute cross-precision for the clusters. 
# 
# If the sampling is exhaustive, then similar clusters should be obtained from each independent set, and the inter-cluster precision between two equivalent clusters should be very low (that is, there should be a 1:1 correspondence between the two sets of clusters, though the ordering may be different).
