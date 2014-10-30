import IMP
import IMP.pmi
import IMP.pmi.macros
import sys

is_mpi=False

model=IMP.Model()

# initialize the macro

mc=IMP.pmi.macros.AnalysisReplicaExchange0(model,
                                        stat_file_name_suffix="stat",     # don't change
                                        merge_directories=["../modeling/"], # change this list splitting the runs or adding new runs
                                        global_output_directory="./output/", # don't change
                                        rmf_dir="rmfs/")  # don't change

# fields that have to be extracted for the stat file

feature_list=["ISDCrossLinkMS_Distance_intrarb",
              "ISDCrossLinkMS_Distance_interrb",
              "ISDCrossLinkMS_Data_Score",
              "GaussianEMRestraint_None",
              "SimplifiedModel_Linker_Score_None",
              "ISDCrossLinkMS_Psi",
              "ISDCrossLinkMS_Sigma"]

# Dictionary of densities to be calculated
# the key is the name of the file and the value if the selection
# example:
#              {"med17-CTD":[(200,300,"med17")],"med17-CTD.med14":[(200,300,"med17"),"med14"]   }

# list of component names needed to calculate the RMSD for the clustering

components_names=["Rpb4",
                  "Rpb7"]

density_names={}
rmsd_names={}
for name in components_names:
    density_names[name]=[name]
    rmsd_names[name]=name

nclusters=1                                        # number of clusters needed by kmeans
mc.clustering("SimplifiedModel_Total_Score_None",  # don't change, field where to find the score
              "rmf_file",                          # don't change, field where to find the path for the rmf_file
              "rmf_frame_index",                   # don't change, field for the frame index
              prefiltervalue=2900.0,              # prefilter the models by score
              number_of_best_scoring_models=50,   # number of models to be clustered
              alignment_components=None,           # don't change, (list of proteins you want to use for structural alignment
              rmsd_calculation_components=rmsd_names, # list of proteins used to calculated the rmsd
              distance_matrix_file="distance.rawmatrix.pkl", # save the distance matrix
              outputdir="kmeans_50_"+str(nclusters)+"/",  # directory name for the clustering
              feature_keys=feature_list,                     # extract these fields from the stat file
              load_distance_matrix_file=False,                # skip the matrix calcuklation and read the precalculated matrix
              skip_clustering=False,                         # skip clustering
              display_plot=False,                            # display the heat map plot of the distance matrix
              exit_after_display=False,                      # exit after having displayed the distance matrix plot
              get_every=1,                                   # skip structures for faster computation
              is_mpi=is_mpi,                                 # mpi enabled
              number_of_clusters=nclusters,                  # number of clusters to be used by kmeans algorithm
              voxel_size=3.0,                                # voxel size of the mrc files
              density_custom_ranges=density_names)    # setup the list of densities to be calculated

