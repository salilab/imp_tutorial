import sys
import matplotlib as mpl
mpl.use('Agg')

import IMP  # noqa: E402
import IMP.pmi  # noqa: E402
import IMP.pmi.macros  # noqa: E402

# most common settings
num_clusters = 1
num_top_models = 5
merge_directories = ["../modeling_em/"]
prefiltervalue = 2900.0
out_dir = "kmeans_%i_%i/" % (num_top_models, num_clusters)
if '--test' in sys.argv:
    prefiltervalue = 8000.0

#################################
# should not have to change below
##################################

model = IMP.Model()

# initialize the macro
mc = IMP.pmi.macros.AnalysisReplicaExchange0(
    model, merge_directories=merge_directories)

# fields that have to be extracted for the stat file
feature_list = ["ISDCrossLinkMS_Distance_intrarb",
                "ISDCrossLinkMS_Distance_interrb",
                "ISDCrossLinkMS_Data_Score",
                "GaussianEMRestraint_None",
                "SimplifiedModel_Linker_Score_None",
                "ISDCrossLinkMS_Psi",
                "ISDCrossLinkMS_Sigma"]

# Dictionary of densities to be calculated
# the key is the name of the file and the value is the selection
# example: {"med17-CTD": [(200,300,"med17")],
#           "med17-CTD.med14": [(200,300,"med17"),"med14"]}
density_names = {"Rpb4": ["Rpb4"],
                 "Rpb7": ["Rpb7"]}

# list of component names needed to calculate the RMSD for the clustering
rmsd_names = {"Rpb4": "Rpb4",
              "Rpb7": "Rpb7"}

# components used for structural alignment
align_names = None  # (None because EM provides reference frame)

mc.clustering(
    # prefilter the models by score
    prefiltervalue=prefiltervalue,
    # number of models to be clustered
    number_of_best_scoring_models=num_top_models,
    # list of proteins you want to use for structural alignment
    alignment_components=None,
    # list of proteins used to calculated the rmsd
    rmsd_calculation_components=rmsd_names,
    # save the distance matrix
    distance_matrix_file="distance.rawmatrix.pkl",
    # location for clustering results
    outputdir=out_dir,
    # extract these fields from the stat file
    feature_keys=feature_list,
    # skip the matrix calculation and read the precalculated matrix
    load_distance_matrix_file=False,
    # display the heat map plot of the distance matrix
    display_plot=True,
    # exit after having displayed the distance matrix plot
    exit_after_display=False,
    # skip structures for faster computation
    get_every=1,
    # number of clusters to be used by kmeans algorithm
    number_of_clusters=num_clusters,
    # voxel size of the mrc files
    voxel_size=3.0,
    # setup the list of densities to be calculated
    density_custom_ranges=density_names)
