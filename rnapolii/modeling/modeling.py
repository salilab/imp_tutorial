"""
#############################################
##  IMP Tutorial Script
##
##  PMI2 
#############################################
#
# Short modeling script combining EM and Crosslinking data
# to localize two domains of RNA Polymerase II
#
# Authors: Riccardo Pellarin, Charles Greenberg, Daniel Saltzberg
#
# References: 
#
# General IMP paper
#
#
"""
import IMP
import IMP.core
import IMP.pmi.restraints.crosslinking
import IMP.pmi.restraints.stereochemistry
import IMP.pmi.restraints.em
import IMP.pmi.representation
import IMP.pmi.tools
import IMP.pmi.macros
import IMP.pmi.topology

import os
import sys


#---------------------------
# Define Input/Output Files
#---------------------------
datadirectory = "../data/"
topology_file = datadirectory+"topology_pmi2.txt" 
target_gmm_file = datadirectory+'emd_1883.map.mrc.gmm.50.txt' # The EM map data
output_directory = "./output"

#--------------------------
# Set MC Sampling Parameters
#--------------------------
num_frames = 1000
if '--test' in sys.argv: num_frames=50
num_mc_steps = 10

#--------------------------
# Mover Parameters
#--------------------------

# rigid body movement params. 
# These should be optimized according to MC acceptance ratios.
rb_max_trans = 4.0
rb_max_rot = 0.04
srb_max_trans = 2.0

# flexible bead movement
bead_max_trans = 4.0

# each list contains list of domain names (from topology) that move together
#  flexible beads are automatically added to missing regions and sampled
#  These are defined in the topology file.
rigid_bodies = [["Rpb4"],
                ["Rpb7"]]
#super_rigid_bodies = [["Rpb4","Rpb7"]]
#chain_of_super_rigid_bodies = [["Rpb4"],
#                               ["Rpb7"]]
#
################################################
#

#--------------------------------
# REPRESENTATION
# Build the Model Representation Using Topology File
#
# The topology file has the format:
#|molecule_name|color|fasta_fn|fasta_id|pdb_fn|chain|residue_range|pdb_offset|bead_size|em_residues_per_gaussian|rigid_body|super_rigid_body|chain_of_super_rigid_bodies|flags|
#|Rpb1   |blue  |1WCM.fasta.txt|1WCM:A|1WCM_map_fitted.pdb|A|1,1140   |0|10|0||||
#|Rpb1   |blue  |1WCM.fasta.txt|1WCM:A|1WCM_map_fitted.pdb|A|1141,1274|0|10|0||||
#...
#|Rpb4   |yellow|1WCM.fasta.txt|1WCM:D|1WCM_map_fitted.pdb|D|all      |0|20|40|1|1,2|||
#--------------------------------

# Initialize model
m = IMP.Model()

# Read in the topology file.  
# Specify the directory wheere the PDB files, fasta files and GMM files are
topology = IMP.pmi.topology.TopologyReader(topology_file, 
                                  pdb_dir=datadirectory, 
                                  fasta_dir=datadirectory, 
                                  gmm_dir=datadirectory)

# Use the BuildSystem macro to build states from the topology file
bs = IMP.pmi.macros.BuildSystem(m)

# Each state can be specified by a topology file.
bs.add_state(topology)

# Once all of your states are added, execute the macro.
# This will return the root hierarchy (root_hier) and the dof object, both of which are used later on.
root_hier, dof = bs.execute_macro(max_rb_trans=rb_max_trans, 
                                  max_rb_rot=rb_max_rot, 
                                  max_bead_trans=bead_max_trans, 
                                  max_srb_trans=srb_max_trans)
sys = bs.system

# This nomenclature to get a list of all the molecule objects is bad...make this simpler.
all_molecules = list(bs.get_molecules()[0].values())

shuffled_molecules=[]
for mol in rigid_bodies:
    shuffled_molecules.append(bs.get_molecule(mol[0]))




# Shuffle the configuration of only the molecules we are interested in (Rpb4 and Rpb7)
IMP.pmi.tools.shuffle_configuration(shuffled_molecules, 
                                    max_translation=50, 
                                    verbose=True,
                                    cutoff=5.0,
                                    hierarchies_included_in_collision=list(bs.get_molecules()[0].values()),
                                    niterations=100)



#-----------------------------------
# SCORING FUNCTION / RESTRAINTS
#
# Here we define a number of restraints on our system. The sum of all of these restraints
# is our scoring function. 
#
#  For all restraints, calling add_to_model() incorporates them into the scoring function
#  Appending the restraints to the outputobjects list reports them in the stat file.
#-----------------------------------
outputobjects = [] # reporter objects...output is included in the stat file


#-------
# Excluded Volume Restraint
#  Since much of our system is fixed, we initialize a bipartite restraint
#  that is only evaluated between included_objects and other_objects. 
#
ev = IMP.pmi.restraints.stereochemistry.ExcludedVolumeSphere(
                                         included_objects=shuffled_molecules, 
                                         other_objects=all_molecules,
                                         resolution=10)
# To apply this restraint to all molecules:
#ev = IMP.pmi.restraints.stereochemistry.ExcludedVolumeSphere(included_objects=all_molecules, resolution=10)
ev.add_to_model()         # add to scoring function
outputobjects.append(ev)  # add to output


#-------
# Crosslinks - dataset 1
#  To use this restraint we have to first define the data format.  
# 
# This data file looks like:
#
# prot1,res1,prot2,res2
# Rpb1,34,Rpb1,49
# Rpb1,101,Rpb1,143
# ...
#
# We then initialize a CrossLinkDataBase that uses a keywords converter to map column to information.
# The required fields are the protein and residue number for each side of the crosslink.
xldbkwc = IMP.pmi.io.crosslink.CrossLinkDataBaseKeywordsConverter()
xldbkwc.set_protein1_key("prot1")
xldbkwc.set_protein2_key("prot2")
xldbkwc.set_residue1_key("res1")
xldbkwc.set_residue2_key("res2")

xl1 = IMP.pmi.io.crosslink.CrossLinkDataBase(xldbkwc)
xl1.create_set_from_file(datadirectory+'polii_juri.csv')

# Now, we set up the restraint.
xl1rest = IMP.pmi.restraints.crosslinking.CrossLinkingMassSpectrometryRestraint(
                                   root_hier=root_hier,  # The root hierarchy
                                   CrossLinkDataBase=xl1,# The XLDB defined above
                                   length=21.0,          # Length of the linker in angstroms
                                   slope=0.01,           # A linear term that biases XLed residues together
                                   resolution=1.0,       # Resolution at which to apply the restraint. Either 1 (residue) or 0 (atomic)
                                   label="Chen",         # Used to label output in the stat file
                                   weight=1.)            # Weight applied to all crosslinks in this dataset
xl1rest.add_to_model()
outputobjects.append(xl1rest)


#-------
# Crosslinks - dataset 2
#  We can easily add a second set of crosslinks.
#  These have a different format and label, but other settings are the same

xldbkwc = IMP.pmi.io.crosslink.CrossLinkDataBaseKeywordsConverter()
xldbkwc.set_protein1_key("pep1.accession")
xldbkwc.set_protein2_key("pep2.accession")
xldbkwc.set_residue1_key("pep1.xlinked_aa")
xldbkwc.set_residue2_key("pep2.xlinked_aa")

xl2 = IMP.pmi.io.crosslink.CrossLinkDataBase(xldbkwc)
xl2.create_set_from_file(datadirectory+'polii_xlinks.csv')

xl2rest = IMP.pmi.restraints.crosslinking.CrossLinkingMassSpectrometryRestraint(
                                   root_hier=root_hier,
                                   CrossLinkDataBase=xl1,
                                   length=21.0,
                                   slope=0.01,
                                   resolution=1.0,
                                   label="Trnka",
                                   weight=1.)

xl2rest.add_to_model()          
outputobjects.append(xl2rest)


# Electron Microscopy Restraint
#  The GaussianEMRestraint uses a density overlap function to compare model to data
#   First the EM map is approximated with a Gaussian Mixture Model (done separately)
#   Second, the components of the model are represented with Gaussians (forming the model GMM)


# First, get the model density objects that will be fitted to the EM density.
em_components = IMP.pmi.tools.get_densities(root_hier)

gemt = IMP.pmi.restraints.em.GaussianEMRestraint(em_components,
                                      target_gmm_file,  # EM map GMM file
                                      scale_target_to_mass=True,  # True if the mass of the map and model are identical.
                                      slope=0.0000001,  # A small funneling force pulling towards the center of the EM density.
                                      weight=80.0)           
gemt.add_to_model()
outputobjects.append(gemt)

#--------------------------
# Monte-Carlo Sampling
#--------------------------
# This object defines all components to be sampled as well as the sampling protocol
mc1=IMP.pmi.macros.ReplicaExchange0(m,
              root_hier=root_hier,                         # The root hierarchy
              monte_carlo_sample_objects=dof.get_movers(), # All moving particles and parameters
              output_objects=outputobjects,                # Objects to put into the stat file
              crosslink_restraints=[xl1rest,xl2rest],      # allows XLs to be drawn in the RMF files
              monte_carlo_temperature=1.0,                 
              simulated_annealing=True,
              simulated_annealing_minimum_temperature=1.0,
              simulated_annealing_maximum_temperature=2.5,
              simulated_annealing_minimum_temperature_nframes=200,
              simulated_annealing_maximum_temperature_nframes=20,
              replica_exchange_minimum_temperature=1.0,
              replica_exchange_maximum_temperature=2.5,
              number_of_best_scoring_models=10,
              monte_carlo_steps=num_mc_steps,
              number_of_frames=num_frames,
              global_output_directory=output_directory)

# Start Sampling
mc1.execute_macro()

