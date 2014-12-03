"""
#############################################
##  IMP Tutorial Script
##
##  Riccardo's awesome seminar: Dec 8, 2014
#############################################
#
# Short modeling script combining EM and Crosslinking data
# to localize two domains of RNA Polymerase II
#
# Authors: Riccardo Pellarin, Charles Greenberg, Daniel Saltzberg
#
# References: Papers where this data is shown...
#
"""
import IMP
import IMP.core
import IMP.base
import IMP.algebra
import IMP.atom
import IMP.container

import IMP.pmi.restraints.crosslinking
import IMP.pmi.restraints.stereochemistry
import IMP.pmi.restraints.em
import IMP.pmi.restraints.basic
import IMP.pmi.representation
import IMP.pmi.tools
import IMP.pmi.samplers
import IMP.pmi.output
import IMP.pmi.macros
import IMP.pmi.topology.topology_io

import os
import sys


#---------------------------
# Define Input Files
#---------------------------
datadirectory = "../data/"
topology_file = datadirectory+"topology.txt"
target_gmm_file = datadirectory+'emd_1883.map.mrc.gmm.50.txt'

#--------------------------
# Set MC Sampling Parameters
#--------------------------
num_frames = 20000
if '--test' in sys.argv: num_frames=50
num_mc_steps = 10

#--------------------------
# Create movers
#--------------------------

# rigid body movement params
rb_max_trans = 2.00
rb_max_rot = 0.04

# flexible bead movement
bead_max_trans = 3.00

# each list contains list of domain names (from topology) that move together
#  flexible beads are automatically added to missing regions and sampled
rigid_bodies = [["Rpb4"],
                ["Rpb7"]]
super_rigid_bodies = [["Rpb4","Rpb7"]]
chain_of_super_rigid_bodies = [["Rpb4"],
                               ["Rpb7"]]
#
################################################
#

#--------------------------------
# Build the Model Representation
#--------------------------------

# Initialize model
m = IMP.Model()

# Create list of components from topology file
topology = IMP.pmi.topology.TopologyReader(topology_file)
domains = topology.component_list


bm = IMP.pmi.macros.BuildModel(m,
                    component_topologies=domains,
                    list_of_rigid_bodies=rigid_bodies,
                    list_of_super_rigid_bodies=super_rigid_bodies,
                    chain_of_super_rigid_bodies=chain_of_super_rigid_bodies)
representation = bm.get_representation()

# Randomize the initial configuration before sampling
representation.shuffle_configuration(50)

#--------------------------
# Define Degrees of Freedom
#--------------------------

# Add default mover parameters to simulation
representation.set_rigid_bodies_max_rot(rb_max_rot)
representation.set_floppy_bodies_max_trans(bead_max_trans)
representation.set_rigid_bodies_max_trans(rb_max_trans)

outputobjects = [] # reporter objects (for stat files)
sampleobjects = [] # sampling objects

# Add the movers to the sample and output object lists
outputobjects.append(representation)
sampleobjects.append(representation)

#-----------------------------------
# Define Scoring Function Components
#-----------------------------------

# Here we are defining a number of restraints on our system.
#  For all of them we call add_to_model() so they are incorporated into scoring
#  We also add them to the outputobjects list, so they are reported in stat files


# Excluded Volume Restraint
#  To speed up this expensive restraint, we operate it at resolution 10
ev = IMP.pmi.restraints.stereochemistry.ExcludedVolumeSphere(representation,resolution=10)
ev.add_to_model()
outputobjects.append(ev)


# Crosslinks - dataset 1
#  To use this restraint we have to first define the data format
#  Here assuming that it's a CSV file with column names that may need to change
#  Other options include the linker length and the slope (for nudging components together)
columnmap={}
columnmap["Protein1"]="pep1.accession"
columnmap["Protein2"]="pep2.accession"
columnmap["Residue1"]="pep1.xlinked_aa"
columnmap["Residue2"]="pep2.xlinked_aa"
columnmap["IDScore"]=None

xl1 = IMP.pmi.restraints.crosslinking.ISDCrossLinkMS(representation,
                                   datadirectory+'polii_xlinks.csv',
                                   length=21.0,
                                   slope=0.01,
                                   columnmapping=columnmap,
                                   resolution=1.0,
                                   label="Trnka",
                                   csvfile=True)

xl1.add_to_model()             # crosslink must be added to the model

sampleobjects.append(xl1) #crosslink restraint is storing a sampled particle
outputobjects.append(xl1)


# Crosslinks - dataset 2
#  We can easily add a second set of crosslinks.
#  These have a different format and label, but other settings are the same
columnmap={}
columnmap["Protein1"]="prot1"
columnmap["Protein2"]="prot2"
columnmap["Residue1"]="res1"
columnmap["Residue2"]="res2"
columnmap["IDScore"]=None

xl2 = IMP.pmi.restraints.crosslinking.ISDCrossLinkMS(representation,
                                   datadirectory+'polii_juri.csv',
                                   length=21.0,
                                   slope=0.01,
                                   columnmapping=columnmap,
                                   resolution=1.0,
                                   label="Chen",
                                   csvfile=True)
xl2.add_to_model()

sampleobjects.append(xl2)
outputobjects.append(xl2)


#### optimize a bit before adding the EM restraint
representation.optimize_floppy_bodies(10)
####

# Electron Microscopy Restraint
#  The GaussianEMRestraint uses a density overlap function to compare model to data
#   First the EM map is approximated with a Gaussian Mixture Model (done separately)
#   Second, the componets of the model are represented with Gaussians (forming the model GMM)
#   Other options: scale_to_target_mass ensures the total mass of model and map are identical
#                  slope: nudge model closer to map when far away
#                  weight: experimental, needed becaues the EM restraint is quasi-bayesian
em_components = bm.get_density_hierarchies([t.domain_name for t in domains])
gemt = IMP.pmi.restraints.em.GaussianEMRestraint(em_components,
                                                 target_gmm_file,
                                                 scale_target_to_mass=True,
                                                 slope=0.000001,
                                                 weight=100.0)
gemt.add_to_model()
outputobjects.append(gemt)

#--------------------------
# Monte-Carlo Sampling
#--------------------------
# This object defines all components to be sampled as well as the sampling protocol
mc1=IMP.pmi.macros.ReplicaExchange0(m,
                                    representation,
                                    monte_carlo_sample_objects=sampleobjects,
                                    output_objects=outputobjects,
                                    crosslink_restraints=[xl1,xl2],    # allows XLs to be drawn in the RMF files
                                    monte_carlo_temperature=1.0,
                                    simulated_annealing=True,
                                    simulated_annealing_minimum_temperature=1.0,
                                    simulated_annealing_maximum_temperature=2.5,
                                    simulated_annealing_minimum_temperature_nframes=200,
                                    simulated_annealing_maximum_temperature_nframes=20,
                                    replica_exchange_minimum_temperature=1.0,
                                    replica_exchange_maximum_temperature=2.5,
                                    number_of_best_scoring_models=100,
                                    monte_carlo_steps=num_mc_steps,
                                    number_of_frames=num_frames,
                                    global_output_directory="output")

# Start Sampling
mc1.execute_macro()
