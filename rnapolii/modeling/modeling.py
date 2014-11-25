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

######################
##  Header
######################


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
rbmaxtrans = 2.00
fbmaxtrans = 3.00
rbmaxrot = 0.04
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

bm.scale_bead_radii(40,0.8)

#----------------------------
# Define Degrees of Freedom
#----------------------------

# Add default mover parameters to simulation
representation.set_rigid_bodies_max_rot(rbmaxrot)
representation.set_floppy_bodies_max_trans(fbmaxtrans)
representation.set_rigid_bodies_max_trans(rbmaxtrans)

# These lists define what objects are sampled and what objects produce output
outputobjects = []
sampleobjects = []

# Add the movers to the sample and output object lists
outputobjects.append(representation)
sampleobjects.append(representation)

#--------------------------
# Define Scoring Function Components
#--------------------------

# Excluded Volume Restraint
ev = IMP.pmi.restraints.stereochemistry.ExcludedVolumeSphere(representation,resolution=10)

# Add to the model and append to output objects
ev.add_to_model()
outputobjects.append(ev)


#------------------
# Crosslink Restraints
#


#### Import Cross-link data:

# Define format of cross-link data files.
columnmap={}
columnmap["Protein1"]="pep1.accession"
columnmap["Protein2"]="pep2.accession"
columnmap["Residue1"]="pep1.xlinked_aa"
columnmap["Residue2"]="pep2.xlinked_aa"
columnmap["IDScore"]=None
columnmap["XLUniqueID"]=None

ids_map=IMP.pmi.tools.map()
ids_map.set_map_element(1.0,1.0)
xl1 = IMP.pmi.restraints.crosslinking.ISDCrossLinkMS(representation,
                                   datadirectory+'polii_xlinks.csv',
                                   length=21.0,             # Cross link length in angstroms
                                   slope=0.01,              
                                   columnmapping=columnmap,
                                   ids_map=ids_map,
                                   resolution=1.0,          # resolution at which restraint is evaluated (1=residue level)
                                   label="Trnka",   
                                   csvfile=True)


xl1.add_to_model()             # crosslink must be added to the model
xl1.set_psi_is_sampled(True)   # For bayesian modeling, we wish to sample uncertainty parameter
psi=xl1.get_psi(1.0)[0]        # create and set range for psi
psi.set_scale(0.05)

# Since we are sampling psi, crosslink restraint must be added to sampleobjects
sampleobjects.append(xl1)
outputobjects.append(xl1)


# crosslinks - dataset 2
columnmap={}
columnmap["Protein1"]="prot1"
columnmap["Protein2"]="prot2"
columnmap["Residue1"]="res1"
columnmap["Residue2"]="res2"
columnmap["IDScore"]=None
columnmap["XLUniqueID"]=None
ids_map=IMP.pmi.tools.map()
ids_map.set_map_element(1.0,1.0)

xl2 = IMP.pmi.restraints.crosslinking.ISDCrossLinkMS(representation,
                                   datadirectory+'polii_juri.csv',
                                   length=21.0,
                                   slope=0.01,
                                   columnmapping=columnmap,
                                   ids_map=ids_map,
                                   resolution=1.0,
                                   label="Chen",
                                   csvfile=True)
xl2.add_to_model()
xl2.set_label("Chen")
xl2.set_psi_is_sampled(True)
psi=xl2.get_psi(1.0)[0]
psi.set_scale(0.05)

sampleobjects.append(xl2)
outputobjects.append(xl2)



# optimize a bit before adding the EM restraint
representation.optimize_floppy_bodies(10)

#----------------------------
# Electron Microscopy Restraint
#----------------------------

# get total mass of system - this can be simplified, no?
#mass=bm.get_density_mass(...)?
resdensities=bm.get_density_hierarchies([t.domain_name for t in domains])
mass=sum((IMP.atom.Mass(p).get_mass() for h in resdensities for p in IMP.atom.get_leaves(h)))

gemt = IMP.pmi.restraints.em.GaussianEMRestraint(resdensities,
                                                 target_gmm_file,
                                                 target_mass_scale=mass,   
                                                 slope=0.000001,           # What is this slope?
                                                 target_radii_scale=3.0)   # What is this scale?
gemt.add_to_model()

# Weight of the EM restraint. -  This is the only restraint with a weight.  Why?
gemt.set_weight(100.0)


outputobjects.append(gemt)

#--------------------------
# Monte-Carlo Sampling
#--------------------------

# ***DS Why must the crosslink restraint be explicitly defined, but the other restraints do not?

# This object defines all components to be sampled as well as the sampling protocol
mc1=IMP.pmi.macros.ReplicaExchange0(m,
                                    representation,
                                    monte_carlo_sample_objects=sampleobjects,
                                    output_objects=outputobjects,
                                    crosslink_restraints=[xl1,xl2],
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
                                    do_clean_first=True,
                                    global_output_directory="output")

# Start Sampling
mc1.execute_macro()
