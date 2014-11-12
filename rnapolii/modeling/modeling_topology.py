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

#---------------------------
# Set up Input Files
#---------------------------
datadirectory="../data/"
topology_file=datadirectory+"topology.txt"

#--------------------------
# Set MC Sampling Parameters
#--------------------------
num_steps=20000

#--------------------------
# Set Mover Parameters
#--------------------------
rbmaxtrans = 2.00
fbmaxtrans = 3.00
rbmaxrot=0.04


# ***DS It would be super nice if the user would not need to change anything below this line

###################################
# Here is where the work begins:
###################################

# ***DS Can we make these lists part Representation?
outputobjects = []
sampleobjects = []

#--------------------------------
# Define Components and Build Model
#--------------------------------

# Set up model and representation.  ***DS This should be condensed to one line or, better yet, included with an inputted topology file
m = IMP.Model()
simo = IMP.pmi.representation.Representation(m,upperharmonic=True,disorderedlength=False)
bm=IMP.pmi.macros.BuildModel1(simo)

# Create list of components from topology file
topology=IMP.pmi.topology.topology_io.TopologyReader(topology_file)
print type(topology)
domains=topology.component_list

# Import list of RB and list of SRB from file
rigid_bodies=None
super_rigid_bodies=None
chain_of_super_rigid_bodies=None

# Build model from components
bm.build_model(domains,
                     list_of_rigid_bodies=rigid_bodies, 
               list_of_super_rigid_bodies=super_rigid_bodies,
              chain_of_super_rigid_bodies=chain_of_super_rigid_bodies)


# ***DS Can we internalize these commands?
bm.scale_bead_radii(40,0.8)
resdensities=bm.get_density_hierarchies([t.domain_name for t in domains])

#----------------------------
# Define Degrees of Freedom
#----------------------------

# Randomize the initial configuration before sampling
# ***DS Can we add this as a parameter to the sampling call?  Like: randomize=True, steps=50
simo.shuffle_configuration(50)

# Add default mover parameters to simulation
simo.set_rigid_bodies_max_rot(rbmaxrot)
simo.set_floppy_bodies_max_trans(fbmaxtrans)
simo.set_rigid_bodies_max_trans(rbmaxtrans)

outputobjects.append(simo)
sampleobjects.append(simo)

#--------------------------
# Define and Import Scoring Function Components
#--------------------------

# Excluded Volume Restraint
ev = IMP.pmi.restraints.stereochemistry.ExcludedVolumeSphere(simo,resolution=10)
ev.add_to_model()
outputobjects.append(ev)


# here we have a protocol for clean datesets, 
# when the model is flexible and there are not large rigid parts


columnmap={}
columnmap["Protein1"]="pep1.accession"
columnmap["Protein2"]="pep2.accession"
columnmap["Residue1"]="pep1.xlinked_aa"
columnmap["Residue2"]="pep2.xlinked_aa"
columnmap["IDScore"]=None
columnmap["XLUniqueID"]=None

ids_map=IMP.pmi.tools.map()
ids_map.set_map_element(1.0,1.0)

xl1 = IMP.pmi.restraints.crosslinking.ISDCrossLinkMS(simo,
                                   datadirectory+'polii_xlinks.csv',
                                   length=21.0,
                                   slope=0.01,
                                   columnmapping=columnmap,
                                   ids_map=ids_map,
                                   resolution=1.0,
                                   label="Trnka",
                                   csvfile=True)
xl1.add_to_model()
xl1.set_label("Trnka")
sampleobjects.append(xl1)
outputobjects.append(xl1)
xl1.set_psi_is_sampled(True)
psi=xl1.get_psi(1.0)[0]
psi.set_scale(0.05)

columnmap={}
columnmap["Protein1"]="prot1"
columnmap["Protein2"]="prot2"
columnmap["Residue1"]="res1"
columnmap["Residue2"]="res2"
columnmap["IDScore"]=None
columnmap["XLUniqueID"]=None

ids_map=IMP.pmi.tools.map()
ids_map.set_map_element(1.0,1.0)

xl2 = IMP.pmi.restraints.crosslinking.ISDCrossLinkMS(simo,
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
sampleobjects.append(xl2)
outputobjects.append(xl2)
xl2.set_psi_is_sampled(True)
psi=xl2.get_psi(1.0)[0]
psi.set_scale(0.05)


simo.optimize_floppy_bodies(10)


mass=sum((IMP.atom.Mass(p).get_mass() for h in resdensities for p in IMP.atom.get_leaves(h)))
gemt = IMP.pmi.restraints.em.GaussianEMRestraint(resdensities,datadirectory+'emd_1883.map.mrc.gmm.50.txt',
                                               target_mass_scale=mass,
                                                slope=0.000001,
                                                target_radii_scale=3.0)
gemt.add_to_model()
gemt.set_weight(100.0)
outputobjects.append(gemt)

#--------------------------
# Set up Sampling Macro
#--------------------------

# ***DS What values here should be "default-ier" than others...i.e. defined in the header of this script?

mc1=IMP.pmi.macros.ReplicaExchange0(m,
                                    simo,
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
                                    number_of_best_scoring_models=0,
                                    monte_carlo_steps=10,
                                    number_of_frames=num_steps,
                                    write_initial_rmf=True,
                                    initial_rmf_name_suffix="initial",
                                    stat_file_name_suffix="stat",
                                    best_pdb_name_suffix="model",
                                    do_clean_first=True,
                                    do_create_directories=True,
                                    global_output_directory="output",
                                    rmf_dir="rmfs/",
                                    best_pdb_dir="pdbs/",
                                    replica_stat_file_suffix="stat_replica")
mc1.execute_macro()


