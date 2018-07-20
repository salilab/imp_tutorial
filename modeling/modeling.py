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

import warnings
warnings.filterwarnings('ignore')


# Then setup the relevant paths of the input files
datadirectory = "../rnapolii/data/"
topology_file = datadirectory+"topology.txt" 
target_gmm_file = datadirectory+'emd_1883.map.mrc.gmm.50.txt' # The EM map data
output_directory = "./output"

# This part of the script defines the topology of the system, including the hierarchy, the representation and the degrees of freedom. This is the content of the file `../data/topology.txt`, which is in a table format:

'''
|molecule_name  |color     |fasta_fn          |fasta_id|pdb_fn             |chain|residue_range|pdb_offset|bead_size|em_residues_per_gaussian|rigid_body|super_rigid_body|chain_of_super_rigid_bodies|
|Rpb1           |blue      |1WCM_new.fasta.txt|1WCM:A  |1WCM_map_fitted.pdb|A    |1,1140       |0         |20       |40                      |1         | 1              |                           |
|Rpb1           |blue      |1WCM_new.fasta.txt|1WCM:A  |1WCM_map_fitted.pdb|A    |1141,1274    |0         |20       |40                      |2         | 1              |                           |
|Rpb1           |blue      |1WCM_new.fasta.txt|1WCM:A  |1WCM_map_fitted.pdb|A    |1275,END     |0         |20       |40                      |3         | 1              |                           |
|Rpb2           |red       |1WCM_new.fasta.txt|1WCM:B  |1WCM_map_fitted.pdb|B    |1,1102       |0         |20       |40                      |4         | 2              |                           |
|Rpb2           |red       |1WCM_new.fasta.txt|1WCM:B  |1WCM_map_fitted.pdb|B    |1103,END     |0         |20       |40                      |5         | 2              |                           |
|Rpb3           |yellow    |1WCM_new.fasta.txt|1WCM:C  |1WCM_map_fitted.pdb|C    |1,END        |0         |20       |40                      |6         | 3              |                           |
|Rpb4           |salmon    |1WCM_new.fasta.txt|1WCM:D  |1WCM_map_fitted.pdb|D    |1,END        |0         |20       |40                      |7         | 4              |                           |
|Rpb5           |gold      |1WCM_new.fasta.txt|1WCM:E  |1WCM_map_fitted.pdb|E    |1,END        |0         |20       |40                      |8         | 5              |                           |
|Rpb6           |pink      |1WCM_new.fasta.txt|1WCM:F  |1WCM_map_fitted.pdb|F    |1,END        |0         |20       |40                      |9         | 6              |                           |
|Rpb7           |gray      |1WCM_new.fasta.txt|1WCM:G  |1WCM_map_fitted.pdb|G    |1,END        |0         |20       |40                      |10        | 7              |                           |
|Rpb8           |orange    |1WCM_new.fasta.txt|1WCM:H  |1WCM_map_fitted.pdb|H    |1,END        |0         |20       |40                      |11        | 8              |                           |
|Rpb9           |tan       |1WCM_new.fasta.txt|1WCM:I  |1WCM_map_fitted.pdb|I    |1,END        |0         |20       |40                      |12        | 9              |                           |
|Rpb10          |brown     |1WCM_new.fasta.txt|1WCM:J  |1WCM_map_fitted.pdb|J    |1,END        |0         |20       |40                      |13        | 10             |                           |
|Rpb11          |purple    |1WCM_new.fasta.txt|1WCM:K  |1WCM_map_fitted.pdb|K    |1,END        |0         |20       |40                      |14        | 11             |                           |
|Rpb12          |cyan      |1WCM_new.fasta.txt|1WCM:L  |1WCM_map_fitted.pdb|L    |1,END        |0         |20       |40                      |15        | 12             |                           |
''';

# Using the table above we define the overall topology: we introduce the molecules with their sequence and their known structure, and define the movers. Each line is a user-defined molecular **Domain**, and each column contains the specifics needed to build the system.
# 
# * `component_name`: Name of the Molecule and the name of the Hierarchy that contains the corresponding Domain.* `color`: The color used in the output coordinates file. Uses Chimera names.
# * `fasta_fn`: Name of the FASTA file containing the sequence for this Molecule.
# * `fasta_id`: header line of FASTA (without the ">" character.
# * `pdb_fn`: Name of PDB file with coordinates (if available). If left empty, will set up as BEADS (you can also specify "BEADS"). Can also write "IDEAL_HELIX".
# * `chain`: Chain ID of this Domain in the PDB file.
# * `residue_range`: Comma delimited pair defining the indexes of the first and the last residue of the Domain. Can leave empty or use 'all' for entire sequence from PDB file. The second item in the pair can be 'END' to select the last residue in the sequence defined in the FASTA file.
# * `pdb_offset`: Offset to sync PDB residue numbering with FASTA numbering.
# * `bead_size`: The size (in residues) of beads used to model Fragments not covered by PDB coordinates.
# * `em_residues`: The number of Gaussians used to model the  density of this domain. Can be set to zero/empty to exclude Domains from the em restraint. The GMM files will be written to `gmm_dir`
# * `rigid_body`: Leave empty if the user does not desire to build a Rigid Body on this Domain. Otherwise, use a unique Rigid Body identifier (an integer). All Domains with the same Rigid Body identifier will be collected in the same Rigid Body.
# * `super_rigid_body`: Like the Rigid Body, the user can specify a unique Super Rigid Body identifier.
# * `chain_of_super_rigid_bodies` Like the Rigid Body, the user can specify a unique Chain of Super Rigid Body identifier.
# 
# The first section defines where input files are located.  The topology file defines how the system components are structurally represented. `target_gmm_file` stores the EM map for the entire complex, which has already been converted into a Gaussian mixture model.


# Initialize IMP model
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


# Here we can set the **Degrees of Freedom** parameters, which should be
# optimized according to MC acceptance ratios. There are three kind of movers: Rigid Body, Bead, and Super Rigid Body. 
# 
# `max_rb_trans` and `max_rb_rot` are the 
# maximum translation and rotation of the Rigid Body mover, `max_srb_trans` and 
# `max_srb_rot` are the maximum translation and rotation of the Super Rigid Body mover
# and `max_bead_trans` is the maximum translation of the Bead Mover.
# 
# The excecution of the macro will return the root hierarchy (`root_hier`) and the degrees of freedom (`dof`) objects, both of which are used later on.
# 

root_hier, dof = bs.execute_macro(max_rb_trans=4.0, 
                                  max_rb_rot=0.3, 
                                  max_bead_trans=4.0, 
                                  max_srb_trans=4.0,
                                  max_srb_rot=0.3)


# At this point we have created the complete representation of the system. If displayed using Chimera the subunits should look like this (where the left and right panels are the assembled complex exploded view, respectively)

# We can display the representation of system the along the sequence. Each color correspond to a domain of the complex assigned to an individual rigid body. White spaces are the beads.

import IMP.pmi.plotting
import IMP.pmi.plotting.topology

IMP.pmi.plotting.topology.draw_component_composition(dof)


# Since we're interested in modelling the stalk, we will fix all subunits except Rpb4 and Rpb7. 

# Fix all rigid bodies but not Rpb4 and Rpb7 (the stalk)
# First select and gather all particles to fix.
fixed_particles=[]
for prot in ["Rpb1","Rpb2","Rpb3","Rpb5","Rpb6","Rpb8","Rpb9","Rpb10","Rpb11","Rpb12"]:
    fixed_particles+=IMP.atom.Selection(root_hier,molecule=prot).get_selected_particles()
    

# Fix the Corresponding Rigid movers and Super Rigid Body movers using dof
# The flexible beads will still be flexible (fixed_beads is an empty list)!
fixed_beads,fixed_rbs=dof.disable_movers(fixed_particles,
                                         [IMP.core.RigidBodyMover,IMP.pmi.TransformMover])


# Finally we randomize the initial configuration to remove any bias from the initial starting configuration read from input files. Since each subunit is composed of rigid bodies (i.e., beads constrained in a structure) and flexible beads, the configuration of the system is initialized by displacing each mobile rigid body and each bead randomly by 50 Angstroms, and rotate them randomly, and far enough from each other to prevent any steric clashes. 
# 
# The `excluded_rigid_bodies=fixed_rbs` will exclude from the randomization everything that was fixed above.

# Shuffle the rigid body and beads configuration of only the molecules we are interested in (Rpb4 and Rpb7)
IMP.pmi.tools.shuffle_configuration(root_hier,
                                    excluded_rigid_bodies=fixed_rbs,
                                    max_translation=50, 
                                    verbose=False,
                                    cutoff=5.0,
                                    niterations=100)


# After defining the representation of the model, we build the **restraints** by which the individual structural models will be scored based on the input data.
# 
# The sum of all of these restraints is our **scoring function**. 
# For all restraints, calling `add_to_model()` incorporates them into the scoring function
# Appending the restraints to the outputobjects list reports them in the log files produced in the sampling.

outputobjects = [] # reporter objects...output is included in the stat file


# Connectivity keeps things connected along the backbone (ignores if inside same rigid body)
mols = IMP.pmi.tools.get_molecules(root_hier)
for mol in mols:
    molname=mol.get_name()        
    IMP.pmi.tools.display_bonds(mol)
    cr = IMP.pmi.restraints.stereochemistry.ConnectivityRestraint(mol,scale=2.0)
    cr.add_to_model()
    cr.set_label(molname)
    outputobjects.append(cr)

ev = IMP.pmi.restraints.stereochemistry.ExcludedVolumeSphere(
                                         included_objects=root_hier,
                                         resolution=10)
ev.add_to_model()         # add to scoring function
outputobjects.append(ev)  # add to output
 
# A crosslinking restraint is implemented as a distance restraint between two residues.  The two residues are each defined by the protein (component) name and the residue number.  The script here extracts the correct four columns that provide this information from the input data file.
# 
# To use this restraint we have to first define the data format.  
#  
# This data file has a csv format, and looks like:
# 
# ```
# prot1,res1,prot2,res2
# Rpb1,34,Rpb1,49
# Rpb1,101,Rpb1,143
# ...
# ```


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
                                   slope=0.01,           # A linear term that biases XLed
                                                         # residues together
                                   resolution=1.0,       # Resolution at which to apply the restraint. 
                                                         # Either 1 (residue) or 0 (atomic)
                                   label="Chen",         # Used to label output in the stat file
                                   weight=1.)            # Weight applied to all crosslinks 
                                                         # in this dataset
xl1rest.add_to_model()
outputobjects.append(xl1rest)

# We can easily add a second set of crosslinks.
# These have a different format and label, but other settings are the same

# In[ ]:


xldbkwc = IMP.pmi.io.crosslink.CrossLinkDataBaseKeywordsConverter()
xldbkwc.set_protein1_key("pep1.accession")
xldbkwc.set_protein2_key("pep2.accession")
xldbkwc.set_residue1_key("pep1.xlinked_aa")
xldbkwc.set_residue2_key("pep2.xlinked_aa")

xl2 = IMP.pmi.io.crosslink.CrossLinkDataBase(xldbkwc)
xl2.create_set_from_file(datadirectory+'polii_xlinks.csv')

xl2rest = IMP.pmi.restraints.crosslinking.CrossLinkingMassSpectrometryRestraint(
                                   root_hier=root_hier,
                                   CrossLinkDataBase=xl2,
                                   length=21.0,
                                   slope=0.01,
                                   resolution=1.0,
                                   label="Trnka",
                                   weight=1.)

xl2rest.add_to_model()          
outputobjects.append(xl2rest)

# The GaussianEMRestraint uses a density overlap function to compare model to data
# First the EM map is approximated with a Gaussian Mixture Model (done separately)
# Second, the components of the model are represented with Gaussians (forming the model GMM)
# 
# The GaussianEMRestraint uses a density overlap function to compare model to data. First the EM map is approximated with a Gaussian Mixture Model (done separately). Second, the components of the model are represented with Gaussians (forming the model GMM)
# 
# * `scale_target_to_mass` ensures the total mass of model and map are identical
# * `slope`: nudge model closer to map when far away
# * `weight`: heuristic, needed to calibrate the EM restraint with the other terms. 
# 
# and then add it to the output object.  Nothing is being sampled, so it does not need to be added to sample objects.

# First, get the model density objects that will be fitted to the EM density.
em_components = IMP.pmi.tools.get_densities(root_hier)

gemt = IMP.pmi.restraints.em.GaussianEMRestraint(em_components,
                                      target_gmm_file,  # EM map GMM file
                                      scale_target_to_mass=True,  # True if the mass of the map and model are identical.
                                      slope=0.0000001,  # A small funneling force pulling towards the center of the EM density.
                                      weight=80.0)           
gemt.add_to_model()
outputobjects.append(gemt)


# With the system representation built and data restraints entered, the system is now ready to sample configurations. A replica exchange run can be set up using the `ReplicaExchange0` macro:
# 
# See the [ReplicaExchange0 documentation](https://integrativemodeling.org/nightly/doc/ref/classIMP_1_1pmi_1_1macros_1_1ReplicaExchange0.html) for a full description of all of the input parameters.
# 
# The sampling is performed by executing the macro:
# 
# ```mc1.execute_macro()```
# 


# total number of saved frames
num_frames = 50

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
              number_of_best_scoring_models=10,
              monte_carlo_steps=10,
              number_of_frames=num_frames,
              global_output_directory=output_directory)

# Start Sampling
mc1.execute_macro()

# The script generates an output directory containing the following:
# 
# * pdbs: a directory containing the 100 best-scoring models (see the number_of_best_scoring_models variable above) from the run, in PDB format.
# * rmfs: a single RMF file containing all the frames. RMF is a file format specially designed to store coarse-grained, multi-resolution and multi-state models such as those generated by %IMP. It is a compact binary format and (as in this case) can also be used to store multiple models or trajectories.
# * Statistics from the sampling, contained in a "statfile", stat.*.out. This file contains information on each restraint, MC acceptance criteria and other things at each step.

