
# coding: utf-8

# # IMP.pmi Tutorial Handout

# ### Integrating EM and Crosslinking data to localize two domains of RNA Polymerase II
# 
# Authors: Riccardo Pellarin, Max Bonomi, Charles Greenberg, Daniel Saltzberg, Ben Webb
# 
# Institut Pasteur, CNRS, C3BI
# USCF, Department of bioengineering and therapeutic sciences

# The Python Modeling Interface (pmi) is a powerful set of tools designed
# to handle all steps of the modeling protocol for
# typical modeling problems. It is designed to be used by writing a set of
# Python scripts.

# IMP.pmi has been used to determine the architecture of several macromolecular complexes, for instance:
# 
# [26S-PIP](https://salilab.org/26S-PIPs), [Yeast 40S-eIF3](https://salilab.org/40S-eIF1-eIF3), [Human Complement](https://salilab.org/Complement), [exosome](https://salilab.org/exosome),
#     [yeast mediator](https://salilab.org/mediator/), [Nup84](https://salilab.org/nup84), [TFIIH](https://salilab.org/tfiih), [Nup82](https://salilab.org/nup82/), [SEA complex](https://salilab.org/sea), and the [Nuclear Pore Complex](https://salilab.org/npc2018)
#     
# Each repository above contains the scripts and the data, as well as all the results, that are needed to reproduce the published results. 
# 
# For a given system, integrative modeling files are stored in different servers:
# 
# - source code is stored in a github repository (e.g., https://github.com/integrativemodeling/npc2018) 
# - data files are stored in the Zenodo data server (e.g., https://zenodo.org/record/1194547#.W02gVq3v5UQ)
# - structures are stored in the pdb-dev server (e.g., https://pdb-dev.wwpdb.org/)

# We will illustrate the use of IMP.pmi by determining the localization of two
# subunits of RNA Polymerase II, utilizing chemical cross-linking coupled with
# mass spectrometry, negative-stain electron microscopy (EM), and x-ray
# crystallography data of the subunits. We will try
# to reconstruct the stalk of the complex, comprising of subunits Rpb4 and Rpb7,
# hypothesizing that we know already the structure of the remaining 10-subunit
# complex. The example can be easily generalized to any other set of subunits.
# 
# IMP.pmi references
# 
# [Saltzberg et al. 2018](https://salilab.org/pdf/Saltzberg_MethodsMolBiol_2019.pdf)
# 
# [Bonomi et al. 2018](https://salilab.org/pdf/Bonomi_Structure_2018.pdf)
# 
# 
# 
# ## Installation
# 
# The current version of the Tutorial is guaranteed to work with IMP 2.11.0. This version can be installed un many plaforms using [anaconda](https://anaconda.org/salilab/imp), which provides all the dependencies.
# 
# To work through the example on your own system, you will need the following
# packages installed in addition to [IMP itself](https://integrativemodeling.org/nightly/doc/manual/installation.html):
# 
# - [numpy and scipy](http://www.scipy.org/scipylib/download.html)
#   for matrix and linear algebra
# 
# - [scikit-learn](http://scikit-learn.org/stable/install.html)
#   for k-means clustering
# 
# - [matplotlib](http://matplotlib.org/downloads.html)
#   for plotting results
# 
# - [Chimera](https://www.cgl.ucsf.edu/chimera/download.html)
#   for visualization of results
# 
# On a Mac you can get them using the
# [pip](https://pypi.python.org/pypi/pip) tool, e.g. by running a command like
# `sudo easy_install pip`, then install the packages with something like
# `sudo pip install scikit-learn; sudo pip install matplotlib`. `numpy` and `scipy` are already installed on modern Macs. Something
# similar may also work on a Linux box, although it's probably better to install
# the packages using the distribution's package manager, such as `yum` or
# `apt-get`.)
# 
# Then download the input files, either by 
# [cloning the GitHub repository](https://github.com/salilab/imp_tutorial/tree/develop)
# or by [downloading the zip file](https://github.com/salilab/imp_tutorial/archive/develop.zip).

# The rnapolii example scripts are contained in the directory `modeling`.

# ## Table of Content
# 
# [//]: # (To compile the Table of Content run `python tools/compile_toc.py Tutorial.ipynb` and paste the output here below)
# 
# &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[ Background of RNA Polymerase II ](#3_Background_of_RNA_Polymerase_II)
# 
# &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[ Integrative Modeling using IMP ](#4_Integrative_Modeling_using_IMP)
# 
# &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[ The four stages of Integrative Modeling ](#3_The_four_stages_of_Integrative_Modeling)
# 
# &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[ Running the script ](#3_Running_the_script)
# 
# &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[ Stage 1 - Gathering of data ](#Stage_1_2)
# 
# &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[ Data for yeast RNA Polymerase II ](#Data_rnapolii_3)
# 
# &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[ Stage 2 - Representation of subunits and translation of the data into spatial restraints ](#Stage_2_2)
# 
# &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[ Setting up Model Representation and Degrees of Freedom in IMP ](#Setting_up_3)
# 
# &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[ Hierarchy ](#Hierarchy_3)
# 
# &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[ Dissecting the script ](#Dissecting_the_script_3)
# 
# &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[ Model Representation Using a Topology File. ](#Topology_file_4)
# 
# &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[ Building the System Representation and Degrees of Freedom ](#Representation_and_DOF_4)
# 
# &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[ Scoring Function ](#Scoring_Function_3) 
# 
# &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[ Connectivity Restraint ](#Connectivity_Restraint_4) 
# 
# &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[ Excluded Volume Restraint ](#Excluded_Volume_Restraint_4) 
# 
# &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[ Crosslinks - dataset 1 ](#Crosslink_1_4)
# 
# &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[ Crosslinks - dataset 2 ](#Crosslink_2_4)
# 
# &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[ Electron Microscopy Restraint ](#EM_4) 
# 
# &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[ Stage 3 - Sampling ](#Sampling_2)
# 
# &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[ Modeling Output ](#Output_3)
# 
# &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[  Using `StatHierarchyHandler` for inline analysis ](#ProcessOutput_3)
# 
# &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[ Stage 4 - Analysis ](#Analysis_3)
# 
# &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[Clustering top models using `analysis.py` ](#Clustering_3)
# 
# &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[ Structural Uncertainty of the solutions ](#uncertainty_3)
# 
# &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[ Accuracy evaluation ](#Accuracy_3)
# 
# &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[ Sampling Exhaustiveness ](#Sampling_Exhaustiveness_3)
# 
# 

# ## Background of RNA Polymerase II <a name="3_Background_of_RNA_Polymerase_II"></a>
# 
# [RNA Pol II](http://en.wikipedia.org/wiki/RNA_polymerase_II) is a eukaryotic complex that catalyzes DNA transcription into mRNA strands.  Eukaryotic RNA polymerase II contains 12 subunits, Rpb1 to Rpb12. The yeast RNA Pol II dissociates into a 10-subunit core and a Rpb4/Rpb7 heterodimer. Rpb4 and Rpb7 are conserved from yeast to humans, and form a stalk-like protrusion extending from the main body of the RNA Pol II complex.
# 
# 
# ### Integrative Modeling using IMP <a name="4_Integrative_Modeling_using_IMP"></a>
# 
# This example will use data from chemical cross linking, EM and x-ray crystallography to localize the two subunits of the RNA Polymerase II stalk (Rpb4, Rpb7) to a static core of the remaining ten subunits.  
# 
# <img src="files/images/rnapolii_integrative.png" alt="Drawing" style="width: 600px;"/>
# 
# ### The four stages of Integrative Modeling <a name="3_The_four_stages_of_Integrative_Modeling"></a>
# 
# Structural modeling using IMP is divided into [four stages](@ref procedure).
# 
# Click the links below to see a breakdown of all the modeling steps.
# 
# 
# * [Stage 1](#Stage_1_2)) Collect biophysical data that can be used as structural restraints and constraints
#   
# * [Stage 2](#Stage_2_2)) Define representations for the RNA Poly II structural model and define each data point as a scoring function.
# 
# * [Stage 3](#Sampling_2)) Run a sampling protocol to find good scoring conformations.  
# 
# * [Stage 4](#Analysis_3)) Analysis of the good scoring conformations.  Clustering; uncertainty; precision; etc...
# 
# ### Running the script <a name="3_Running_the_script"></a>
# 
# The first three modeling stages are all contained within one script, `modeling.py`. You can get started by simply changing into the `rnapolii/modeling` directory and then running the script with Python:
# 
# ```
# python modeling.py
# ```
# 
# It will take a very long time to complete the sampling; to get an idea of what's going on you can run it with only 50 output frames by adding the `--test` option:
# 
# ```
# python modeling.py --test
# ```

# ## Stage 1 - Gathering of data <a name="Stage_1_2"></a>
# 
# In this stage, we find all available experimental data that we wish to utilize in structural modeling.  In theory, any method that provides information about absolute or relative structural information can be used.
# 
# ### Data for yeast RNA Polymerase II <a name="Data_rnapolii_3"></a>
# The `rnapolii/data` folder in the tutorial input files contains the data included in this example:
# 
# * Sequence information (FASTA files for each subunit)
# * [Electron density maps](http://www.ebi.ac.uk/pdbe/entry/EMD-1883/visualization) (`.mrc`, `.txt` files)
# * [High resolution structure from x-ray crystallography](http://www.rcsb.org/pdb/explore/explore.do?structureId=1WCM) (PDB file)
# * Chemical crosslinking datasets (we use two data sets, one from [Al Burlingame's lab](http://www.mcponline.org/content/13/2/420.long), and another from [Juri Rappsilber's lab](http://emboj.embopress.org/content/29/4/717))
# 
# 
# **FASTA File**  
# Each residue included in modeling must be explicitly defined in the FASTA text file.  Each individual component (i.e., a protein chain) is identified by a string in the FASTA header line.  From `1WCM.fasta.txt`:
# 
#     >1WCM:A
#     MVGQQYSSAPLRTVKEVQFGLFSPEEVRAISVAKIRFPETMDETQTRAKIGGLNDPRLGSIDRNLKCQTCQEGMNECPGH
#     FGHIDLAKPVFHVGFIAKIKKVCECVCMHCGKLLLDEHNELMRQALAIKDSKKRFAAIWTLCKTKMVCETDVPSEDDPTQ  
#     ...
#     >1WCM:B
#     MSDLANSEKYYDEDPYGFEDESAPITAEDSWAVISAFFREKGLVSQQLDSFNQFVDYTLQDIICEDSTLILEQLAQHTTE
#     SDNISRKYEISFGKIYVTKPMVNESDGVTHALYPQEARLRNLTYSSGLFVDVKKRTYEAIDVPGRELKYELIAEESEDDS  
#     ...
# 
# defines two chains with unique IDs of 1WCM:A and 1WCM:B respectively.  The entire complex is 12 chains and 4582 residues.
# 
# **Electron Density Map**  
# The electron density map of the entire RNA Poly II complex is at 20.9 Angstrom resolution.  The raw data file for this is stored in `emd_1883.map.mrc`.
# 
# <figure><img src="files/images/rnapolii_em_raw.png" width="200px" />
# <figcaption>_Electron microscopy density map for yeast RNA Polymerase II_</figcaption></figure>
# 
# **Electron Density as Gaussian Mixture Models**  
# Gaussian mixture models (GMMs) are used to greatly speed up scoring by approximating the electron density of individual subunits and experimental EM maps.  A GMM has been created for the experimental density map, and is stored in `emd_1883.map.mrc.gmm.50.mrc`.  The weight, center, and covariance matrix of each Gaussian used to approximate the original EM density can be seen in the corresponding `.txt` file.  
# 
# <figure><img src="files/images/rnapolii_em_gmm_50.png" width="200px" />
# <figcaption>_The EM data represented as a 50 Gaussian mixture model_</figcaption></figure>
# 
# 
# **PDB File**  
# High resolution coordinates for all 12 chains of RNA Pol II are found in `1WCM.pdb`.  
# 
# <figure><img src="files/images/rnapolii_all_1wc4.png" width="200px" />
# <figcaption>_Coordinates from PDBID [1WCM](http://www.rcsb.org/pdb/explore.do?structureId=1wcm)_</figcaption></figure>
# 
# **Chemical Cross-Links**  
# All chemical cross-linking data is located in `polii_xlinks.csv` and `polii_juri.csv`.  These files contain multiple comma-separated columns; four of these specify the protein and residue number for each of the two linker residues.
# 
#     prot1,res1,prot2,res2
#     Rpb1,34,Rpb1,49
#     Rpb1,101,Rpb1,143
#     Rpb1,101,Rpb1,176
# 
# The length of the DSS/BS3 cross-linker reagent, 21 angstroms, will be specified later in the modeling script.  

# ## Stage 2 - Representation of subunits and translation of the data into spatial restraints <a name="Stage_2_2"></a>
# 
# 
# In this stage, we will initially define a representation of the system. Afterwards, we will convert the data into spatial restraints.  This is performed using the script `modeling/modeling.py` and uses the
# topology file, `topology.txt`, to define the system components and their representation
# parameters.
# 
# ### Setting up Model Representation and Degrees of Freedom in IMP <a name="Setting_up_3"></a>
# 
# Very generally, the *representation* of a system is defined by all the variables that need to be determined based on input information, including the assignment of the system components to geometric objects (e.g. points, spheres, ellipsoids, and 3D Gaussian density functions). 
# 
# Our RNA Pol II representation employs **spherical beads** of varying sizes and **3D Gaussians**, which coarsen domains of the complex using several resolution scales simultaneously. 
# 
# <figure><img src="files/images/rnapolii_Multi-scale_representation.png" width="600px" />
# <figcaption>_Multi-scale representation of Rpb1 subunit of RNA Pol II_</figcaption></figure>
# 
# The **spatial restraints** will be applied to individual resolution scales as appropriate. 
# 
# Beads and Gaussians of a given domain are arranged into either a rigid body or a flexible string, based on the crystallographic structures. 
# 
# The GMM of a subunit is the set of all 3D Gaussians used to represent it; it will be used to calculate the EM score. The calculation of the GMM of a subunit can be done automatically in the **topology file**.
# For the purposes of this tutorial, we already created these for Rpb4 and Rpb7 and placed them in the `rnapolii/data` directory in their respective `.mrc` and `.txt` files. 
# 
# In a **rigid body**, all the beads and the Gaussians of a given domain have their relative distances constrained during configurational sampling, while in a **flexible string** the beads and the Gaussians are restrained by the sequence connectivity. 
# 
# 
# <figure><img src="files/images/rnapolii_rb.png" width="300px" />
# <figcaption>_Rigid Bodies and beads_</figcaption></figure>
# 
# **super rigid bodies** are sets of rigid bodies and beads that will move together in an additional Monte Carlo move.
# 
# <figure><img src="files/images/rnapolii_srb.png" width="300px" />
# <figcaption>_Super Rigid Bodies_</figcaption></figure>
# 
# **chain_of_super_rigid_bodies** are additional degrees of freedom along the connectivity chain of a subunit. It groups sequence-connected rigid domains and/or beads into overlapping pairs and triplets. Each of these groups will be moved rigidly. This mover helps to sample more efficiently complex topologies, made of several rigid bodies, connected by flexible linkers.
# 
# <figure><img src="files/images/rnapolii_cosrb.png" width="300px" />
# <figcaption>_Chain of Super Rigid Bodies_</figcaption></figure>
# 
# 
# ### Hierarchy <a name="Hierarchy_3"></a>
# 
# A hierarchy in IMP is a tree that stores information on molecules, residues, atoms, etc., where the resolution of the representation increases as you move further from the root. IMP.pmi was designed to support a specialised multi-state/multi-copy/multi-resolution hierarchy
# 
# <figure><img src="files/images/rnapolii_hierarchy.png" width="600px" />
# <figcaption>_PMI hierarchy_</figcaption></figure>
# 
# The **States** are used as putative structural and compositional alternatives of the system. 
# 
# Each **State** contains the **Molecules**, and each Molecule can occur in different stochiometric **Copies** (eg. here MolA has three identical copies: MolA.0, MolA.1, and MolA.2). 
# 
# The **Molecules** contains structures (ie, particles with coordinates, masses and radii) classified by several **resolutions**: Atomic (Resolution 0), Residues (Resolution 1), Fragments (Resolution > 1), and the Gaussians (Densities). 
# 
# All structures (except the densities) are represented by Spheres with appropriate radius and mass. The resolutions concur simultanously, therefore the same part of the molecule can be represented by several resolutions.

# ### Dissecting the script <a name="Dissecting_the_script_3"></a>
# 
# The script `modeling/modeling.py` sets up the representation of the system and the restraint.
# 
# The first part of the script import the necessary libraries.

# In[ ]:


from __future__ import print_function

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

# In[ ]:


datadirectory = "../rnapolii/data/"
topology_file = datadirectory+"topology.txt" 
target_gmm_file = datadirectory+'emd_1883.map.mrc.gmm.50.txt' # The EM map data
output_directory = "./output"


# #### Model Representation Using a Topology File. <a name="Topology_file_4"></a>
# 
# This part of the script defines the topology of the system, including the hierarchy, the representation and the degrees of freedom. This is the content of the file `../data/topology.txt`, which is in a table format:

# In[ ]:


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

# In[ ]:


# Initialize IMP model
m = IMP.Model()

# Read in the topology file.  
# Specify the directory wheere the PDB files, fasta files and GMM files are
topology = IMP.pmi.topology.TopologyReader(topology_file, 
                                  pdb_dir=datadirectory, 
                                  fasta_dir=datadirectory, 
                                  gmm_dir=datadirectory)


# In[ ]:


# Use the BuildSystem macro to build states from the topology file
bs = IMP.pmi.macros.BuildSystem(m)


# In[ ]:


# Each state can be specified by a topology file.
bs.add_state(topology)


# #### Building the System Representation and Degrees of Freedom <a name="Representation_and_DOF_4"></a>
# 
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

# In[ ]:


root_hier, dof = bs.execute_macro(max_rb_trans=4.0, 
                                  max_rb_rot=0.3, 
                                  max_bead_trans=4.0, 
                                  max_srb_trans=4.0,
                                  max_srb_rot=0.3)


# At this point we have created the complete representation of the system. If displayed using Chimera the subunits should look like this (where the left and right panels are the assembled complex exploded view, respectively)

# <figure><img src="files/images/rnapolii_domain_representation.png" width="600px" />
# <figcaption>_Domain Representation_</figcaption></figure>

# We can display the representation of system the along the sequence. Each color correspond to a domain of the complex assigned to an individual rigid body. White spaces are the beads.

# In[ ]:



import IMP.pmi.plotting
import IMP.pmi.plotting.topology

IMP.pmi.plotting.topology.draw_component_composition(dof)


# Since we're interested in modelling the stalk, we will fix all subunits except Rpb4 and Rpb7. Note that we are using [IMP.atom.Selection](https://integrativemodeling.org/nightly/doc/ref/classIMP_1_1atom_1_1Selection.html) to get the particles that correspond to the fixed Molecules.

# In[ ]:


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

# In[ ]:


# Shuffle the rigid body and beads configuration of only the molecules we are interested in (Rpb4 and Rpb7)
IMP.pmi.tools.shuffle_configuration(root_hier,
                                    excluded_rigid_bodies=fixed_rbs,
                                    max_translation=50, 
                                    verbose=False,
                                    cutoff=5.0,
                                    niterations=100)


# ### Scoring Function <a name="Scoring_Function_3"></a>

# After defining the representation of the model, we build the **restraints** by which the individual structural models will be scored based on the input data.
# 
# The sum of all of these restraints is our **scoring function**. 
# For all restraints, calling `add_to_model()` incorporates them into the scoring function
# Appending the restraints to the outputobjects list reports them in the log files produced in the sampling.

# In[ ]:


outputobjects = [] # reporter objects...output is included in the stat file


# #### Connectivity Restraint <a name="Connectivity_Restraint_4"></a>

# #### Excluded Volume Restraint <a name="Excluded_Volume_Restraint_4"></a>

# In[ ]:


# Connectivity keeps things connected along the backbone (ignores if inside same rigid body)
mols = IMP.pmi.tools.get_molecules(root_hier)
for mol in mols:
    molname=mol.get_name()        
    IMP.pmi.tools.display_bonds(mol)
    cr = IMP.pmi.restraints.stereochemistry.ConnectivityRestraint(mol,scale=2.0)
    cr.add_to_model()
    cr.set_label(molname)
    outputobjects.append(cr)


# In[ ]:


ev = IMP.pmi.restraints.stereochemistry.ExcludedVolumeSphere(
                                         included_objects=root_hier,
                                         resolution=10)
ev.add_to_model()         # add to scoring function
outputobjects.append(ev)  # add to output


# #### Crosslinks - dataset 1 <a name="Crosslink_1_4"></a>
# 
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

# In[ ]:


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


# #### Crosslinks - dataset 2 <a name="Crosslink_2_4"></a>
# 
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


# #### Electron Microscopy Restraint <a name="EM_4"></a>

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
# 

# In[ ]:


# First, get the model density objects that will be fitted to the EM density.
em_components = IMP.pmi.tools.get_densities(root_hier)

gemt = IMP.pmi.restraints.em.GaussianEMRestraint(em_components,
                                      target_gmm_file,  # EM map GMM file
                                      scale_target_to_mass=True,  # True if the mass of the map and model are identical.
                                      slope=0.0000001,  # A small funneling force pulling towards the center of the EM density.
                                      weight=80.0)           
gemt.add_to_model()
outputobjects.append(gemt)


# ## Stage 3 - Sampling <a name="Sampling_2"></a>
# 
# With the system representation built and data restraints entered, the system is now ready to sample configurations. A replica exchange run can be set up using the `ReplicaExchange0` macro:
# 
# See the [ReplicaExchange0 documentation](https://integrativemodeling.org/nightly/doc/ref/classIMP_1_1pmi_1_1macros_1_1ReplicaExchange0.html) for a full description of all of the input parameters.
# 
# The sampling is performed by executing the macro:
# 
# ```mc1.execute_macro()```
# 

# In[ ]:


# total number of saved frames
num_frames = 5000

# This object defines all components to be sampled as well as the sampling protocol
mc1=IMP.pmi.macros.ReplicaExchange0(m,
              root_hier=root_hier,                         # The root hierarchy
              monte_carlo_sample_objects=dof.get_movers(), # All moving particles and parameters
              rmf_output_objects=outputobjects,            # Objects to put into the rmf file
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


# ### Modeling Output <a name="Output_3"></a>
# 
# The script generates an output directory containing the following:
# 
# * pdbs: a directory containing the 100 best-scoring models (see the number_of_best_scoring_models variable above) from the run, in PDB format.
# * rmfs: a single RMF file containing all the frames. RMF is a file format specially designed to store coarse-grained, multi-resolution and multi-state models such as those generated by %IMP. It is a compact binary format and (as in this case) can also be used to store multiple models or trajectories. It stores the hierarchy and the coordinates of the particles, as well as information on each restraint, MC acceptance criteria and other things at each step.
# * Statistics from the sampling, contained in a "statfile", stat.*.out. 

# ### Using `StatHierarchyHandler` for inline analysis <a name="ProcessOutput_3"></a>
