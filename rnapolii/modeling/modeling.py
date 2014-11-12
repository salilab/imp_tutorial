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

import os

# setting up parameters

rbmaxtrans = 2.00
fbmaxtrans = 3.00
rbmaxrot=0.04
outputobjects = []
sampleobjects = []

# setting up topology

m = IMP.Model()
simo = IMP.pmi.representation.Representation(m,upperharmonic=True,disorderedlength=False)

datadirectory="../data/"

'''
       # compname  hier_name    color         fastafile              fastaid          pdbname      chain    resrange                                     read    "BEADS"ize rigid_body super_rigid_body emnum_components emtxtfilename  emmrcfilename chain of super rigid bodies
domains=[("Rpb1",  "Rpb1_1",       0.0,     datadirectory+"1WCM.fasta.txt", "1WCM:A|PDBID|CHAIN|SEQUENCE",  datadirectory+"1WCM_map_fitted.pdb", "A",  (1,1140),    None,          5,       0,                 [0,1,2],     -40,   datadirectory+"Rpb1_1.txt",  datadirectory+"Rpb1_1.mrc",   [0]),
         ("Rpb1",  "Rpb1_2",       0.0,     datadirectory+"1WCM.fasta.txt", "1WCM:A|PDBID|CHAIN|SEQUENCE",  datadirectory+"1WCM_map_fitted.pdb", "A",  (1141,1274), None,          5,       1,                 [0,1,2],     -40,   datadirectory+"Rpb1_2.txt",  datadirectory+"Rpb1_2.mrc",   [0]),
         ("Rpb1",  "Rpb1_3",       0.0,     datadirectory+"1WCM.fasta.txt", "1WCM:A|PDBID|CHAIN|SEQUENCE",  datadirectory+"1WCM_map_fitted.pdb", "A",  (1275,-1),   None,          5,       0,                 [0,1,2],     -40,   datadirectory+"Rpb1_3.txt",  datadirectory+"Rpb1_3.mrc",   [0]),
         ("Rpb2",  "Rpb2_1",       1.0,     datadirectory+"1WCM.fasta.txt", "1WCM:B|PDBID|CHAIN|SEQUENCE",  datadirectory+"1WCM_map_fitted.pdb", "B",  (1,1102),    None,          5,       2,                 [0,1,3],     -40,   datadirectory+"Rpb2_1.txt",  datadirectory+"Rpb2_1.mrc",   [1]),
         ("Rpb2",  "Rpb2_2",       1.0,     datadirectory+"1WCM.fasta.txt", "1WCM:B|PDBID|CHAIN|SEQUENCE",  datadirectory+"1WCM_map_fitted.pdb", "B",  (1103,-1),   None,          5,       0,                 [0,1,3],     -40,   datadirectory+"Rpb2_2.txt",  datadirectory+"Rpb2_2.mrc",   [1]),
         ("Rpb3",  "Rpb3",         0.2,     datadirectory+"1WCM.fasta.txt", "1WCM:C|PDBID|CHAIN|SEQUENCE",  datadirectory+"1WCM_map_fitted.pdb", "C",  (1,-1),      None,          5,       3,                 [0,4],       -40,     datadirectory+"Rpb3.txt",  datadirectory+"Rpb3.mrc",   [2]),
         ("Rpb4",  "Rpb4",         0.3,     datadirectory+"1WCM.fasta.txt", "1WCM:D|PDBID|CHAIN|SEQUENCE",  datadirectory+"1WCM_map_fitted.pdb", "D",  (1,-1),      True,          5,       4,                 [0,5],       -40,     datadirectory+"Rpb4.txt",  datadirectory+"Rpb4.mrc",   [3]),
         ("Rpb5",  "Rpb5",         0.4,     datadirectory+"1WCM.fasta.txt", "1WCM:E|PDBID|CHAIN|SEQUENCE",  datadirectory+"1WCM_map_fitted.pdb", "E",  (1,-1),      None,          5,       5,                 [0,6],       -40,     datadirectory+"Rpb5.txt",  datadirectory+"Rpb5.mrc",   [4]),
         ("Rpb6",  "Rpb6",         0.5,     datadirectory+"1WCM.fasta.txt", "1WCM:F|PDBID|CHAIN|SEQUENCE",  datadirectory+"1WCM_map_fitted.pdb", "F",  (1,-1),      None,          5,       6,                 [0,7],       -40,     datadirectory+"Rpb6.txt",  datadirectory+"Rpb6.mrc",   [5]),
         ("Rpb7",  "Rpb7",         0.6,     datadirectory+"1WCM.fasta.txt", "1WCM:G|PDBID|CHAIN|SEQUENCE",  datadirectory+"1WCM_map_fitted.pdb", "G",  (1,-1),      True,          5,       7,                 [0,8],       -40,     datadirectory+"Rpb7.txt",  datadirectory+"Rpb7.mrc",   [6]),
         ("Rpb8",  "Rpb8",         0.7,     datadirectory+"1WCM.fasta.txt", "1WCM:H|PDBID|CHAIN|SEQUENCE",  datadirectory+"1WCM_map_fitted.pdb", "H",  (1,-1),      None,          5,       8,                 [0,9],       -40,     datadirectory+"Rpb8.txt",  datadirectory+"Rpb8.mrc",   [7]),
         ("Rpb9",  "Rpb9",         0.8,     datadirectory+"1WCM.fasta.txt", "1WCM:I|PDBID|CHAIN|SEQUENCE",  datadirectory+"1WCM_map_fitted.pdb", "I",  (1,-1),      None,          5,       9,                 [0,10],      -40,    datadirectory+"Rpb9.txt",  datadirectory+"Rpb9.mrc",   [8]),
         ("Rpb10", "Rpb10",        0.9,     datadirectory+"1WCM.fasta.txt", "1WCM:J|PDBID|CHAIN|SEQUENCE",  datadirectory+"1WCM_map_fitted.pdb", "J",  (1,-1),      None,          5,       10,                [0,11],      -40,    datadirectory+"Rpb10.txt",  datadirectory+"Rpb10.mrc",   [9]),
         ("Rpb11", "Rpb11",        0.1,     datadirectory+"1WCM.fasta.txt", "1WCM:K|PDBID|CHAIN|SEQUENCE",  datadirectory+"1WCM_map_fitted.pdb", "K",  (1,-1),      None,          5,       11,                [0,12],      -40,    datadirectory+"Rpb11.txt",  datadirectory+"Rpb11.mrc",   [10]),
         ("Rpb12", "Rpb12",        0.35,    datadirectory+"1WCM.fasta.txt", "1WCM:L|PDBID|CHAIN|SEQUENCE",  datadirectory+"1WCM_map_fitted.pdb", "L",  (1,-1),      None,          5,       12,                [0,13],      -40,    datadirectory+"Rpb12.txt",  datadirectory+"Rpb12.mrc",   [12])]
'''


       # compname  hier_name    color         fastafile              fastaid          pdbname      chain    resrange                                     read    "BEADS"ize rigid_body super_rigid_body emnum_components emtxtfilename  emmrcfilename chain of super rigid bodies
domains=[("Rpb1",  "Rpb1_1",       0.0,     datadirectory+"1WCM.fasta.txt", "1WCM:A|PDBID|CHAIN|SEQUENCE",  datadirectory+"1WCM_map_fitted.pdb", "A",  (1,1140),    None,          10,       None,                 None,     -40,   datadirectory+"Rpb1_1.txt",  datadirectory+"Rpb1_1.mrc",   None),
         ("Rpb1",  "Rpb1_2",       0.0,     datadirectory+"1WCM.fasta.txt", "1WCM:A|PDBID|CHAIN|SEQUENCE",  datadirectory+"1WCM_map_fitted.pdb", "A",  (1141,1274), None,          10,       None,                 None,     -40,   datadirectory+"Rpb1_2.txt",  datadirectory+"Rpb1_2.mrc",   None),
         ("Rpb1",  "Rpb1_3",       0.0,     datadirectory+"1WCM.fasta.txt", "1WCM:A|PDBID|CHAIN|SEQUENCE",  datadirectory+"1WCM_map_fitted.pdb", "A",  (1275,-1),   None,          10,       None,                 None,     -40,   datadirectory+"Rpb1_3.txt",  datadirectory+"Rpb1_3.mrc",   None),
         ("Rpb2",  "Rpb2_1",       1.0,     datadirectory+"1WCM.fasta.txt", "1WCM:B|PDBID|CHAIN|SEQUENCE",  datadirectory+"1WCM_map_fitted.pdb", "B",  (1,1102),    None,          10,       None,                 None,     -40,   datadirectory+"Rpb2_1.txt",  datadirectory+"Rpb2_1.mrc",   None),
         ("Rpb2",  "Rpb2_2",       1.0,     datadirectory+"1WCM.fasta.txt", "1WCM:B|PDBID|CHAIN|SEQUENCE",  datadirectory+"1WCM_map_fitted.pdb", "B",  (1103,-1),   None,          10,       None,                 None,     -40,   datadirectory+"Rpb2_2.txt",  datadirectory+"Rpb2_2.mrc",   None),
         ("Rpb3",  "Rpb3",         0.2,     datadirectory+"1WCM.fasta.txt", "1WCM:C|PDBID|CHAIN|SEQUENCE",  datadirectory+"1WCM_map_fitted.pdb", "C",  (1,-1),      None,          10,       None,                 None,       -40,     datadirectory+"Rpb3.txt",  datadirectory+"Rpb3.mrc",   None),
         ("Rpb4",  "Rpb4",         0.3,     datadirectory+"1WCM.fasta.txt", "1WCM:D|PDBID|CHAIN|SEQUENCE",  datadirectory+"1WCM_map_fitted.pdb", "D",  (1,-1),      True,          10,       4,                 [0,5],       -40,     datadirectory+"Rpb4.txt",  datadirectory+"Rpb4.mrc",     [3]),
         ("Rpb5",  "Rpb5",         0.4,     datadirectory+"1WCM.fasta.txt", "1WCM:E|PDBID|CHAIN|SEQUENCE",  datadirectory+"1WCM_map_fitted.pdb", "E",  (1,-1),      None,          10,       None,                 None,       -40,     datadirectory+"Rpb5.txt",  datadirectory+"Rpb5.mrc",   None),
         ("Rpb6",  "Rpb6",         0.5,     datadirectory+"1WCM.fasta.txt", "1WCM:F|PDBID|CHAIN|SEQUENCE",  datadirectory+"1WCM_map_fitted.pdb", "F",  (1,-1),      None,          10,       None,                 None,       -40,     datadirectory+"Rpb6.txt",  datadirectory+"Rpb6.mrc",   None),
         ("Rpb7",  "Rpb7",         0.6,     datadirectory+"1WCM.fasta.txt", "1WCM:G|PDBID|CHAIN|SEQUENCE",  datadirectory+"1WCM_map_fitted.pdb", "G",  (1,-1),      True,          10,       7,                 [0,8],       -40,     datadirectory+"Rpb7.txt",  datadirectory+"Rpb7.mrc",     [6]),
         ("Rpb8",  "Rpb8",         0.7,     datadirectory+"1WCM.fasta.txt", "1WCM:H|PDBID|CHAIN|SEQUENCE",  datadirectory+"1WCM_map_fitted.pdb", "H",  (1,-1),      None,          10,       None,                 None,       -40,     datadirectory+"Rpb8.txt",  datadirectory+"Rpb8.mrc",   None),
         ("Rpb9",  "Rpb9",         0.8,     datadirectory+"1WCM.fasta.txt", "1WCM:I|PDBID|CHAIN|SEQUENCE",  datadirectory+"1WCM_map_fitted.pdb", "I",  (1,-1),      None,          10,       None,                 None,      -40,    datadirectory+"Rpb9.txt",  datadirectory+"Rpb9.mrc",     None),
         ("Rpb10", "Rpb10",        0.9,     datadirectory+"1WCM.fasta.txt", "1WCM:J|PDBID|CHAIN|SEQUENCE",  datadirectory+"1WCM_map_fitted.pdb", "J",  (1,-1),      None,          10,       None,                 None,      -40,    datadirectory+"Rpb10.txt",  datadirectory+"Rpb10.mrc",   None),
         ("Rpb11", "Rpb11",        0.1,     datadirectory+"1WCM.fasta.txt", "1WCM:K|PDBID|CHAIN|SEQUENCE",  datadirectory+"1WCM_map_fitted.pdb", "K",  (1,-1),      None,          10,       None,                 None,      -40,    datadirectory+"Rpb11.txt",  datadirectory+"Rpb11.mrc",   None),
         ("Rpb12", "Rpb12",        0.35,    datadirectory+"1WCM.fasta.txt", "1WCM:L|PDBID|CHAIN|SEQUENCE",  datadirectory+"1WCM_map_fitted.pdb", "L",  (1,-1),      None,          10,       None,                 None,      -40,    datadirectory+"Rpb12.txt",  datadirectory+"Rpb12.mrc",   None)]





bm=IMP.pmi.macros.BuildModel1(simo)
bm.build_model(data_structure=domains)
bm.scale_bead_radii(40,0.8)

resdensities=bm.get_density_hierarchies([t[1] for t in domains])

# randomize the initial configuration

simo.shuffle_configuration(50)

# defines the movers

simo.set_rigid_bodies_max_rot(rbmaxrot)
simo.set_floppy_bodies_max_trans(fbmaxtrans)
simo.set_rigid_bodies_max_trans(rbmaxtrans)

outputobjects.append(simo)
sampleobjects.append(simo)

# scoring function

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






# sampling


simo.optimize_floppy_bodies(10)


mass=sum((IMP.atom.Mass(p).get_mass() for h in resdensities for p in IMP.atom.get_leaves(h)))
gemt = IMP.pmi.restraints.em.GaussianEMRestraint(resdensities,datadirectory+'emd_1883.map.mrc.gmm.50.txt',
                                               target_mass_scale=mass,
                                                slope=0.000001,
                                                target_radii_scale=3.0)
gemt.add_to_model()
gemt.set_weight(100.0)
outputobjects.append(gemt)

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
                                    number_of_frames=20000,
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


