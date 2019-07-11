from __future__ import print_function
import IMP
import IMP.pmi
import IMP.pmi.output
import matplotlib.pyplot as plt
plot = plt.plot
figure = plt.figure
boxplot = plt.boxplot
savefig= plt.savefig
m=IMP.Model()
import warnings
warnings.filterwarnings('ignore')
# ### Using `StatHierarchyHandler` for inline analysis <a name="ProcessOutput_3"></a>
# 
# We can use the class StatHierarchyHandler to analyse and plot the content of the rmf files.
# First, we print all the keywords. This class coordinates the structures that have been generated 
# and all the features that have been saved during the modeling run. It is a Hierarchy object, but it works like a list.
# The python script can be found in `modeling/short_analysis.py`.

# In[ ]:


import IMP.pmi.output

hh=IMP.pmi.output.StatHierarchyHandler(m,"./output/rmfs/0.rmf3")

print("number of frames",len(hh))

print("describe the content of the stat file", hh[1])

#list down all the feature names
for k in hh[1].features.keys(): print(k)
    


# We can use the class IMP.atom.Selection to analyse the structures generated. 

# In[ ]:


# For instance we can compute the distance between two residues


p0=IMP.atom.Selection(hh,molecule="Rpb4",residue_index=10).get_selected_particles()[0]
p1=IMP.atom.Selection(hh,molecule="Rpb7",residue_index=10).get_selected_particles()[0]

d0=IMP.core.XYZ(p0)
d1=IMP.core.XYZ(p1)

plot([IMP.core.get_distance(d0,d1) for h in hh]);

figure()

# Or we can get the radius of gyration of the whole complex

ps=IMP.atom.Selection(hh).get_selected_particles()
plot([IMP.atom.get_radius_of_gyration(ps) for h in hh])


# Next, we plot the time series of selected features stored in the rmf file

# In[ ]:


#first we store the data internal to hh, so that it is now read from the files
# and it is faster

data=hh.data

# then we plot the scores
plot([x.score for x in data])

figure() 

# finally we plot distances of two crosslinks
plot([float(x.features["CrossLinkingMassSpectrometryRestraint_Data_Score_Chen"]) for x in data]);
plot([float(x.features["CrossLinkingMassSpectrometryRestraint_Distance_|Trnka|103.1|Rpb1|1|Rpb1|343|0|PSI|"]) for x in data]);



# Additionally, we can draw the box-plot for all crosslink distances involving Rpb7 and Rpb4

# In[ ]:


values=[]
for k in data[0].features.keys():
    if "Distance" in k and ("Rpb7" in  k or "Rpb4" in k):
        print(k)
        values.append([float(x.features[k]) for x in data])

boxplot(values);


# ## Stage 4 - Analysis <a name="Analysis_3"></a>
