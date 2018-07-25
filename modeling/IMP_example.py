
# coding: utf-8

# # IMP example

# ##What is IMP?
# 
# IMP is a C++/Python library. 
# 
# IMP provides building blocks and tools to allow methods developers 
# 
# * to convert data from new experimental methods into spatial restraints, 
# * to implement optimization and analysis techniques, 
# * to implement an integrative modeling procedure from scratch; 
# 
# the developer can use the C++ and Python programming languages to achieve these tasks.

# ## IMP implementation
# 
# IMP is primarily implemented in C++ for speed; 
# 
# each of the classes is wrapped so that it can also be used from Python. 

# IMP is organized into module. Each module contains classes, methods and data which are related. 

# In[ ]:


from __future__ import print_function
import IMP
import IMP.atom
import IMP.core


# ## The Model
# 
# In IMP, the system is represented by the `IMP.Model` class, which stores a collection of "particles".

# In[ ]:


# create an IMP model
m=IMP.Model()


# ## The Particle
# 
# A `IMP.Particle` is a flexible and abstract data container, able to hold whatever information is necessary to represent the system.

# In[ ]:


# create a new particle
pa=IMP.Particle(m)

# set the name
pa.set_name("My Particle A")


# ## Decoration
# 
# Decorators allows to access and assing to data associated to particles.

# In[ ]:


# decorate it as a sphere
dr=IMP.core.XYZR.setup_particle(pa)

# set the coordinates
dr.set_coordinates((0,0,0))

# set the radius
dr.set_radius(1.0)

# set the mass
IMP.atom.Mass.setup_particle(pa,1.0)

# set the optimization of the coordinates to True
dr.set_coordinates_are_optimized(True)


# ## Hierarchy
# 
# Biological modules are represented hierarchically 

# In[ ]:


# create a hierarchy
ha=IMP.atom.Hierarchy(pa)


# In[ ]:


# set color
ca=IMP.display.Color(1,1,0)
c=IMP.display.Colored.setup_particle(pa,ca)


# In[ ]:


# create a second particle
pb=IMP.Particle(m)
pb.set_name("My Particle B")
dr=IMP.core.XYZR.setup_particle(pb)
dr.set_coordinates((0,0,0))
dr.set_radius(1.0)
IMP.atom.Mass.setup_particle(pb,1.0)
dr.set_coordinates_are_optimized(True)
hb=IMP.atom.Hierarchy(pb)
cb=IMP.display.Color(0,1,1)
c=IMP.display.Colored.setup_particle(pb,cb)


# ## The Movers
# 
# Move particle attributes with specific rules

# In[ ]:


#now create the movers
mva=IMP.core.BallMover(m,pa,1.0)
mvb=IMP.core.BallMover(m,pb,1.0)


# ## Restraints
# 
# Every Restraint in IMP is implemented as a function that returns a score for some subset of the Model. 
# 
# A Restraint is satisfied by modifying the system to minimize its score.

# In[ ]:


# create an harmonic restraint
hf = IMP.core.Harmonic(4.0,1.0)
dr=IMP.core.DistanceRestraint(m,hf,pa,pb)


# ## Containers
# 
# Behind the scenes, IMP maintains an IMP::DependencyGraph that tracks how information flows between the particles and the containers, based on the constraints. 
# 
# It is used to optimize the efficiency of the computation

# In[ ]:


# and a restraint that contains the two particles
center = IMP.algebra.Vector3D(0,0,0)
ub = IMP.core.HarmonicUpperBound(20.0,1.0)
ss = IMP.core.DistanceToSingletonScore(ub, center)
lsc = IMP.container.ListSingletonContainer(m)
lsc.add([pa,pb])
rest = IMP.container.SingletonsRestraint(ss, lsc)


# ## Scoring Function 
# 
# This is used to calculate how well a configuration in the Model satisfies the Restraints 

# In[ ]:


# wrap the restraints in a Scoring Function
sf = IMP.core.RestraintsScoringFunction([rest,dr])


# ## Monte Carlo Sampling

# In[ ]:


# Build the Monte Carlo Sampler
mc = IMP.core.MonteCarlo(m)
mc.set_scoring_function(sf)
sm = IMP.core.SerialMover([mva,mvb])
mc.add_mover(sm)
mc.set_return_best(False)
mc.set_kt(1.0)


# ## RMF trajectories
# 
# Configurations of the model can be saved and visualized using RMF (Rich Molecular Format) files.

# In[ ]:


# Prepare the trajectory file
import IMP.rmf
import RMF

rh = RMF.create_rmf_file("out.rmf")
IMP.rmf.add_hierarchies(rh, [ha,hb])
IMP.rmf.add_restraints(rh,[dr])
IMP.rmf.save_frame(rh)


# ## Run the sampling

# In[ ]:


# run the sampling
for i in range(1000):
    mc.optimize(1)
    IMP.rmf.save_frame(rh)
    print(sf.evaluate(False))
    
del rh

