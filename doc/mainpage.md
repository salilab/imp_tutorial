Modeling of complexes using IMP::pmi {#mainpage}
====================================

[TOC]

The _Python Modeling Interface_ (PMI) is a powerful set of tools designed
to handle all [steps of the modeling protocol](@ref procedure) for
typical modeling problems. It is designed to be used by writing a set of
Python scripts.

# Introduction {#introduction}

We will illustrate the use of PMI by determining the localization of two
subunits of RNA Polymerase II, utilizing chemical cross-linking coupled with
mass spectrometry, negative-stain electron microscopy (EM), and x-ray
crystallography data. We will try
to reconstruct the stalk of the complex, comprising of subunits Rpb4 and Rpb7,
hypothesizing that we know already the structure of the remaining 10-subunit
complex. The example can be easily generalized to any other set of subunits.

To work through the example on your own system, you will need the following
packages installed in addition to [IMP itself](@ref installation):

- [numpy and scipy](http://www.scipy.org/scipylib/download.html)
  for matrix and linear algebra

- [scikit-learn](http://scikit-learn.org/stable/install.html)
  for k-means clustering

- [matplotlib](http://matplotlib.org/downloads.html)
  for plotting results

- [Chimera](https://www.cgl.ucsf.edu/chimera/download.html)
  for visualization of results

(If you are using [Anaconda Python](https://www.anaconda.com/download/),
you can get the Python packages above by simply running
`conda install numpy scipy scikit-learn matplotlib`.
On a Mac you can get them using the
[pip](https://pypi.python.org/pypi/pip) tool, e.g. by running a command like
`sudo easy_install pip`, then install the packages with something like
`sudo pip install scikit-learn; sudo pip install matplotlib`. `numpy` and `scipy` are already installed on modern Macs. Something
similar may also work on a Linux box, although it's probably better to install
the packages using the distribution's package manager, such as `yum` or
`apt-get`.)

Then download the input files, either by 
[cloning the GitHub repository](https://github.com/salilab/imp_tutorial/tree/master)
or by [downloading the zip file](https://github.com/salilab/imp_tutorial/archive/master.zip).

The rnapolii example contains three directories: `analysis`, `data` and
`modeling`.

# Background of RNA Polymerase II {#background}

[RNA Pol II](http://en.wikipedia.org/wiki/RNA_polymerase_II) is a eukaryotic complex that catalyzes DNA transcription to synthesize mRNA strands.  Eukaryotic RNA polymerase II contains 12 subunits, Rpb1 to Rpb12. The yeast RNA Pol II dissociates into a 10-subunit core and a Rpb4/Rpb7 heterodimer. Rpb4 and Rpb7 are conserved from yeast to humans, and form a stalk-like protrusion extending from the main body of the RNA Pol II complex.


# Integrative Modeling using IMP {#usingimp}

This example will use data from chemical cross linking, EM and x-ray crystallography to localize the two subunits of the RNA Polymerase II stalk (Rpb4, Rpb7) to a static core of the remaining ten subunits.  

\image html rnapolii_integrative.png width=600px

# The four stages of Integrative Modeling {#stages}

Structural modeling using IMP is divided into [four stages](@ref procedure).

Click the links below to see a breakdown of all the modeling steps.

- \subpage gatherdata
  Collect biophysical data that can be used as structural restraints and constraints

- \subpage representation
  Define representations for the RNA Poly II structural model and define each data point as a scoring function.

- \subpage sampling
  Run a sampling protocol to find good scoring conformations.  

- \subpage analysis1 and \subpage analysis2
  Analysis of the good scoring conformations.  Clustering; uncertainty; precision; etc...


# Running the script {#script}

The first three modeling stages are all contained within one script, `modeling.py`. You can get started by simply changing into the `rnapolii/modeling` directory and then running the script with Python:

\code{.sh}
python modeling.py
\endcode

It will take a very long time to complete the sampling; to get an idea of what's going on you can run it with only 100 output frames by adding the `--test` option:

\code{.sh}
python modeling.py --test
\endcode

[On to stage 1...](@ref gatherdata)
