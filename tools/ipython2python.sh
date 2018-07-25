#!/bin/bash
jupyter nbconvert --to script ../doc/Tutorial.ipynb
grep -v "get_ipython" ../doc/Tutorial.py > tmp.py
sed "/# ### Using \`StatHierarchyHandler\`/q" tmp.py > modeling.py

echo 'from __future__ import print_function' > short_analysis.py
echo 'import IMP' >> short_analysis.py
echo 'import IMP.pmi' >> short_analysis.py
echo 'import IMP.pmi.output' >> short_analysis.py
echo 'import matplotlib.pyplot as plt' >> short_analysis.py
echo 'plot = plt.plot' >> short_analysis.py
echo 'figure = plt.figure' >> short_analysis.py
echo 'boxplot = plt.boxplot' >> short_analysis.py
echo 'savefig= plt.savefig' >> short_analysis.py

echo 'm=IMP.Model()' >> short_analysis.py
echo 'import warnings' >> short_analysis.py
echo "warnings.filterwarnings('ignore')"  >> short_analysis.py


sed -n "/# ### Using \`StatHierarchyHandler\`/,/# ## Stage 4 - Analysis/p" tmp.py >> short_analysis.py

echo 'from __future__ import print_function' > long_analysis.py
echo 'import IMP' >> long_analysis.py
echo 'import IMP.pmi' >> long_analysis.py
echo 'import IMP.pmi.macros' >> long_analysis.py
echo 'import matplotlib.pyplot as plt' >> long_analysis.py
echo 'plot = plt.plot' >> long_analysis.py 
echo 'figure = plt.figure' >> long_analysis.py
echo 'boxplot = plt.boxplot' >> long_analysis.py
echo 'savefig= plt.savefig' >> long_analysis.py

echo 'm=IMP.Model()' >> long_analysis.py
echo 'import warnings' >> long_analysis.py
echo "warnings.filterwarnings('ignore')"  >> long_analysis.py



sed -n "/# ## Stage 4 - Analysis/,//p" tmp.py >> long_analysis.py

jupyter nbconvert --to script ../doc/IMP_example.ipynb

mv ../doc/IMP_example.py ../modeling/
mv modeling.py ../modeling/
mv short_analysis.py ../modeling/
mv long_analysis.py ../modeling/

