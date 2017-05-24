
	    
	    
import re
import sys

search_term = "<a name="
f = sys.argv[1]

for line in open(f, 'r'):
    if re.search(search_term, line):
        #print line
        t=line.replace("#","&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[").replace(";[&",";&").replace('<a name=\\"',"](#").replace('\\"></a>\\n",',")").replace('"',' ').replace("\ ></a>",")").lstrip(' ')
	print t
