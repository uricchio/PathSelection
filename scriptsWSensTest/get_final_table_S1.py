
import sys
import numpy as np
from collections import defaultdict


def tree(): return defaultdict(tree)

f = sys.argv[1]

fh = open(f,'r')

var = tree()

for line in fh:
    data = line.strip().split()
    var[data[0]][data[2]][data[4]] = []
fh.close()

fh = open(f,'r')

for line in fh:
    data = line.strip().split()
    var[data[0]][data[2]][data[4]].append([data[1],data[3]])

for pos in var:
    for exp in var[pos]:
        for typ in var[pos][exp]:
            if len(var[pos][exp][typ]) > 1:
                printStr = pos+","+typ+","
                for thing in var[pos][exp][typ]:
                    printStr += str(round(float(thing[0]),1))+" "
                printStr += " ,"+exp+" ,"
                signs = 0
                for thing in var[pos][exp][typ]:
                    printStr += thing[1]+" "
                    signs += np.sign(float(thing[0]))
                
                if signs == 0:
                    continue
                print(printStr) 


   
