import os
import numpy as np
import sys
from collections import defaultdict

samps = {}

# open data file
f = sys.argv[1]
fh = open(f, 'r')

for line in fh:
    break

for line in fh:
    data = line.strip().split()
    samps[data[0]] = [data[1],data[2]]


# open SNP frequencies
f = sys.argv[2]
fh = open(f, 'r')

col = {}
freqs = defaultdict(dict)
for line in fh:
    fields = line.strip().split(',')
    i = 0
    for field in fields:
        col[field] = i
        i += 1
    break

for line in fh:
    data = line.strip().split(',')
    if data[col["FREQ"]] == ".":
        data[col["FREQ"]] = "-1%"
    if data[col["POS"]] not in freqs:
        freqs[data[col["POS"]]][samps[data[col["SAMPLE"]]][1]] = [[float(data[col["FREQ"]][:-1]),samps[data[col["SAMPLE"]]][0]]]
    elif samps[data[col["SAMPLE"]]][1] not in freqs[data[col["POS"]]]:
        freqs[data[col["POS"]]][samps[data[col["SAMPLE"]]][1]] = [[float(data[col["FREQ"]][:-1]),samps[data[col["SAMPLE"]]][0]]]
    else:
        freqs[data[col["POS"]]][samps[data[col["SAMPLE"]]][1]].append([float(data[col["FREQ"]][:-1]),samps[data[col["SAMPLE"]]][0]])   

highfreqs = defaultdict(dict)
for pos in freqs:
    addkey = 0
    #if freqs[pos]["0"][0][0]/100. > 0.:
    #    continue
    for samp in freqs[pos]:
        for freq in freqs[pos][samp]:
            if freq[0]/100. > 0.3 and freq[0]/100. < 0.95:
                addkey += 1 
                break
        
    if addkey > 0:
        for samp in freqs[pos]:
            highfreqs[pos][samp] = freqs[pos][samp]
 
for pos in highfreqs:  
    for samp in highfreqs[pos]:
        sys.stdout.write(pos+' '+samp+',')
        highfreqs[pos][samp].sort(key=lambda x: x[1])
        for freq in highfreqs[pos][samp]:
            sys.stdout.write(str(freq[0])+',')
        #sys.stdout.write(' ')
        sys.stdout.write('\n')
