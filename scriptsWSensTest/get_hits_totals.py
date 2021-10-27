
import sys
from collections import defaultdict

f = "g_init_N_table.txt"

fh = open(f, 'r')

N_g = defaultdict(dict)
numHits = defaultdict(dict)

for line in fh:
    data = line.strip().split()
    N_g[data[0]][data[1]] = data[2] 
    numHits[data[2]] = 0

fh.close()

f = "comp_num_hits.txt"

fh = open(f, 'r')

for line in fh:
    data = line.strip().split()
    numHits[N_g[data[1]][data[2]]] += float(data[0])
fh.close()

for N in numHits:
    print (N, numHits[N])


