import sys
import numpy as np
import gzip

f = sys.argv[1]
fields0 = f.split('/')
fields = fields0[3].split('.')
pos = int(fields[0][3:])

fh = gzip.open(f, 'r')
alphas = []

for line in fh:
    break
for line in fh:
    data = line.strip().split()
    alpha1 = float(data[3])
    alpha2 = float(data[4])

    alphas.append((alpha1+alpha2)/2.)

if len(alphas) == 0:
    print (pos,'nan','nan','nan','nan')
    exit()

si  = 0
if np.sign(np.percentile(alphas,5)) == np.sign(np.percentile(alphas,95)) and np.sign(np.percentile(alphas,5)) > 0:
   si = 1

if np.sign(np.percentile(alphas,5)) == np.sign(np.percentile(alphas,95)) and np.sign(np.percentile(alphas,5)) < 0:
   si =2
print (pos,np.mean(alphas),np.percentile(alphas,5),np.percentile(alphas,95),si)
