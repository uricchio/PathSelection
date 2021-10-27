#
import sys
import numpy as np

g0 = sys.argv[1]
init0 = sys.argv[2]

samp = sys.argv[3]
rep = sys.argv[4]


typ = sys.argv[5]

#if g1 == "5.0" and (init1 == "200" or init1 == "50"):
#@    exit()

# first file
f = "/Users/uricchio/projects/elisa/data/analysis/MCMCinit"+init0+".rate"+g0+".mutType"+typ+".samp"+samp+"."+rep+".traj.param.data"
fh = open(f, 'r')

sig = {}

for line in fh:
    data = [float(x) for x in line.strip().split()]
    if np.sign(data[2]) == np.sign(data[3]):
        sig[int(data[0])] = data[1]

fh.close()

for site in sig:
    print (site, sig[site], samp, rep, typ)

