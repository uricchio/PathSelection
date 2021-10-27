#
import sys
import numpy as np

g0 = sys.argv[1]
init0 = sys.argv[2]

samp = sys.argv[3]
rep = sys.argv[4]

g1 = sys.argv[5]
init1 = sys.argv[6]

typ = sys.argv[7]

if g1 == "5.0" and (init1 == "200" or init1 == "50"):
    exit()

# first file
f = "/Users/uricchio/projects/elisa/data/analysis/MCMCinit"+init0+".rate"+g0+".mutType"+typ+".samp"+samp+"."+rep+".traj.param.data"
fh = open(f, 'r')

sig = {}

for line in fh:
    data = [float(x) for x in line.strip().split()]
    if np.sign(data[2]) == np.sign(data[3]):
        sig[int(data[0])] = data[1]

fh.close()

f = "/Users/uricchio/projects/elisa/data/analysis/MCMCinit"+init1+".rate"+g1+".mutType"+typ+".samp"+samp+"."+rep+".traj.param.data"
fh = open(f, 'r')

sig1 = {}

for line in fh:
    data = [float(x) for x in line.strip().split()]
    if data[0] in sig:
        sig1[int(data[0])] = [data[1],data[4]]

for site in sig:
    if site in sig1:
        print (site, sig[site], sig1[site][0], samp, rep, sig1[site][1])

