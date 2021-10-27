import numpy as np
import sys
import os

n = sys.argv[1]
init= sys.argv[2]
mutType= sys.argv[3]
rate = sys.argv[4]


def ret_same_sign(s0,s1,s2):
    
    myVars = []
    for var in s0:
        if np.sign(s0[var][1]) == np.sign(s0[var][2]) and np.sign(s1[var][1]) == np.sign(s1[var][2]) and np.sign(s2[var][1]) == np.sign(s2[var][2]):
            myVars.append(var)
    return(myVars)
 
def get_signs(t,init,rate,mutType):
    fh = open("../data/analysis/MCMCinit"+init+".rate"+rate+".mutType"+mutType+".samp"+t+".traj.param.data","r") 
    
    signs = {}
    for line in fh:
        data = [x for x in line.strip().split()]
        signs[data[0]] = [float(x) for x in data[1:]]
    fh.close()
    return (signs)

s0 = get_signs(n+".1",init,rate,mutType)
s1 = get_signs(n+".2",init,rate,mutType)
s2 = get_signs(n+".3",init,rate,mutType)

svals = [s0,s1,s2]

myVars = ret_same_sign(s0,s1,s2)
for var in myVars:
    for s in svals:
        sys.stdout.write(var+" ")
        for datum in s[var]:
            sys.stdout.write(str(datum)+" ")
        sys.stdout.write("\n")
    



