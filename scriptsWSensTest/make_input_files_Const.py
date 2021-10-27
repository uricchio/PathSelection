import sys
import os
import math
import numpy as np
from collections import defaultdict

N0 = sys.argv[1]
model = sys.argv[2]
rate = float(sys.argv[3])
mutType = sys.argv[4]

# make a "viral" pop history file assuming exp growth and final pop size
# of 10**10-ish particles

def get_harm_mean(f,Nf,rate):
  
    fh   = open(f,'w')

    t = 0
    N0 = 10**10
    N = N0

    Ns = [1./N]

    ngens = 0
    while(N >= Nf):
        N /= rate
        Ns.append(1./N)
        ngens += 1

    N = 1./np.mean(Ns)
    out = str(N)+"\t0\t"+str("-Inf")+"\n"
    fh.write(out)
    fh.close()

    return(N,ngens)

# read in the data from Elisa's experiment
def read_traj(f0,f1):
    samps = {}

    # open data file
    fh = open(f0, 'r')

    for line in fh:
        break

    for line in fh:
        data = line.strip().split()
        samps[data[0]] = [data[1],data[2]]
    fh.close()

    # open SNP frequencies
    fh = open(f1, 'r')

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
        #if data[col["FREQ"]] == ".":
        #    data[col["FREQ"]] = "-1%"
        #print (data[col["FREQ"]],data[col["AD"]],data[col["DP"]])
        if data[col["POS"]] not in freqs:
            freqs[data[col["POS"]]][samps[data[col["SAMPLE"]]][1]] = [["-1","-1","-1"] for x in range(0,5) ] #[[data[col["AD"]],data[col["RD"]],data[col["DP"]]]]
        elif samps[data[col["SAMPLE"]]][1] not in freqs[data[col["POS"]]]:
            freqs[data[col["POS"]]][samps[data[col["SAMPLE"]]][1]] =[["-1","-1","-1"] for x in range(0,5) ] # [[data[col["AD"]],data[col["RD"]],data[col["DP"]]]]
         
        if data[col["AD"]] == "." or data[col["RD"]] == ".":
            freqs[data[col["POS"]]][samps[data[col["SAMPLE"]]][1]][int(samps[data[col["SAMPLE"]]][0])] =[".",".","."]
        else:
            freqs[data[col["POS"]]][samps[data[col["SAMPLE"]]][1]][int(samps[data[col["SAMPLE"]]][0])] = [data[col["AD"]],data[col["RD"]],str(int(data[col["AD"]])+int(data[col["RD"]]))]
        # how come AD +RD is not always equal to DP or SDP? especially indels?

    for pos in freqs:
        for samp in freqs[pos]:
            if samp == "0":
                continue
            freqs[pos][samp][0] = freqs[pos]["0"][0]
    return (freqs) 

# print out trajectory files for all alleles
def print_traj(freqs,thresh,pos,sample,N0,ngens,Ninit):

    if freqs[pos][sample][0][0] == ".":
        return 

    def t_of_i(N0,rate,ngens):

        t0 = ngens/N0
        t = 0

        i = 9
        times = {}
        while i >=0:
            times[i] = t
            t += t0
            i -= 1

        return times

    # check to make all the calls exist so freq can be calculated
    def not_int(ar):
        for item in ar:
            try:
                i = int(item)
            except:
                return 1
        return 0
  
    times = t_of_i(N0,rate,ngens)
    try:
        os.stat("../data/traj"+model+str(rate)+".init"+str(Ninit)+".mutType"+mutType)
    except:
        os.mkdir("../data/traj"+model+str(rate)+".init"+str(Ninit)+".mutType"+mutType)
 

    for samp in freqs[pos]:
        if samp != sample:
            continue
        ts_vals = [0,1,4,6,9]
    
        # if there are no calls above 0 or only a small number of calls, skip
        # 
        badCall = 0
        noVar = 0
        lowCall = 0
        for t in freqs[pos][samp]:
            if t[0] == ".":
                badCall += 1
            if t[0] != "." and (t[0] == 0 or t[0] == t[2]):
                noVar += 1
            if t[0] != "." and (int(t[0]) < 2 or (int(t[2]) - int(t[0])) < 2):            
                lowCall +=1

        if noVar + badCall > 3:
            continue    

        if lowCall + noVar + badCall > 3 and mutType == "SNP":
            continue
        
        fh = open("../data/traj"+model+str(rate)+".init"+str(Ninit)+".mutType"+mutType+"/pos"+pos+".samp"+sample+".traj",'w')

        i =4
        for t in freqs[pos][samp][::-1]:
            ts_val = str(-1*times[ts_vals[i]])
            pr = t[0]
            if pr == "." or t[2] == ".":
                i -= 1
                continue
            if float(freqs[pos][samp][0][0])/float(freqs[pos][samp][0][2]) > 0.5:
                pr = int(t[2])-int(t[0])
            fh.write(str(pr)+' '+str(t[2])+' '+ts_val+' '+ts_val+'\n')
            i -= 1
        fh.close()

datafile = ''
if mutType == "INDEL":
    datafile ="../data/SpecializationIndelDP20.csv"
else:
    datafile = "../data/SpecializationSNPDP20.csv"

(harm_mean, ngens) = get_harm_mean("../data/pop.hist."+model+str(rate)+".init"+str(N0)+".txt",int(N0),rate)
freqs = read_traj("../data/sample_by_rep.txt",datafile) 

for pos in freqs:
    for sample in freqs[pos]:
        if sample == "0" or sample == 0:
            continue
        print_traj(freqs, 0.0, pos,sample,harm_mean,ngens,N0)
