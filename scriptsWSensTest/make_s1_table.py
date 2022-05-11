import sys
from collections import defaultdict

fh = open(sys.argv[1],'r')

sig_pos_ind = defaultdict(dict)
sig_pos_snp = defaultdict(dict)

for line in fh:
    data = line.strip().split()
    if data[4] == "INDEL":  
        if data[0] not in sig_pos_ind or data[2] not in sig_pos_ind[data[0]]:
           sig_pos_ind[data[0]][data[2]] = []
        sig_pos_ind[data[0]][data[2]].append([data[1],data[3],data[4]])
    if data[4] == "SNP":  
        if data[0] not in sig_pos_snp or data[2] not in sig_pos_snp[data[0]]:
           sig_pos_snp[data[0]][data[2]] = []
        sig_pos_snp[data[0]][data[2]].append([data[1],data[3],data[4]])

for pos in sig_pos_ind:
    for exp in sig_pos_ind[pos]:
        if len(sig_pos_ind[pos][exp]) > 0:
            pr_str = pos+"\tINDEL\t"
            num = len(sig_pos_ind[pos][exp]) 
            for thing in sig_pos_ind[pos][exp]:
                pr_str += str(round(float(thing[0]),1))+","
            pr_str = pr_str[:-1]
            pr_str += "\t"+exp+"\t"+str(num)
            print(pr_str)

for pos in sig_pos_snp:
    for exp in sig_pos_snp[pos]:
        if len(sig_pos_snp[pos][exp]) > 0:
            pr_str = pos+"\tSNP\t"
            num = len(sig_pos_snp[pos][exp])
            for thing in sig_pos_snp[pos][exp]:
                pr_str += str(round(float(thing[0]),1))+","
            pr_str = pr_str[:-1]
            pr_str += "\t"+exp+"\t"+str(num)
            print(pr_str)



