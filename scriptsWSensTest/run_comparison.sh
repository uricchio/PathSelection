rm comp_across_reps.txt
rm comp_num_hits.txt
rm table_s1_data.txt

#rm g_init_N_table.txt
#for N in {35,50,200}; do for r in {1.2,2.0,5.0}; do echo -n $r; echo -n " "; echo -n $N; echo -n " "; cat ../data/pop.hist.exp${r}.init${N}.txt; done; done >> g_init_N_table.txt

for typ in {"INDEL","SNP"}; do for i in {1,2,3}; do for exp in {2,9,17}; do for rate in {1.2,2.0,5.0}; do for init in {35,50,200}; do python compare_sig.py 5.0 35 $exp $i $rate $init $typ; done; done; done; done; done >> comp_across_reps.txt

for typ in {"INDEL","SNP"}; do for i in {1,2,3}; do for exp in {2,9,17}; do for rate in {1.2,2.0,5.0}; do for init in {35,50,200}; do myout=$(cat /Users/uricchio/projects/elisa/data/analysis/MCMCinit${init}.rate${rate}.mutType${typ}.samp${exp}.${i}.traj.param.data | awk '{if ($5 != 0) print}' | wc -l); echo $myout $rate $init $typ; done; done; done; done; done >> comp_num_hits.txt

python get_hits_totals.py > total_hits.txt

for var in {"INDEL","SNP"}; do for exp in {2,9,17}; do for i in {1,2,3}; do python compare_sig_table.py 5.0 35 $exp $i $var; done; done; done >> table_s1_data.txt

python make_s1_table.py table_s1_data.txt > table_S1_data_formatted.tsv
