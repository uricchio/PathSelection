rm comp_across_reps.txt
rm comp_num_hits.txt

for typ in {"INDEL","SNP"}; do for i in {1,2,3}; do for exp in {2,9,17}; do for rate in {1.2,2.0,5.0}; do for init in {50,200,500}; do python compare_sig.py 5.0 50 $exp $i $rate $init $typ; done; done; done; done; done >> comp_across_reps.txt

for typ in {"INDEL","SNP"}; do for i in {1,2,3}; do for exp in {2,9,17}; do for rate in {1.2,2.0,5.0}; do for init in {50,200,500}; do myout=$(cat /Users/uricchio/projects/elisa/data/analysis/MCMCinit${init}.rate${rate}.mutType${typ}.samp${exp}.${i}.traj.param.data | awk '{if ($5 != 0) print}' | wc -l); echo $myout $rate $init $typ; done; done; done; done; done >> comp_num_hits.txt

python get_hits_totals.py > total_hits.txt

for var in {"INDEL","SNP"}; do for exp in {2,9,17}; do for i in {1,2,3}; do python compare_sig_table.py 2.0 50 $exp $i $var; done; done; done >> table_s1_data.txt
