#MCMCinit200rate1.2mutTypeINDEL

rm ../data/analysis/MCMCinit*.rate*.mut*.samp*traj.param.data
for init in {50,200,500}; do for rate in {1.2,"2.0","5.0"}; do for mut in {"INDEL","SNP"}; do for samp in {2,9,17}; do for i in {1..3}; do for file in ../MCMCout/MCMCinit${init}rate${rate}mutType${mut}/pos*samp${samp}.$i.traj.param*; do python plot_sel.py $file >> ../data/analysis/MCMCinit${init}.rate${rate}.mutType${mut}.samp${samp}.$i.traj.param.data; done; done; done; done; done; done
