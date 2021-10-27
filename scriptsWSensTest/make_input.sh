for init in {50,200,500}; do for rate in {1.2,2,5}; do python3 make_input_files_Const.py $init exp $rate INDEL ;  done; done
for init in {50,200,500}; do for rate in {1.2,2,5}; do python3 make_input_files_Const.py $init exp $rate SNP ;  done; done
for init in {50,200,500}; do for rate in {1.2,"2.0","5.0"}; do for mutType in {"SNP","INDEL"}; do mkdir ../MCMCout/MCMCinit${init}rate${rate}mutType$mutType; done; done; done
