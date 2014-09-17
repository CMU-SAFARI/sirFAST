#./generator (Path to  reference genome) (# of errors) (# of reads to generate) (fractions of reads having error [%])
#./generator /home/donghyu1/Projects/DNA/SirFAST/reference/build36.fa 0 10 50 > ../../reference/reads_100_e0
#./generator /home/donghyu1/Projects/DNA/SirFAST/reference/build36.fa 1 10 50 > ../../reference/reads_100_e1
#./generator /home/donghyu1/Projects/DNA/SirFAST/reference/build36.fa 2 10 50 > ../../reference/reads_100_e2

#./generator /home/donghyu1/Projects/DNA/SirFAST/reference/build37.fa 0 1000 50 > ../../reference/reads_1k_e0_verify
#./generator /home/donghyu1/Projects/DNA/SirFAST/reference/build37.fa 1 1000 50 > ../../reference/reads_1k_e1_verify
#./generator /home/donghyu1/Projects/DNA/SirFAST/reference/build37.fa 2 1000 50 > ../../reference/reads_1k_e2_verify

#./generator /home/donghyu1/Projects/DNA/SirFAST/reference/build37.fa 1 1000 50 > ../../reference/reads_1k_e1
#./generator /home/donghyu1/Projects/DNA/SirFAST/reference/build37.fa 2 1000 50 > ../../reference/reads_1k_e
#./generator ../../reference/build37.fa 1 10 100 > ../../reference/reads_1k_e1
#./generator ../../reference/build37.fa 0 100 100 > ../../reference/reads_100_e0

#./generator ../../reference/build37.fa 0 1000 50 > ../../reference/reads_1k_e0
./generator ../../reference/build37.fa 1 1000 50 > ../../reference/reads_1k_e1
./generator ../../reference/build37.fa 2 1000 50 > ../../reference/reads_1k_e2
