## Main Test

pushd .
cd ../sirfast_dev4
make clean
make
popd
clear

#../sirfast_dev4/sirfast 	--search ../reference/build37_chr1.fa  --seq ../reference/reads.tsv -e 0 -t 4	--min 250 --max 350 -u ../output/sing4.unmap -o ../output/sing4.out	> output_sing4
#../sirfast/sirfast			--search ../reference/build37_chr1.fa  --seq ../reference/reads.tsv -e 0		--min 250 --max 350 -u ../output/singR.unmap -o ../output/singR.out	#> output_singR

#../sirfast_dev4/sirfast 	--search ../reference/build37_chr1.fa --pe --seq ../reference/reads.tsv -e 0 -t 4	--min 250 --max 350 -u ../output/bug4.unmap -o ../output/bug4.out	#> output_bug4
#../sirfast/sirfast			--search ../reference/build37_chr1.fa --pe --seq ../reference/reads.tsv -e 0		--min 250 --max 350 -u ../output/bugR.unmap -o ../output/bugR.out	#> output_bugR

### 100K Single
../sirfast_dev4/sirfast 	--search ../reference/build37_chr1.fa --seq ../reference/reads100k.tsv -e 0 -t 200000	-u ../output/sing100k200k.unmap -o ../output/sing100k200k.out	> output_sing100k200k
#../sirfast_dev4/sirfast 	--search ../reference/build37_chr1.fa --seq ../reference/reads100k.tsv -e 0 -t 100000	-u ../output/sing100k100k.unmap -o ../output/sing100k100k.out	> output_sing100k100k
../sirfast_dev4/sirfast 	--search ../reference/build37_chr1.fa --seq ../reference/reads100k.tsv -e 0 -t 80000	-u ../output/sing100k080k.unmap -o ../output/sing100k080k.out	> output_sing100k080k
#../sirfast_dev4/sirfast 	--search ../reference/build37_chr1.fa --seq ../reference/reads100k.tsv -e 0 -t 40000	-u ../output/sing100k040k.unmap -o ../output/sing100k040k.out	> output_sing100k040k
../sirfast_dev4/sirfast 	--search ../reference/build37_chr1.fa --seq ../reference/reads100k.tsv -e 0 -t 20000	-u ../output/sing100k020k.unmap -o ../output/sing100k020k.out	> output_sing100k020k
../sirfast/sirfast			--search ../reference/build37_chr1.fa --seq ../reference/reads100k.tsv -e 0			-u ../output/sing100kRRRR.unmap -o ../output/sing100kRRRR.out	> output_sing100kRRRR

#sort sing100k100k.out > a1
#mv a1 sing100k100k.out
#sort sing100k080k.out > a2
#mv a2 sing100k080k.out
#sort sing100k040k.out > a3
#mv a3 sing100k040k.out
#sort sing100k020k.out > a4
#mv a4 sing100k020k.out
#sort sing100kRRRR.out > a5
#mv a5 sing100kRRRR.out

#diff sing100k100k.out sing100kRRRR.out >  diff_sing100k.out
#diff sing100k080k.out sing100kRRRR.out >> diff_sing100k.out
#diff sing100k040k.out sing100kRRRR.out >> diff_sing100k.out
#diff sing100k020k.out sing100kRRRR.out >> diff_sing100k.out

#diff sing100k100k.unmap sing100k100k.unmap >  diff_sing100k.unmap
#diff sing100k080k.unmap sing100k080k.unmap >> diff_sing100k.unmap
#diff sing100k040k.unmap sing100k040k.unmap >> diff_sing100k.unmap
#diff sing100k020k.unmap sing100k020k.unmap >> diff_sing100k.unmap


### 100K Pair
#../sirfast_dev4/sirfast 	--search ../reference/build37_chr1.fa --pe --seq ./reads100k.tsv -e 0 -t 200000	--min 250 --max 350 -u ../output/bug100k200k.unmap -o ../output/bug100k200k.out	> output_bug100k200k
#../sirfast_dev4/sirfast 	--search ../reference/build37_chr1.fa --pe --seq ./reads100k.tsv -e 0 -t 150000	--min 250 --max 350 -u ../output/bug100k150k.unmap -o ../output/bug100k150k.out	#> output_bug100k150k
#../sirfast_dev4/sirfast 	--search ../reference/build37_chr1.fa --pe --seq ./reads100k.tsv -e 0 -t 100000	--min 250 --max 350 -u ../output/bug100k100k.unmap -o ../output/bug100k100k.out	> output_bug100k100k
#../sirfast_dev4/sirfast 	--search ../reference/build37_chr1.fa --pe --seq ./reads100k.tsv -e 0 -t 80000	--min 250 --max 350 -u ../output/bug100k080k.unmap -o ../output/bug100k080k.out	#> output_bug100k080k
#../sirfast_dev4/sirfast 	--search ../reference/build37_chr1.fa --pe --seq ./reads100k.tsv -e 0 -t 40000	--min 250 --max 350 -u ../output/bug100k040k.unmap -o ../output/bug100k040k.out	#> output_bug100k040k
#../sirfast_dev4/sirfast 	--search ../reference/build37_chr1.fa --pe --seq ./reads100k.tsv -e 0 -t 20000	--min 250 --max 350 -u ../output/bug100k020k.unmap -o ../output/bug100k020k.out	> output_bug100k020k
#../sirfast/sirfast			--search ../reference/build37_chr1.fa --pe --seq ./reads100k.tsv -e 0			--min 250 --max 350 -u ../output/bug100kRRRR.unmap -o ../output/bug100kRRRR.out	> output_bug100kRRRR

#diff bug100k100k.out bug100kRRRR.out >  diff_pair100k.out
#diff bug100k080k.out bug100kRRRR.out >> diff_pair100k.out
#diff bug100k040k.out bug100kRRRR.out >> diff_pair100k.out
#diff bug100k020k.out bug100kRRRR.out >> diff_pair100k.out

#diff bug100k100k.out_OEA.sam bug100kRRRR.out_OEA.sam >	diff_pair100k.oea
#diff bug100k080k.out_OEA.sam bug100kRRRR.out_OEA.sam >> diff_pair100k.oea
#diff bug100k040k.out_OEA.sam bug100kRRRR.out_OEA.sam >> diff_pair100k.oea
#diff bug100k020k.out_OEA.sam bug100kRRRR.out_OEA.sam >> diff_pair100k.oea

#diff bug100k100k.unmap bug100k100k.unmap >  diff_pair100k.unmap
#diff bug100k080k.unmap bug100k080k.unmap >> diff_pair100k.unmap
#diff bug100k040k.unmap bug100k040k.unmap >> diff_pair100k.unmap
#diff bug100k020k.unmap bug100k020k.unmap >> diff_pair100k.unmap

