#!/bin/bash

# Build Hashtable
./sirfast --index ../reference/build37.fa --ws 10

# Single End Mode with first read format (-f 1), 1k read slice (-t 1000)
./sirfast --search ../reference/build37.fa      --seq ../reference/read_form1 -e 0 -f 1 -t 1000                     -o out_se_form1.sam -u out_se_form1.unmap

# Mate Pair Mode with first read format (-f 1), 1k read slice (-t 1000)
./sirfast --search ../reference/build37.fa --mp --seq ../reference/read_form1 -e 0 -f 1 -t 1000 --min 200 --max 600 -o out_pe_form1.sam -u out_pe_form1.unmap

# Single End Mode with second read format (-f 1), 1k read slice (-t 1000)
./sirfast --search ../reference/build37.fa      --seq ../reference/read_form2 -e 0 -f 2 -t 1000                     -o out_se_form2.sam -u out_se_form2.unmap

# Mate Pair Mode with second read format (-f 1), 1k read slice (-t 1000)
./sirfast --search ../reference/build37.fa --mp --seq ../reference/read_form2 -e 0 -f 2 -t 1000 --min 200 --max 600 -o out_pe_form2.sam -u out_pe_form2.unmap
