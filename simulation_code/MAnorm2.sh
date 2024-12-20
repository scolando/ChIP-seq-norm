#!/bin/bash

#uses bedtools package which not available through pip but is in homebrew
#uses MAnorm2_utils package which is available through pip


for file in ./sim_folder_$1_do_$2/MAnorm2/*.bam; do
    echo "converting bam reads to bed files"
    bedtools bamtobed -i ${file} > ${file}.bed
    echo "finished!"
done

echo "done!"

n=$3
combinedpeaks=""

for i in $(seq 1 $(($n > 0? $n: 0))); do
  token="./sim_folder_$1_do_$2/MAnorm2/tf_out_$i.narrowPeak"
  combinedpeaks="${combinedpeaks}${combinedpeaks:+,}$token"
done

combinedreads=""

for j in $(seq 1 $(($n > 0? $n: 0))); do
  token="./sim_folder_$1_do_$2/MAnorm2/tf_out_$j.bam.bed"
  combinedreads="${combinedreads}${combinedreads:+,}$token"
done

labs=""

for k in $(seq 1 $(($n > 0? $n/2: 0))); do
  token="s1"
  labs="${labs}${labs:+,}$token"
done
echo $labs

for k in $(seq 1 $(($n > 0? $n/2: 0))); do
  token="s2"
  labs="${labs}${labs:+,}$token"
done

echo $labs

echo "creating bin profile"
profile_bins --peaks=$combinedpeaks --typical-bin-size=1000 --reads=$combinedreads --labs=$labs -n ./sim_folder_$1_do_$2/MAnorm2/cond_A_and_B
echo "finished profiling"

### all conditions in one profile_bin

## reads and peaks are for each sample (so replicates but not between-samples)
## specify single end reads
## they recommend bin size of 1000 for TFs
## keep duplicates all is set by default