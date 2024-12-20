#!/bin/bash

## running macs2 peak-calling with shell script...

for file in ./sim_folder_$1_do_$2/*.bam; do
  name=$(basename $file) 
  echo "running..."
  macs2 callpeak -t ./sim_folder_$1_do_$2/${name} --gsize 300734518 -m 2 50 --keep-dup=all -f BAM --outdir ./sim_folder_$1_do_$2/${name}_macs_peak -n ${name}_macs_peak
  echo "done!"
done