#!/bin/bash
## running macs2 with shell script...
# pip install macs2  # (do we need to repeatedly install?)

for file in ./sim_folder_$1_do_$2/*.bam; do
  name=$(basename $file) 
  echo "running..."
  macs2 callpeak -t ./sim_folder_$1_do_$2/${name} --gsize 300734518 -m 2 50 --keep-dup=all -f BAM --outdir ./sim_folder_$1_do_$2/${name}_macs_peak -n ${name}_macs_peak
  echo "done!"
done