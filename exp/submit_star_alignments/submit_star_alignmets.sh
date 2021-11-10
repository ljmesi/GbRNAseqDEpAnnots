#!/bin/bash -l

module load bioinfo-tools gnuparallel/20180822
cat cuttrim2_SF-2243-GB-12-13_S47_L002.txt | \
parallel \
--interactive \
--verbose \
--jobs 1 \
--max-replace-args 1 \
"sbatch STARmap.sh {}"
