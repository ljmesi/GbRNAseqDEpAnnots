#!/bin/bash -l

module load bioinfo-tools gnuparallel/20180822

cat sample_names.txt | \
parallel \
--interactive \
--verbose \
--jobs 1 \
--max-replace-args 1 \
"sbatch STARmap.sh {}"
