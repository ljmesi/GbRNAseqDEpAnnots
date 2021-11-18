#!/bin/bash -l

set -e
set -uo pipefail

module load conda
export CONDA_ENVS_PATH=/crex/proj/snic2019-35-58/water_strider/RNA-seq_data_analysis/GbRNAseqDEpAnnots/envs/
source conda_init.sh
conda init bash
conda activate collate_counts

STAR_MAPPING_OUTDIRS=submit_star_alignments
OUTFILE=all_ReadsPerGene.out.tsv

find ../"$STAR_MAPPING_OUTDIRS" -name "ReadsPerGene.out.tab" -print0 \
| xargs -0 ./collate_counts.py -o "$OUTFILE"

