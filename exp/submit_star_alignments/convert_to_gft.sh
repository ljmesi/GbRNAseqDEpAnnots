#!/bin/bash

set -e
set -uo pipefail

singularity exec agat.sif agat_convert_sp_gff2gtf.pl --gff genome.genes.flybasewithcurated.gff3 -o genome.genes.flybasewithcurated.gtf

