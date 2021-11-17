#!/bin/bash

set -e
set -uo pipefail

singularity exec gffread.sif gffread -E -F -T genome.genes.gff3 > genome.genes.gtf

