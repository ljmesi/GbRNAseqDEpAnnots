#!/bin/bash

set -e
set -uo pipefail

#./pull_simgs_docker.sh fastqc.sif quay.io/biocontainers/fastqc:0.11.9--hdfd78af_1
#./pull_simgs_docker.sh seqkit.sif quay.io/biocontainers/seqkit:2.0.0--h9ee0642_0
#./pull_simgs_docker.sh gffread.sif quay.io/biocontainers/gffread:0.12.7--hd03093a_1
# ./pull_simgs_docker.sh agat.sif quay.io/biocontainers/agat:0.9.2--pl5321hdfd78af_1
./pull_simgs_docker.sh multiqc.sif quay.io/biocontainers/multiqc:1.13a--pyhdfd78af_1

