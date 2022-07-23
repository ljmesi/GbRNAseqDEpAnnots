#!/bin/bash

set -e
set -uo pipefail

singularity exec multiqc.sif multiqc .
