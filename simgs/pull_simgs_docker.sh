#!/bin/bash
set -e
set -uo pipefail

NAME=$1
URL=$2

singularity pull "$NAME" docker://"$URL"

