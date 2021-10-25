from pathlib import Path
from snakemake.utils import min_version
min_version("6.10.0")

configfile: "config.yaml"

report: "report/workflow.rst"

include: "rules/common.smk"