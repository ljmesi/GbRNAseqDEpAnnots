# Directory contents

In this directory both STAR genome index is created with `STARindex.sh` into 
`starIndex/` directory and STAR alignment of RNA-seq reads to genome are 
executed with `submit_star_alignmets.sh` script. The results of the alignemnt 
are output to directories named `cuttrim2_SF-2243-GB-1*`. Furthermore 
multiqc summary of STAR alignment statistics can be run with: 
`run_multiqc_onl_alns_stats.sh`. 
