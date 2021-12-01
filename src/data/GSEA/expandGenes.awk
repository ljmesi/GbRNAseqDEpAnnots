#!/usr/bin/awk -f
BEGIN{
    OFS="\t"
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
    "Genes",
    "GO.ID_"GOcat,
    "Term_"GOcat,
    "Annotated_"GOcat,
    "Significant_"GOcat,
    "Expected_"GOcat,
    "Rank in classicFisher_"GOcat,
    "classicFisher_"GOcat,
    "elimFisher_"GOcat,
    "weight01Fisher_"GOcat,
    "parentchildFisher_"GOcat,
    "weightFisher_"GOcat,
    "leaFisher_"GOcat;
}
{
    # Skip the header in the input
    if ($0 ~ /^GO.ID/) { 
        next;
    }
    $0=$0
    genesUnparsed = $13
    numGenes = split(genesUnparsed,genes2GO,",")
    for (gene in genes2GO) {
        print genes2GO[gene],$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12
    }
}