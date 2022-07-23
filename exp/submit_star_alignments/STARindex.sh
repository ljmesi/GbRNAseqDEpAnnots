#!/bin/bash -l
#SBATCH -A snic2022-5-328
#SBATCH -o star_index_softmasked_%j.out
#SBATCH -J star_index_softmasked
#SBATCH -p core -c 12
#SBATCH -t 10:00:00
#SBATCH --mail-user=Lauri.Mesilaakso.5423@student.uu.se
#SBATCH --mail-type=ALL

module load bioinfo-tools
module load star/2.7.9a


STAR --runThreadN 12 \
--runMode genomeGenerate \
--genomeDir /proj/snic2019-35-58/water_strider/RNA-seq_data_analysis/GbRNAseqDEpAnnots/exp/submit_star_alignments/starIndex \
--genomeFastaFiles /proj/snic2019-35-58/water_strider/final_files/genome.softmasked.fa \
--sjdbGTFfile /proj/snic2019-35-58/water_strider/final_files/genome.genes.flybasewithcurated.gff3 \
--sjdbGTFtagExonParentTranscript Parent \
--genomeSAindexNbases 13


