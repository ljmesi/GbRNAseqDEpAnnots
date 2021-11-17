#!/bin/bash -l
#SBATCH -A snic2021-5-195
#SBATCH -o star_mapping_softmasked_%j.out
#SBATCH -J star_mapping_softmasked
#SBATCH -p devel -c 12
#SBATCH -t 01:00:00
#SBATCH --mail-user=Lauri.Mesilaakso.5423@student.uu.se
#SBATCH --mail-type=ALL

module load bioinfo-tools
module load star/2.7.9a

SAMPLE_PREFIX=$1

REFERENCE=/proj/snic2019-35-58/water_strider/final_files/indices/GenomeDir-star-2.7.9a
ANNOTATION=/proj/snic2019-35-58/water_strider/RNA-seq_data_analysis/GbRNAseqDEpAnnots/exp/submit_star_alignments/genome.genes.gtf
THREADS=${SLURM_NTASKS:-12}
READS_DIR=/proj/snic2019-35-58/water_strider/RNA-seq_data_analysis/GbRNAseqDEpAnnots/data/raw

R1=$READS_DIR/${SAMPLE_PREFIX}_R1_001.fastq.gz
R2=$READS_DIR/${SAMPLE_PREFIX}_R2_001.fastq.gz

[[ -f $R1 && -f $R2 ]] || { echo "either '$R1' or '$R2' do not exist"; exit 1; }

mkdir ${SAMPLE_PREFIX}
cd ${SAMPLE_PREFIX}

STAR \
--runThreadN $THREADS \
--genomeDir "$REFERENCE" \
--readFilesIn "$R1" "$R2"  \
--readFilesCommand zcat \
--outSAMtype BAM SortedByCoordinate \
--quantMode GeneCounts \
--sjdbGTFfile "$ANNOTATION" \
--twopassMode Basic
