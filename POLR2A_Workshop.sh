#!/bin/bash

#SBATCH --account=def-training-wa_cpu --reservation=wgss1-wr_cpu

## Mail Options
#SBATCH --mail-user=lncao@sfu.ca
#SBATCH --mail-type=ALL

## CPU Usage
#SBATCH --mem-per-cpu=4G
#SBATCH --cpus-per-task=8
#SBATCH --time=3:00:00
#SBATCH --nodes=1

## Output and Stderr
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.error


# Requirements
module load bwa/0.7.15
module load samtools/1.9

# Conda activate
source /scratch/richmonp/TRAINING/JUNE2020/TOOLS/opt/miniconda3/etc/profile.d/conda.sh
conda activate /scratch/richmonp/TRAINING/JUNE2020/TOOLS/opt/SummerSchool/

# Genome
GENOME="/cvmfs/ref.mugqic/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa"
INDEX="/cvmfs/ref.mugqic/genomes/species/Homo_sapiens.GRCh38/genome/bwa_index/Homo_sapiens.GRCh38.fa"

# Globals
NAME="LNCAO"
ID="POLR2A"
WORKDIR=/scratch/richmonp/TRAINING/JUNE2020/$NAME/$ID
PLATFORM="Illumina"
SAMPLE="heart_left_ventricle"
THREADS=8
mkdir -p $WORKDIR
cd $WORKDIR

# Files
FASTQ_R1=$WORKDIR/$SAMPLE.$ID.R1.fastq.gz
FASTQ_R2=$WORKDIR/$SAMPLE.$ID.R2.fastq.gz
SAM_FILE=$WORKDIR/$SAMPLE.$ID.sam
BAM_FILE=$WORKDIR/$SAMPLE.$ID.bam
MACS2_DIR=$WORKDIR/MACS2
PEAKS_FILE=$MACS2_DIR/${SAMPLE}.${ID}.peaks.narrowPeak
CONTROL=$WORKDIR/${SAMPLE}.$ID.Control.bam
BIGWIG_FILE=$WORKDIR/$ID.$SAMPLE.bw

# Download ENCODE data
wget https://www.encodeproject.org/files/ENCFF999VOH/@@download/ENCFF999VOH.fastq.gz
wget https://www.encodeproject.org/files/ENCFF435KUB/@@download/ENCFF435KUB.fastq.gz
wget https://www.encodeproject.org/files/ENCFF925AMH/@@download/ENCFF925AMH.bam
mv $WORKDIR/ENCFF999VOH.fastq.gz $FASTQ_R1
mv $WORKDIR/ENCFF435KUB.fastq.gz $FASTQ_R2
mv $WORKDIR/ENCFF925AMH.bam $CONTROL

# Map reads to genome
if [ ! -f $SAM_FILE ]; then
	# bwa mem [options] <idxbase> <in1.fq> [in2.fq]
	bwa mem -t $THREADS -M \
	-R "@RG\tID:$ID\tSM:$SAMPLE\tPL:$PLATFORM" \
	$INDEX $FASTQ_R1 $FASTQ_R2 > $SAM_FILE
fi

# Convert SAM to sorted BAM
if [ ! -f $BAM_FILE ]; then
# create sorted BAM
	samtools view -Shu --threads `expr $THREADS - 1` $SAM_FILE | \
	samtools sort --threads `expr $THREADS - 1` -o $BAM_FILE -
	samtools index $BAM_FILE
fi

# BAM to BigWig
if [ ! -f $BIGWIG_FILE ]; then
	bamCoverage -b $BAM_FILE -of bigwig -p $THREADS -o $BIGWIG_FILE
fi
