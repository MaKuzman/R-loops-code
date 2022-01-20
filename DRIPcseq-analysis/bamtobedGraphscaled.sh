#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -M maja.fabijanic@gmail.com
#PBS -m n
#PBS -N bam2bedgraphscaled
#PBS -l select=ncpus=1:mem=10g
#PBS -J 0-6
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Loading variables------------------- #
CHR_LENGTH=/common/WORK/mfabijanic/GENOMES/$GENOME.chrom.sizes

INPUT_DIR=.
IN_SEQ=(`find $INPUT_DIR -name "*sorted.nodup.bam"`)
FILE=${IN_SEQ[$PBS_ARRAY_INDEX]}
BASE=${FILE#${INPUT_DIR}/}
BASE=${BASE%.bam}

# ----------------Commands------------------- #
# get number of mapped reads in millions
RPM=`cat *.flagstat.txt | grep "mapped (" | awk '{print $1}'`
RPM=`echo "scale=6; 1000000.0/"$RPM | bc`
BASE=${BASE}.scaled

# bam to scaled bedGraph, bedGraph to bigWig
genomeCoverageBed -ibam $FILE -bg -scale $RPM -split -g $CHR_LENGTH > ${BASE}.bedGraph
#wigToBigWig ${BASE}.bedGraph $CHR_LENGTH ${BASE}.bw
#[ -f "${BASE}.bw" ] && rm ${BASE}.bedGraph
