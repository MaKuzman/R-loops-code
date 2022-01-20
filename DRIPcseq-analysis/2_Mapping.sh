#!/bin/bash
#PBS -N map

#PBS -m a
#PBS -M maja@kuzman.org 

#PBS -q MASTER
#PBS -l select=1:ncpus=24:mem=120G

cd $PBS_O_WORKDIR

BWA=/common/WORK/mfabijanic/programs/miniconda3/bin/bwa
$BWA mem -t 22 ${GENOME}.fa $IN1 $IN2 > $OUT.sam