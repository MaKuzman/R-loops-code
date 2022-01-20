#!/bin/bash
#PBS -N bbdukdedup

#PBS -m a
#PBS -M maja@kuzman.org 

#PBS -q MASTER
#PBS -l select=1:ncpus=6:mem=200G

cd $PBS_O_WORKDIR

/common/WORK/mfabijanic/programs/bbmap/bbduk.sh in=$IN1 in2=$IN2 out=$OUT1 out2=$OUT2 outs=$OUT_SINGLE threads=6 qtrim=rl trimq=10 minlength=50 minavgquality=20 -Xmx180g ref=/common/WORK/mfabijanic/Marina/MarinaRloops/scripts/adapters.fa ktrim=l 
