#!/bin/bash
#PBS -N peakDRIPc

#PBS -m a
#PBS -M maja@kuzman.org

#PBS -q MASTER
#PBS -l select=1:ncpus=2:mem=120G

cd $PBS_O_WORKDIR

## IN je iz ls DRIPcA*wig 


/common/WORK/mfabijanic/programs/DRIPc/peak_calling/DRIPc/normalize.pl $IN ${IN%%.wig}_normalized.wig
/common/WORK/mfabijanic/programs/DRIPc/peak_calling/DRIPc/wig2fa.pl -i ${IN%%.wig}_normalized.wig -o ${IN%%.wig}_normalized.customfa
/common/WORK/mfabijanic/programs/StochHMM/stochhmm -seq ${IN%%.wig}_normalized.customfa -model /common/WORK/mfabijanic/programs/DRIPc/peak_calling/DRIPc/DRIPc.hmm -posterior -threshold 0.9 -gff > ${IN%%.wig}.peaks.txt 


