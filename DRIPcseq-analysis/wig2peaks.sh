#!/bin/bash
#PBS -N strands

#PBS -m a
#PBS -M maja@kuzman.org

#PBS -q MASTER
#PBS -l select=1:ncpus=4:mem=50G

cd $PBS_O_WORKDIR

x=${x%%.forward.wig}

/common/WORK/mfabijanic/programs/DRIPc/peak_calling/DRIPc/normalize.pl ${x}.forward.wig ${x}.forward.norm.wig
/common/WORK/mfabijanic/programs/DRIPc/peak_calling/DRIPc/wig2fa.pl -i ${x}.forward.norm.wig -o ${x}.forward.norm.customfa

/common/WORK/mfabijanic/programs/StochHMM/stochhmm -seq ${x}.forward.norm.customfa -model /common/WORK/mfabijanic/programs/DRIPc/peak_calling/DRIPc/DRIPc.hmm -posterior -threshold 0.9 -gff > ${x}.forward.peaks.txt

/common/WORK/mfabijanic/programs/DRIPc/peak_calling/DRIPc/normalize.pl ${x}.reverse.wig ${x}.reverse.norm.wig
/common/WORK/mfabijanic/programs/DRIPc/peak_calling/DRIPc/wig2fa.pl -i ${x}.reverse.norm.wig -o ${x}.reverse.norm.customfa
/common/WORK/mfabijanic/programs/StochHMM/stochhmm -seq ${x}.reverse.norm.customfa -model /common/WORK/mfabijanic/programs/DRIPc/peak_calling/DRIPc/DRIPc.hmm -posterior -threshold 0.9 -gff > ${x}.reverse.peaks.txt

