#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -M maja.fabijanic@gmail.com
#PBS -m n
#PBS -N PEAKCALLINGPIPELINE
#PBS -l select=ncpus=1:mem=10g
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------Run this script from a directory containing links READS of DRIPc SINGLE end sequencing experiment----------------- #
# make sure this directory contains only one DRIPc experiment!
# !!!!!!!!!!! call this script with qsub -v GENOME=sth processOneSample_singleEnd.sh
# GENOME variable contains the genome you want to map to. hg19 or hg38
# !!!! DO NOT RUN FOR MULTIPLE GENOME VERSIONS IN THE SAME TIME - it will not work because of dependencies in qsub jobs - to make it work 
# change names of jobs and pids in this script to include genome names but check if they are too long ...
# you can run for multiple different samples at the same time this is ok. SAMPLES SHOULD BE NAMED DRIPsomething
 
# ----------------Trim the adapters, quality trimming and filtering with bbduk----------------- #
# Requires 1_bbduk.sh file (in the scripts directory)
# trimming adapters + bases with quality lower then 10
# filtering all reads with avgquality <20
y=`ls *fastq.gz | sort |uniq`
y=${y%%.fastq.gz}
qsub -v IN1=${y}.fastq.gz,OUT1=${y}_TRIMMED.fastq -N trim${y##DRIP} 1_bbduk_single.sh
pid=`qstat | grep trim${y##DRIP} | awk  {'print $1'}`

# ----------------After the trimming is done----------------- #
# ----------------Mapping to indexed genome ----------------- #

# GENOME variable contains the genome you want to map to. hg19 or hg38
#GENOME=hg19
ln -s /common/WORK/mfabijanic/GENOMES/$GENOME* . 
qsub -W depend=afterok:$pid -v GENOME=$GENOME,IN1=${y}_TRIMMED.fastq,OUT=${y} -N map${y##DRIP} 2_Mapping_single.sh
pid2=`qstat | grep map${y##DRIP} | awk  {'print $1'}`


# ----------------Sam to sorted bam removed duplicates----------------- #
qsub -W depend=afterok:$pid2 -v IN=$y -N s2sb${y##DRIP} /common/WORK/mfabijanic/programs/samToSortedBamSingle.sh
pid3=`qstat | grep s2sb${y##DRIP} | awk  {'print $1'}`


#----------------separate strands---------------------------------------#

qsub -W depend=afterok:$pid3 -v GENOME=$GENOME,i=${y}.sorted.nodup.bam -N strands${y##DRIP} 3_strands.sh 
pid4=`qstat | grep strands${y##DRIP} | awk  {'print $1'}`

 
#---------------------------bedgraph2wig------------------------------------#

qsub -W depend=afterok:$pid4 -v x=${y}.sorted.nodup.reverse.q30.bedgraph.sortd -N rev${y##DRIP} bedgraph2wigreverse.sh
qsub -W depend=afterok:$pid4 -v x=${y}.sorted.nodup.forward.q30.bedgraph.sortd -N fw${y##DRIP} bedgraph2wigforward.sh

pid5=`qstat | grep rev${y##DRIP} | awk  {'print $1'}`
pid6=`qstat | grep fw${y##DRIP} | awk  {'print $1'}`


########## peak calling
qsub -W depend=afterok:$pid5 -v IN=${y}.reverse.wig -N revPeak${y##DRIP} peakCalling.sh
qsub -W depend=afterok:$pid6 -v IN=${y}.forward.wig -N fwPeak${y##DRIP} peakCalling.sh
