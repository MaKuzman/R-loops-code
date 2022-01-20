## This folder contains all analysis of DRIP-seq and DRIPc-seq datasets.

Please use master files processOneSample_pairedEnd.sh and processOneSample_singleEnd.sh to execute the analysis.  

instructions:  

Run the script from a directory containing links to both or single ends of DRIPc paired or single end sequencing experiment
- make sure this directory contains only one DRIPc experiment!
- call this script with qsub -v GENOME=sth processOneSample_pairedEnd.sh
- GENOME variable contains the genome you want to map to. hg19 or hg38
-Please be careful not to run for multiple genome versions at the same time - it will not work because of dependencies in qsub jobs. If you wish to do this anyway, change the names of jobs and pids in this script to include genome names and make sure they are not too long 
- You can run for multiple different samples at the same time this is ok. Samples should be named DRIPsomething


