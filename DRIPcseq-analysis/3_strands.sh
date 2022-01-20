#!/bin/bash
#PBS -N strands

#PBS -m a
#PBS -M maja@kuzman.org 

#PBS -q MASTER
#PBS -l select=1:ncpus=12:mem=120G

cd $PBS_O_WORKDIR
CHRSIZES=/common/WORK/mfabijanic/GENOMES/$GENOME.chrom.sizes
#for i in `ls *sorted.nodup.bam`
#do
bamCoverage --filterRNAstrand forward --binSize 10 --normalizeUsing RPKM --outFileFormat bedgraph -b $i -o ${i%%bam}forward.bedgraph -p 10 --minMappingQuality 30
bamCoverage --filterRNAstrand reverse --binSize 10 --normalizeUsing RPKM --outFileFormat bedgraph -b $i -o ${i%%bam}reverse.bedgraph -p 10 --minMappingQuality 30
bedSort ${i%%bam}reverse.bedgraph ${i%%bam}reverse.q30.bedgraph.sortd
bedSort ${i%%bam}forward.bedgraph ${i%%bam}forward.q30.bedgraph.sortd
bedGraphToBigWig ${i%%bam}reverse.q30.bedgraph.sortd $CHRSIZES ${i%%bam}reverse.q30.bw
bedGraphToBigWig ${i%%bam}forward.q30.bedgraph.sortd $CHRSIZES ${i%%bam}forward.q30.bw
#done

