#!/bin/bash
# ----------------QSUB Parameters----------------- #
#PBS -q MASTER
#PBS -M maja.fabijanic@gmail.com
#PBS -m n
#PBS -N bam2bedgraphscaled
#PBS -l select=ncpus=1:mem=10g
#PBS -j oe
cd $PBS_O_WORKDIR

# ----------------THIS IS NOT A PIPELINE; NEEDS TO BE MANUALLY SUBMITTED AFTER EACH STEP----------------- #

# ----------------create links for files----------------- #
for x in `ls /common/RAW/Marina_DRIPseq/HY5JMBGXC_BL327_19s005254-1-1_Ibberson_lane1*`
do 
y=${x##/common/RAW/Marina_DRIPseq/HY5JMBGXC_BL327_19s005254-1-1_Ibberson_lane1}
ln -s $x ${y%%_sequence.txt.gz}.fastq.gz
done
# ----------------Trim the adapters, quality trimming and filtering with bbduk----------------- #
# Needs the adapters.fa file in the same directory
# Also requires 1_bbduk.sh file 
# trimming adapters + bases with quality lower then 10
# filtering all reads with avgquality <20
for y in `ls *fastq.gz | cut -d "_" -f1 | sort |uniq`
do 
qsub -v IN1=${y}_1.fastq.gz,IN2=${y}_2.fastq.gz,OUT1=${y}_1_TRIMMED.fastq,OUT2=${y}_2_TRIMMED.fastq,OUT_SINGLE=${y}_singles_TRIMMED.fastq 1_bbduk.sh
done
# ----------------After the trimming is done----------------- #
# ----------------FastQC of trimmed reads----------------- #

# ----------------Mapping----------------- #
bwa index hg19.fa

for y in `ls *fastq.gz | cut -d "_" -f1 | sort |uniq`
do 
qsub -v IN1=${y}_1_TRIMMED.fastq,IN2=${y}_2_TRIMMED.fastq,OUT=${y} -N ${y}_map 2_Mapping.sh
done


# ----------------Sam to sorted bam removed duplicates----------------- #
for y in `ls *sam | cut -d "." -f1 | sort |uniq`
do 
qsub -v IN=$y -N ${y}_samtosortbam /common/WORK/mfabijanic/programs/samToSortedBam.sh
done


#----------------separate strands---------------------------------------#
for i in `ls *nodup.bam`; do qsub -v i=$i -N ${i%%.sorted.nodup.bam}strands 3_strands.sh ; done
#### do separate mapping for single stranded sequencing (first round, sample 177)
ln -s ../../hg19/InputA_hg19.sorted.nodup.bam InputAD1770h.sorted.nodup.bam
ln -s ../../hg19/InputA_hg19.sorted.nodup.bam.bai InputAD1770h.sorted.nodup.bam.bai
qsub -v i=InputAD1770h.sorted.nodup.bam -N InputAD1770hstrands 3_strands.sh
ln -s ../../hg19/DRIPcA_hg19.sorted.nodup.bam DRIPcAD1770h.sorted.nodup.bam
ln -s ../../hg19/DRIPcA_hg19.sorted.nodup.bam.bai DRIPcAD1770h.sorted.nodup.bam.bai
qsub -v i=DRIPcAD1770h.sorted.nodup.bam -N DRIPcAD1770strands 3_strands.sh
 
#---------------------------bedgraph2wig------------------------------------#
for x in `ls *.sorted.nodup.forward.q30.bedgraph.sortd`
do

SPAN=10
S=$x
R=${x%%.sorted.nodup.forward.q30.bedgraph.sortd}.forward.wig
rm -f ${R}
head -45 ${S} | egrep "^browser|^track" > ${R}
grep "^chr" ${S} | cut -f1 | sort -u > chr.list
cat chr.list | while read C
do
    echo "variableStep chrom=${C} span=${SPAN}" >> ${R}
    awk '{if (match($1,"^'"${C}"'$")) { print } }' ${S} | sort -k2n | awk '
{
    printf "%d\t%g\n", $2+1, $4
}
' >> ${R}
done
done

for x in `ls *.sorted.nodup.reverse.q30.bedgraph.sortd`
do

SPAN=10
S=$x
R=${x%%.sorted.nodup.reverse.q30.bedgraph.sortd}.reverse.wig
rm -f ${R}
head -45 ${S} | egrep "^browser|^track" > ${R}
grep "^chr" ${S} | cut -f1 | sort -u > chr.list
cat chr.list | while read C
do
    echo "variableStep chrom=${C} span=${SPAN}" >> ${R}
    awk '{if (match($1,"^'"${C}"'$")) { print } }' ${S} | sort -k2n | awk '
{
    printf "%d\t%g\n", $2+1, $4
}
' >> ${R}
done
done

########## peak calling
for IN in `ls *wig | grep -v "norma"`
do
qsub -v IN=$IN -N ${IN%%.wig} peakCalling.sh
done

# ----------------Link it to hex----------------- #



http://hex.bioinfo.hr/~mfabijanic/is_hg19.bed
track type=bigWig name="177_0h_forward_DRIPcA" color=255,0,0 visibility=2 bigDataUrl="http://hex.bioinfo.hr/~mfabijanic/DRIPcAD1770h.sorted.nodup.forward.q30.bw"
track type=bigWig name="177_0h_reverse_DRIPcA" color=255,102,102 visibility=2 bigDataUrl="http://hex.bioinfo.hr/~mfabijanic/DRIPcAD1770h.sorted.nodup.reverse.q30.bw"
track type=bigWig name="178_0h_forward_DRIPcA" color=52,102,0 visibility=2 bigDataUrl="http://hex.bioinfo.hr/~mfabijanic/DRIPcAD1780h.sorted.nodup.forward.q30.bw"
track type=bigWig name="178_0h_reverse_DRIPcA" color=102,204,0 visibility=2 bigDataUrl="http://hex.bioinfo.hr/~mfabijanic/DRIPcAD1780h.sorted.nodup.reverse.q30.bw"
track type=bigWig name="179_0h_forward_DRIPcA" color=255,128,0 visibility=2 bigDataUrl="http://hex.bioinfo.hr/~mfabijanic/DRIPcAD1790h.sorted.nodup.forward.q30.bw"
track type=bigWig name="179_0h_reverse_DRIPcA" color=255,178,102 visibility=2 bigDataUrl="http://hex.bioinfo.hr/~mfabijanic/DRIPcAD1790h.sorted.nodup.reverse.q30.bw"
http://hex.bioinfo.hr/~mfabijanic/ConsensusPeaksFromTriplicates.txt
http://hex.bioinfo.hr/~mfabijanic/peaks177.txt
http://hex.bioinfo.hr/~mfabijanic/peaks178.txt
http://hex.bioinfo.hr/~mfabijanic/peaks179.txt

track type=bam name="RNA1_activatedctrl" bigDataUrl="http://hex.bioinfo.hr/~mfabijanic/mappedSRR8209209.sorted.nodup.bam"
track type=bam name="RNA2_activatedctrl" bigDataUrl="http://hex.bioinfo.hr/~mfabijanic/mappedSRR8209205.sorted.nodup.bam"
track type=bam name="RNA3_activatedctrl" bigDataUrl="http://hex.bioinfo.hr/~mfabijanic/mappedSRR8209207.sorted.nodup.bam"

track type=bam name="H3K36me3_198" bigDataUrl="http://hex.bioinfo.hr/~mfabijanic/mappedH3K36me3_198.sorted.nodup.bam"
track type=bam name="H3K9me2_197" bigDataUrl="http://hex.bioinfo.hr/~mfabijanic/mappedH3K9me2_197.sorted.nodup.bam"
track type=bam name="H3K36me3_197" bigDataUrl="http://hex.bioinfo.hr/~mfabijanic/mappedH3K36me3_197.sorted.nodup.bam"
track type=bam name="H3K9me2_198" bigDataUrl="http://hex.bioinfo.hr/~mfabijanic/mappedH3K9me2_198.sorted.nodup.bam"
track type=bam name="H3K27Ac_1" bigDataUrl="http://hex.bioinfo.hr/~mfabijanic/mappedSRR8235444.sorted.nodup.bam"
track type=bam name="H3K27Ac_2" bigDataUrl="http://hex.bioinfo.hr/~mfabijanic/mappedSRR8235445.sorted.nodup.bam"
track type=bam name="H3K27Ac_3" bigDataUrl="http://hex.bioinfo.hr/~mfabijanic/mappedSRR8235446.sorted.nodup.bam"

track type=bigWig name="177_18h_forward_DRIPcA" color=255,0,0 visibility=0 bigDataUrl="http://hex.bioinfo.hr/~mfabijanic/DRIPcAD17718h.sorted.nodup.forward.q30.bw"
track type=bigWig name="177_18h_reverse_DRIPcA" color=255,102,102 visibility=0 bigDataUrl="http://hex.bioinfo.hr/~mfabijanic/DRIPcAD17718h.sorted.nodup.reverse.q30.bw"
track type=bigWig name="177_18h_forward_DRIPcH" visibility=0 color=255,0,0 bigDataUrl="http://hex.bioinfo.hr/~mfabijanic/DRIPcHD17718h.sorted.nodup.forward.q30.bw"
track type=bigWig name="177_18h_reverse_DRIPcH" visibility=0 color=255,102,102 bigDataUrl="http://hex.bioinfo.hr/~mfabijanic/DRIPcHD17718h.sorted.nodup.reverse.q30.bw"
track type=bigWig name="178_0h_reverse_DRIPcH" visibility=0 color=192,192,192 bigDataUrl="http://hex.bioinfo.hr/~mfabijanic/DRIPcHD1780h.sorted.nodup.reverse.q30.bw"
track type=bigWig name="178_0h_forward_DRIPcH" visibility=0 color=192,192,192 bigDataUrl="http://hex.bioinfo.hr/~mfabijanic/DRIPcHD1780h.sorted.nodup.forward.q30.bw"
track type=bam name="177_0h_InputA" color=192,192,192 visibility=2 doWiggle=on bigDataUrl="http://hex.bioinfo.hr/~mfabijanic/InputA.sorted.nodup.bam"
track type=bam name="177_0h_InputH" color=192,192,192 visibility=2 doWiggle=on bigDataUrl="http://hex.bioinfo.hr/~mfabijanic/InputH.sorted.nodup.bam"
track type=bigWig name="179_0h_reverse_InputA" color=192,192,192 visibility=2 bigDataUrl="http://hex.bioinfo.hr/~mfabijanic/InputAD1790h.sorted.nodup.reverse.q30.bw"





track type=bigWig name="177_0h_DRIPA" bigDataUrl="http://hex.bioinfo.hr/~mfabijanic/DRIPA.sorted.nodup.reverse.bw"
track type=bigWig name="177_0h_DRIPH" bigDataUrl="http://hex.bioinfo.hr/~mfabijanic/DRIPH.sorted.nodup.reverse.bw"


track type=bigWig name="177_0h_DRIPcH" bigDataUrl="http://hex.bioinfo.hr/~mfabijanic/dripcH.sorted.nodup.reverse.bw"

http://hex.bioinfo.hr/~mfabijanic/DRIPcA_overInputA_peaks.gappedPeak.ucsc














# ----------------After mapping is done----------------- #
# ----------------MACS2----------------- #
















