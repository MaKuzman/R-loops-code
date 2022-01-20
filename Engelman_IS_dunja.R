library(data.table)
library(ggplot2)
library(GenomicRanges)
setwd("/common/WORK/mfabijanic/MarinaRloops/EngelmanIS")

geni <- readRDS(file="../HEK293Rloops/geni_hg19_with_numberofexons_rigs_rloops_cpsf6_ledgf_rloopsHEK.RDS")
genigr <- GRanges(geni)
hiv1 <- fread("/common/WORK/dglavas/Rloops/IS_sites/HIV1_IS_clean.txt")
hiv2 <- fread("/common/WORK/dglavas/Rloops/IS_sites/HIV2_IS_clean.txt")
mlv <- fread("/common/WORK/dglavas/Rloops/IS_sites/MLV_IS_clean.txt")
int <- fread("/common/WORK/dglavas/Rloops/IS_sites/INT_IS_clean.txt")
hiv2[,virustype:="HIV2"]
hiv1[,virustype:="HIV2"]
hiv1[,virustype:="HIV1"]
mlv[,virustype:="MLV"]
int[,virustype:="INT"]
allis <- rbind(hiv1,hiv2,mlv,int)
rloops <- readRDS("../NewSequences/hg19/consensusPeaks_withMeanSignal.RDS")
grrloops<-GRanges(rloops)
gris <- GRanges(allis$seqnames,IRanges(allis$IS_position,allis$IS_position))
allis$OverlapsRloop <- countOverlaps(gris,grrloops, ignore.strand=T)
allis$OverlapsgenethatovlsRloop <- countOverlaps(gris,genigr[genigr$rloop=="rloop"], ignore.strand=T)



ress <- allis[,.N,.(virustype,virus,celltype,OverlapsgenethatovlsRloop>0)]                
ress[,OverlapsgenethatovlsRloop:=ifelse(OverlapsgenethatovlsRloop==TRUE,"overlaps gene with CD4 R loop", "no overlap with gene with CD4 R loop")]   
ress[,type:=paste(paste(virustype,virus,sep="_"),celltype,sep=";")]
ggplot(ress, aes(type,N, fill=OverlapsgenethatovlsRloop))+theme_light()+geom_bar(stat="identity", position="fill")+ theme(axis.text.x = element_text(angle = 90), legend.position="bottom")

ress <- allis[,.N,.(virustype,virus,celltype,OverlapsRloop>0)]                
ress[,OverlapsRloop:=ifelse(OverlapsRloop==TRUE,"overlaps gene with CD4 R loop", "no overlap with gene with CD4 R loop")]   
ress[,type:=paste(paste(virustype,virus,sep="_"),celltype,sep=";")]
ggplot(ress, aes(type,N, fill=OverlapsRloop))+theme_light()+geom_bar(stat="identity", position="fill")+ theme(axis.text.x = element_text(angle = 90), legend.position="bottom")
dev.off()