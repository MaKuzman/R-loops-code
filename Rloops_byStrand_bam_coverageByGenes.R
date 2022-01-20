library(data.table)
setwd("/common/WORK/mfabijanic/Marina/MarinaRloops/NewSequences")
library(liftOver)
library(biomaRt)        
library(GenomicDistributions)
library(ggplot2)

xforward <- fread("DRIPcAD1770h.sorted.nodup.forward.coverage1000bases.txt")
xreverse <- fread("DRIPcAD1770h.sorted.nodup.reverse.coverage1000bases.txt")
is <- fread("../is_hg38.txt")


load("../../MarinaIntegrationSites_nih/genidf.Robj")

nih <- readRDS("../../MarinaIntegrationSites_nih/IS_condensedSites.RDS")
path <- "/common/WORK/mfabijanic/Marina/hg19ToHg38.over.chain"
ch <- import.chain(path)
seqlevelsStyle(nih) = "UCSC"  # necessary
nih38 = liftOver(nih, ch)
nih38 <- as.data.table(unlist(nih38))


ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
genes <- as.data.table(getBM(attributes=c('ensembl_gene_id','hgnc_symbol',
                  'chromosome_name','start_position','end_position','strand','gene_biotype',
                  'ensembl_transcript_id','transcription_start_site',
                  'transcript_start','transcript_end'), mart = ensembl))
numexonsintranscripts <- as.data.table(getBM(attributes=c('ensembl_transcript_id','ensembl_exon_id'), mart = ensembl))
numexonsintranscripts<- numexonsintranscripts[,.N,ensembl_transcript_id]
setkey(genes, ensembl_transcript_id)
setkey(numexonsintranscripts,ensembl_transcript_id)
genes<-numexonsintranscripts[genes]
setnames(genes,"N","numberOfExons")
genes[,strand:=ifelse(strand==1,"+","-")]
genes <- genes[chromosome_name%in%as.character(1:24,"X","Y")]
genes[,chromosome_name:=paste0("chr",chromosome_name)]

peaksreverse <- fread("DRIPcAD1770h.sorted.nodup.reverse.overInput178reverse_peaks.narrowPeak")
peaksforward <- fread("DRIPcAD1770h.7orted.nodup.forward.overInput178forward_peaks.narrowPeak")
grreversepeaksRloops<-GRanges(peaksreverse$V1,
                              IRanges(peaksreverse$V2,peaksreverse$V3),
                              strand="-",
                              score=peaksreverse$V5)
grforwardpeaksRloops<-GRanges(peaksforward$V1,
                              IRanges(peaksforward$V2,peaksforward$V3),
                              strand="+",
                              score=peaksforward$V5)
grpeaksRloops <- c(grreversepeaksRloops,grforwardpeaksRloops)
grpeaksRloops[seqnames(grpeaksRloops)%in%paste("chr",c(1:23,"X","Y"), sep="")]
grpeaksRloops <- grpeaksRloops[seqnames(grpeaksRloops)%in%paste("chr",c(1:23,"X","Y"), sep="")]
grpeaksRloopsminus <- grpeaksRloops[strand(grpeaksRloops)=="-"]
grpeaksRloopsplus <- grpeaksRloops[strand(grpeaksRloops)=="+"]
grpeaksRloops <- GRangesList(grpeaksRloopsminus,grpeaksRloopsplus)
names(grpeaksRloops) <- c("Rloops minus strand","Rloops plus strand")
x = calcChromBinsRef(grpeaksRloops,"hg38")
plotChromBins(x)




getDistanceFor <- function(tss){
    tssdistanceplusRloops <- as.data.table(distanceToNearest(grpeaksRloopsplus, tss, ignore.strand=TRUE))
    tssdistanceminusRloops <- as.data.table(distanceToNearest(grpeaksRloopsminus, tss, ignore.strand=TRUE))
    tssdistanceminusRloops$geneOrientation <- as.character(strand(tss[tssdistanceminusRloops$subjectHits]))
    tssdistanceplusRloops$geneOrientation <- as.character(strand(tss[tssdistanceplusRloops$subjectHits]))
    tssdistanceminusRloops$TSS <- start(tss[tssdistanceminusRloops$subjectHits])
    tssdistanceplusRloops$TSS <- start(tss[tssdistanceplusRloops$subjectHits])
    tssdistanceplusRloops$Rloop <- start(grpeaksRloopsplus[tssdistanceplusRloops$queryHits])
    tssdistanceminusRloops$Rloop <- start(grpeaksRloopsminus[tssdistanceminusRloops$queryHits])
    
    
    tssdistanceminusRloops[,distancWithOrientationInMind:=distance*ifelse(geneOrientation=="+",
                                                                          ifelse(TSS<Rloop,+1,-1),
                                                                          ifelse(TSS>Rloop,+1,-1))]
    tssdistanceplusRloops[,distancWithOrientationInMind:=distance*ifelse(geneOrientation=="+",
                                                                         ifelse(TSS<Rloop,+1,-1),
                                                                         ifelse(TSS>Rloop,+1,-1))]
    tssdistanceminusRloops[,RloopOrientation:="-"]
    tssdistanceplusRloops[,RloopOrientation:="+"]
    tssdistdt <- rbind(tssdistanceminusRloops,tssdistanceplusRloops)
    TSSdist = list(`Rloops minus strand`=tssdistanceminusRloops$distancWithOrientationInMind,
                   `Rloops plus strand`=tssdistanceplusRloops$distancWithOrientationInMind)
    list(TSSdist,tssdistdt)
}



transcriptEnd <- GRanges(genes$chromosome_name, 
               #               IRanges(genes$transcription_start_site,genes$transcription_start_site), 
               IRanges(genes$transcript_end,genes$transcript_end), 
               strand=genes$strand,
               gene_biotype=genes$gene_biotype,
               ensembl_gene_id = genes$ensembl_gene_id,
               hgnc_symbol=genes$hgnc_symbol)
Tend <- getDistanceFor(transcriptEnd)
tss <- GRanges(genes$chromosome_name, 
                         #               IRanges(genes$transcription_start_site,genes$transcription_start_site), 
                         IRanges(genes$transcript_end,genes$transcript_end), 
                         strand=genes$strand,
                         gene_biotype=genes$gene_biotype,
                         ensembl_gene_id = genes$ensembl_gene_id,
                         hgnc_symbol=genes$hgnc_symbol)
TSSdist <- getDistanceFor(transcriptEnd)


plotFeatureDist(Tend[[1]], featureName="Transcript end")
plotFeatureDist(TSSdist[[1]], featureName="TSS")
plotFeatureDist(Tend[[1]], featureName="Transcript end", tile=TRUE)
plotFeatureDist(TSSdist[[1]], featureName="TSS", tile=TRUE)
d <- Tend[[2]]
d[,RloopOrientation:=ifelse(RloopOrientation=="-","Rloop on - strand","R loop on + strand")]
d[,geneOrientation:=ifelse(geneOrientation=="-","Gene on - strand","Gene on + strand")]

ggplot(d[distance<10000], aes(distancWithOrientationInMind))+
    geom_vline(xintercept=0, linetype="dashed", color="red") +
    geom_histogram(binwidth=150) + 
    xlab("Distance to TSS") + 
    facet_grid(geneOrientation~RloopOrientation)+
    theme_light()+
    coord_cartesian(xlim=c(-5000,5000)) 

dev.off()



######################### for each transcript
grtranscripts <- GRanges(genes$chromosome_name, 
                         IRanges(genes$transcript_start, 
                                 genes$transcript_end), 
                         strand=genes$strand,
                         ensembl_gene_id = genes$ensembl_gene_id,
                         hgnc_symbol = genes$hgnc_symbol,
                         gene_start = genes$start_position,
                         gene_end = genes$end_position,
                         gene_biotype = genes$gene_biotype,
                         ensembl_transcript_id = genes$ensembl_transcript_id,
                         numberOfExons = genes$numberOfExons,
                         transcription_start_site = genes$transcription_start_site
                         )


nn <- nearest(grtranscripts, grpeaksRloopsplus, ignore.strand=T)
transcripts <- cbind(as.data.table(grtranscripts), as.data.table(grpeaksRloopsplus)[nn])
colnames(transcripts)[14:19] <- paste("Rloops_plus", colnames(transcripts)[14:19], sep="_")
nn <- nearest(grtranscripts, grpeaksRloopsminus, ignore.strand=T)
transcripts <- cbind(transcripts, as.data.table(grpeaksRloopsminus)[nn])
colnames(transcripts)[20:25] <- paste("Rloops_minus", colnames(transcripts)[20:25], sep="_")


transcripts[,bothStrandRloop:=ifelse( 
    (Rloops_minus_start<Rloops_plus_end & Rloops_minus_start>Rloops_plus_start)|
        (Rloops_plus_start<Rloops_minus_end & Rloops_plus_start>Rloops_minus_start),
    "overlap", "no overlap")]


setkey(transcripts, seqnames,start,end)
is <- rbind(is,nih38[,.(seqnames,start,end,width=end-start+1,strand,Dataset=paste("nih",therapy,sep="_"))])
is[Dataset%in%c("nih_onART1","nih_onART2"),Dataset:="nih_onART"]

setkey(is, seqnames,start,end)
ovls <- foverlaps(transcripts, is)
ovls[!is.na(start),NisTranscript:=uniqueN(start),.(ensembl_transcript_id)]
ovls[is.na(start),NisTranscript:=0]
ovls[!is.na(start),NisGene:=uniqueN(start), ensembl_gene_id]
ovls[is.na(start),NisGene:=0]
ovls[!is.na(start),NlistsTranscript:=uniqueN(Dataset),.(ensembl_transcript_id)]
ovls[is.na(start),NlistsTranscript:=0]
ovls[!is.na(start),NlistsGene:=uniqueN(Dataset), ensembl_gene_id]
ovls[is.na(start),NlistsGene:=0]
ovls[,Ntranscripts:=uniqueN(paste(i.start, i.end,sep="-")),.(ensembl_gene_id)]
dcast(ovls[,.(.N,mean(Ntranscripts)),
           .(bothStrandRloop,NlistsGene)],NlistsGene+V2~bothStrandRloop,value.var=N)

ovls <- unique(ovls[,-c("start","end","width","strand","Dataset")])
genes_onlymaxistranscript <- ovls[order(-NisTranscript),.SD[1],ensembl_gene_id]
genes_onlymaxistranscript[,densityIDTranscript:=NisTranscript/i.width]

ggplot(genes_onlymaxistranscript,aes(RloopsMaxScore,densityIDTranscript))+
    geom_point()+stat_smooth()+theme_light()+facet_grid(RloopsMaxScoreStrand~i.strand)

ggplot(genes_onlymaxistranscript,aes(densityIDTranscript))+
    geom_histogram()+theme_light()+facet_grid(RloopsMaxScoreStrand~i.strand)
dev.off()



















































