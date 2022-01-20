library(data.table)
library(GenomicRanges)
library(ggplot2)
setDTthreads(6)

lapply_pb <- function(X, FUN, ...)
{
  env <- environment()
  pb_Total <- length(X)
  counter <- 0
  pb <- txtProgressBar(min = 0, max = pb_Total, style = 3)   
  
  # wrapper around FUN
  wrapper <- function(...){
    curVal <- get("counter", envir = env)
    assign("counter", curVal +1 ,envir=env)
    setTxtProgressBar(get("pb", envir=env), curVal +1)
    FUN(...)
  }
  res <- lapply(X, wrapper, ...)
  close(pb)
  res
}


setwd("/common/WORK/mfabijanic/Marina/MarinaRloops/NewSequences/hg19")
x177f <- fread('grep \"chr\" DRIPcAD1770h/DRIPcAD1770h.forward.peaks.txt')
x177f <- x177f[,.(chr=V1,start=V4,end=V5,strand="+",type=V3,sample="177")]
x177r <- fread('grep \"chr\" DRIPcAD1770h/DRIPcAD1770h.reverse.peaks.txt')
x177r <- x177r[,.(chr=V1,start=V4,end=V5,strand="-",type=V3,sample="177")]
x178f <- fread('grep \"chr\" DRIPcAD1780h/DRIPcAD1780h.forward.peaks.txt')
x178f <- x178f[,.(chr=V1,start=V4,end=V5,strand="+",type=V3,sample="178")]
x178r <- fread('grep \"chr\" DRIPcAD1780h/DRIPcAD1780h.reverse.peaks.txt')
x178r <- x178r[,.(chr=V1,start=V4,end=V5,strand="-",type=V3,sample="178")]
x179f <- fread('grep \"chr\" DRIPcAD1790h/DRIPcAD1790h.forward.peaks.txt')
x179f <- x179f[,.(chr=V1,start=V4,end=V5,strand="+",type=V3,sample="179")]
x179r <- fread('grep \"chr\" DRIPcAD1790h/DRIPcAD1790h.reverse.peaks.txt')
x179r <- x179r[,.(chr=V1,start=V4,end=V5,strand="-",type=V3,sample="179")]

allsamples <- rbind(x177f,x177r, x178f,x178r,x179f,x179r)
allsamples <- allsamples[(end-start+1)>100]
ass <- allsamples[1:nrow(allsamples),.(chr,start,end,name=paste(paste(type,sample,sep="_"),1:nrow(allsamples),sep="_"), score=1000,strand)]
split(ass, allsamples$sample)
ll <- split(ass, allsamples$sample)

fwrite(list('track name="DRIPcAD177 peaks" color=255,0,0 visibility=2'), file="peaks177.txt", quote=FALSE)
fwrite(list('track name="DRIPcAD178 peaks" color=52,102,0 visibility=2'), file="peaks178.txt", quote=FALSE)
fwrite(list('track name="DRIPcAD179 peaks" color=255,128,0 visibility=2'), file="peaks179.txt", quote=FALSE)
fwrite(ll[[1]], file="peaks177.txt", sep="\t", col.names=FALSE, append=T)
fwrite(ll[[2]], file="peaks178.txt", sep="\t", col.names=FALSE, append=T)
fwrite(ll[[3]], file="peaks179.txt", sep="\t", col.names=FALSE, append=T)

rm(list=grep("x17",ls(), value=T))
rm(ls, ass)

##########################################
# Defining consensus peaks from replicates
##########################################
gr <- GRanges(allsamples)
reduced <- reduce(gr)
x <- split(gr, gr$sample)
ols <- sapply(x,function(y)countOverlaps(reduced,y)>0)
reduced$s177 <- ols[,1]
reduced$s178 <- ols[,2]
reduced$s179 <- ols[,3]
reduced$totalSupport <- rowSums(ols)
rm(x,ols)


rdt <- as.data.table(reduced)
saveRDS(rdt, file="allConsensusPeaks_withsampleinfo.RDS")
goodones <- rdt[totalSupport>1]
consensusPeaks <- goodones[,.(seqnames,start,end,name=paste(totalSupport,1:nrow(goodones),sep="_"), 
                              score=ifelse(totalSupport==2,666,1000),strand)]
fwrite(list('track name="Consensus peaks" color=0,0,0 visibility=2'), file="ConsensusPeaksFromTriplicates.txt", quote=FALSE)
fwrite(consensusPeaks, file="ConsensusPeaksFromTriplicates.txt", sep="\t", col.names=FALSE, append=T)


#################################################################################
# Once we have consensus peaks: 
#  consensusPeaks <- fread("ConsensusPeaksFromTriplicates.txt", skip=1)
# setnames(consensusPeaks, 1:6, c("seqnames","start","end","name","score","strand"))
consensusPeaks[,totalSupport:=stringr::str_extract(name,"\\d")]
consensusPeaks[,width:=end-start+1]

pdf(file = "PeakSizeDistribution.pdf", width=12, height=12)
ggplot( consensusPeaks, aes(width/1000)) + 
  geom_histogram(alpha=0.8, binwidth=0.5, position = "stack", color="black") + 
  theme_light()+
  coord_cartesian(xlim=c(0,15.5))+
  theme(legend.position="bottom",text = element_text(size=30)) +
  xlab("DRIPc peak size (kb)")+
  ylab("Count")
dev.off()

pdf(file = "PeakSizeDistributionByReplicates.pdf", width=12, height=12)
ggplot( consensusPeaks, aes(width/1000, group=totalSupport, fill=factor(totalSupport))) + 
  geom_histogram(alpha=0.8, binwidth=0.5, position = "stack") + 
  theme_light()+
  coord_cartesian(xlim=c(0,15.5))+
  theme(legend.position="bottom",text = element_text(size=30)) +
#  scale_fill_discrete(name = "How many replicates support the peak") +
  scale_fill_brewer(palette="Dark2", name = "Number of replicates supporting the peak") + 
  xlab("DRIPc peak size (kb)")+
  ylab("Count")
dev.off()



#----------------------------------------#
# How correlated are the samples?
#----------------------------------------#

allfiles <- list.files(pattern="DRIPcA.*0h.sorted.nodup.*q30.bedgraph.sortd", recursive=T)
allfs <- lapply(allfiles, fread)
names(allfs) <- stringr::str_remove(stringr::str_remove(allfiles,".sorted.nodup"),".q30.bedgraph.sortd")
brred <- sapply(allfs, nrow)
allfiles <- do.call("rbind", allfs)
allfiles[,sample:=rep(names(brred), brred)]
allfiles[,scaleFactor:=sum(V4*(V3-V2+1))/8797954785,sample]
allfiles[,scaledRPKM:=V4/scaleFactor]
setnames(allfiles, 1:4,c("seqnames","start","end","signal"))
samples <- split(allfiles,allfiles$sample)



#-----------------------------------------------------------------------#
# merging with consensusPeaks and calculating mean signal over each peak
setkey(consensusPeaks,seqnames,start, end)

calcMeanSignal <- function(x){
  y <- x[,.(seqnames,start,end,scaledRPKM,sample)]
  setkey(y,seqnames,start, end)
  ols <- foverlaps(consensusPeaks,y)
  ols[,':='(modstart=ifelse(start<i.start,i.start,start), 
            modend=ifelse(end>i.end,i.end,end))]
  ols[,':='(doit=scaledRPKM*(modend-modstart+1),totlen=i.end-i.start+1)]
  ols[,.(meanSignalOnPeak=sum(doit/totlen)),.(seqnames, start=i.start,end=i.end)]
  
}
allsignals <- lapply_pb(samples, calcMeanSignal)
names(allsignals) <- stringr::str_replace(names(allsignals),".*/","")
# saveRDS(allsignals, file="allsignals.RDS")
# allsignals <- readRDS("allsignals.RDS")
sigs <- sapply(allsignals,function(x)x$meanSignalOnPeak)

sigs <- as.data.table(sigs)
sigs[,':='(seqnames=allsignals$DRIPcAD1790h.reverse$seqnames, start=allsignals$DRIPcAD1790h.reverse$start, end=allsignals$DRIPcAD1790h.reverse$end)]
setkey(sigs,seqnames,start,end)
setkey(consensusPeaks,seqnames,start,end)
consensusPeaks <- sigs[consensusPeaks]
consensusPeaks[,':='(s177=ifelse(strand=="+",DRIPcAD1770h.forward,DRIPcAD1770h.reverse),
                     s178=ifelse(strand=="+",DRIPcAD1780h.forward,DRIPcAD1780h.reverse),
                     s179=ifelse(strand=="+",DRIPcAD1790h.forward,DRIPcAD1790h.reverse))]
consensusPeaks <- consensusPeaks[,.(seqnames,start,end,name, score,strand,s177,s178,s179)]
saveRDS(consensusPeaks, file="consensusPeaks_withMeanSignal.RDS")

library(corrplot)
colnames(sigs)[1:6] <- c("D1 +", "D1 -", "D2 +", "D2 -","D3 +", "D3 -")
M <- cor(sigs[,1:6], method="spearman")
res1 <- cor.mtest(sigs[,1:6], conf.level = .95, method="spearman")
# everything is super significant p<10e-16
pdf("CorrelationBetweenReplicates.pdf", height=6, width=6)
corrplot(M, p.mat = res1$p, sig.level = .2, order="hclust", type="upper",addCoef.col = "black")
dev.off()


library(UpSetR)
rdt[,.N,.(s177,s178,s179)][order(-s177,-s178,-s179)]
#     s177  s178  s179     N
# 1:  TRUE  TRUE  TRUE 12608
# 2:  TRUE  TRUE FALSE  3758
# 3:  TRUE FALSE  TRUE   835
# 4:  TRUE FALSE FALSE  3465
# 5: FALSE  TRUE  TRUE  7171
# 6: FALSE  TRUE FALSE 98786
# 7: FALSE FALSE  TRUE  5427

pdf("UpSetPlot_peak_overlaps.pdf", height=5, width=8)
expressionInput <- c(D1 = 3465, D2 = 98786, D3 = 5427, `D1&D2` = 3758, `D1&D3` = 835, 
                     `D2&D3` = 7171, `D1&D2&D3` = 12608)
upset(fromExpression(expressionInput))
dev.off()


##########################################
# Getting genes from ensembl biomart
##########################################

# GRCh38.p13

library(biomaRt)
ensembl <- useEnsembl(biomart = "ensembl",
                      dataset = "hsapiens_gene_ensembl")
listAttributes(ensembl)[1:40,1:2]
atts <- c("ensembl_gene_id", "ensembl_transcript_id", "description", "chromosome_name", "start_position", 
          "end_position", "strand", "transcript_start", "transcript_end", "transcription_start_site", 
          "gene_biotype", "transcript_biotype",
          "hgnc_symbol")
genes <- as.data.table(getBM(attributes = atts, mart = ensembl))
goodchroms <- genes[,.N,chromosome_name][order(-N)][1:23,chromosome_name]
genes <- genes[chromosome_name%in%goodchroms]
genes[,seqnames:=paste("chr",chromosome_name,sep="")]
genes[,strand:=ifelse(strand==1,"+","-")]

## biotype grouping together, take only those transcripts which have same biotype as gene
genes[, biotype:=gene_biotype]
genes[grepl("IG_",gene_biotype), biotype:="IG"]
genes[grepl("TR",gene_biotype), biotype:="TR"]
genes[grepl("pseudogene",gene_biotype), biotype:="pseudogene"]
genes[gene_biotype%in%c("rRNA","ribozyme","scaRNA","sRNA","scRNA","vault_RNA"), biotype:="misc_RNA"]
#genes <- genes[gene_biotype==transcript_biotype]
genes[,start:=transcript_start]
genes[,end:=transcript_end]
setkey(genes, seqnames, start,end)
### samo geni 
geni <- genes[,.N,.(ensembl_gene_id,seqnames, start=start_position,end=end_position,biotype, strand, symbol=hgnc_symbol)]


################################################################################
# spajanje s genidf hg19 tablicom
# 
genidf <- readRDS("genidf_hg19.RDS")
setnames(genidf,  names(genidf), paste("hg19", names(genidf), sep="_"))
setnames(genidf,"hg19_names","ensembl_gene_id")
setkey(genidf,ensembl_gene_id)
setkey(geni,ensembl_gene_id)
genidfhg38 <- merge(geni, genidf)
notinhg38 <- genidf[!geni]
cn <- colnames(genidfhg38)
genidf[,symbol:=hg19_symbol]
setkey(genidf,symbol)
setkey(geni,symbol)
missinginhg38_symbols <- notinhg38[,.N,hg19_symbol][N==1,hg19_symbol]
missinginhg38_bysymbols <- unique(merge(geni[symbol%in%missinginhg38_symbols],genidf[symbol%in%missinginhg38_symbols]))
setnames(missinginhg38_bysymbols,"ensembl_gene_id.x","ensembl_gene_id")
genidfhg38 <- rbind(genidfhg38,missinginhg38_bysymbols[,..cn])
limits <- quantile(genidfhg38[biotype=="protein_coding"][hg19_i.Mean>0, hg19_i.Mean], c(0.1, 0.9))
genidfhg38[, ExpressionGroups:=factor(ifelse(hg19_i.Mean<=0,"Not expressed", 
                                         ifelse(hg19_i.Mean<=limits[1], "Low10%",
                                                ifelse(hg19_i.Mean<=limits[2], "Mid","High10%")
                                         )
),
levels=c("Not expressed", "Low10%", "Mid", "High10%")     
)]
#saveRDS(genidfhg38, file="genidf_hg38.RDS")

# 
################################################################################



setkey(consensusPeaks,seqnames,start,end)
setkey(genidfhg38,seqnames,start,end)
ols <- foverlaps(genidfhg38, consensusPeaks, mult="all")
ols[,direction:=ifelse(strand==i.strand,"sense","antisense")]
ols[,isRIG:=ifelse(hg19_RIGSLISTE=="2 or more","RIG","not RIG")]
ols[,RIGSliste:=ifelse(hg19_RIGSLISTE=="2 or more","RIG","other protein coding")]
xrig <- ols[!is.na(direction)&biotype=="protein_coding",.N,.(direction,biotype,RIGSliste)][order(-N)]
xrig$direction <- factor(xrig$direction, levels=c("sense", "antisense"))

x <- ols[!is.na(direction),.N,.(direction,biotype)][order(-N)]
fo <- unique(ols[!is.na(direction),.N,.(direction,biotype)][order(-N),biotype])
x$biotype <- factor(x$biotype, levels = fo) 
x[,biotype2:=biotype]
x[!biotype%in%c("protein_coding","lncRNA"),biotype2:="other"]
library(RColorBrewer)
x[,N2:=sum(N),.(direction,biotype2)]

x$direction <- factor(x$direction, levels=c("sense", "antisense"))
x$direction2 <- factor(as.character(x$direction))
pdf("Sense_antisense_Rloops_RIGS_protcodgenes.pdf", height=6, width=8)
ggplot(xrig, aes(RIGSliste,N, fill=direction)) + 
  geom_bar(stat="identity", position = "dodge") + 
  theme_light() + 
  theme(axis.text.x = element_text(angle = 90), text = element_text(size=20)) +
  scale_fill_brewer(palette="Dark2") +
  ylab("Number of R loops") + xlab("Gene type") +
  geom_text(aes(label=N), position=position_dodge(width=0.9), vjust=-0.25)
dev.off()

pdf("Sense_antisense_Rloops.pdf", height=8, width=12)
ggplot(x, aes(biotype2,N2, fill=direction)) + 
  geom_bar(stat="identity", position = "dodge") + 
  theme_light() + 
  theme(axis.text.x = element_text(angle = 90), text = element_text(size=20)) +
  scale_fill_brewer(palette="Dark2") +
  ylab("Number of R loops") + xlab("Gene type") +
  geom_text(aes(label=N2), position=position_dodge(width=0.9), vjust=-0.25)
ggplot(xrig, aes(RIGSliste,N, fill=direction)) + 
  geom_bar(stat="identity", position = "fill") + 
  theme_light() + 
  theme(axis.text.x = element_text(angle = 90), text = element_text(size=20)) +
  scale_fill_brewer(palette="Dark2") +
  coord_polar(theta = "y")
  
ggplot(x, aes(biotype,N, fill=direction2)) + 
  geom_bar(stat="identity", position="fill") + 
  theme_light() + 
  theme(axis.text.x = element_text(angle = 90), text = element_text(size=20)) +
  ylab("Direction of R loops as percent of total R loops
             overlapping gene type") + xlab("Gene type")+
  scale_y_continuous(labels = scales::percent)+
  scale_fill_manual(name = "direction", values = c(rev(brewer.pal(n = 3, name = "Dark2")[1:2])))

dev.off()

