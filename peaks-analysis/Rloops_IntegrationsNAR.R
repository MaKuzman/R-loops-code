library(data.table)
library(liftOver)
library(GenomicRanges)

#integrations in hg19: 

cur <- fread("/common/WORK/mfabijanic/Marina/hiv_allsamples_NAR_Guill_Eduard_HC2020_hg19.txt")
curgr <- GRanges(cur$chrom, IRanges(cur$pos, cur$pos), 
               strand=cur$strand, 
               sampleID=cur$rep, 
               barcode=cur$brcd,
               mapq=cur$mapq,
               cat=cur$cat,
               dna= cur$dna,
               rna=cur$rna,
               exprscore=cur$exprscore)

## liftover to hg38: 
path <- "/common/WORK/mfabijanic/Marina/hg19ToHg38.over.chain"
ch <- import.chain(path)
seqlevelsStyle(curgr) = "UCSC"  # necessary
cur38 = liftOver(curgr, ch)
integrations <- as.data.table(unlist(cur38))

# adding sample info
sampleInfo <- fread("/common/WORK/mfabijanic/Marina/hiv_allsamples_NAR_Guill_Eduard_HC2020_hg19_sampleinfo.csv")
setkey(sampleInfo, SampleID)
setkey(integrations, sampleID)
integrations <- sampleInfo[integrations]
fwrite(integrations, file="hiv_allsamples_NAR_Guill_Eduard_HC2020_hg38_withsampleinfo.txt")


# load R loops info:  
rloopsallfiles <- list.files(path="/common/WORK/mfabijanic/Marina/MarinaRloops/Macs2peaks/", 
                             pattern="xls", 
                             full.names=TRUE)

rloops <- lapply(rloopsallfiles,fread)

rloops <- do.call("rbind", rloops)
rloops[,c("drip","sample","time","input","peak","number"):=tstrsplit(name,"_",fixed=TRUE)]
rloopsgr <- GRanges(rloops$chr, IRanges(rloops$start, rloops$end), 
                    sample=rloops$sample, 
                    peak=rloops$number, 
                    pileup=rloops$pileup,
                    neglog10pval=rloops$`-log10(pvalue)`,
                    neglog10qval=rloops$`-log10(qvalue)`,
                    fold_enrichment= rloops$fold_enrichment
                    )

# Check GC content of Rloop sequences:  
library("BSgenome.Hsapiens.UCSC.hg38")
rloopseqs <- getSeq(BSgenome.Hsapiens.UCSC.hg38,rloopsgr)
onf1 <- oligonucleotideFrequency(rloopseqs,1)
onf2 <- oligonucleotideFrequency(rloopseqs,2)

rloops$GC_content <- (onf1[,"C"]+onf1[,"G"])/rowSums(onf1)
rloops$CpG_dinucl_content <- (onf2[,"CG"])/rowSums(onf2)


rloops[chr%in%paste("chr",c(1:22,"X","Y","x","y"), sep="")][time=="0h"][,
       .(meanGC=mean(na.omit(GC_content)),numberOfPeaks=.N,meanLength=mean(length)),
       .(sample,time,input)]

# distance to nearest for each sample:  
grintegrations<-GRanges(integrations)
grrloops <- GRanges(rloops[time=="0h"])
grrloops <- split(grrloops, grrloops$sample)

getdisttonear <- function(x,y){
    d <- as.data.frame(distanceToNearest(x,y, ignore.strand=TRUE))
    dist <- d$distance
    names(dist) <- d$queryHits
    dist
    }

dtn <- getdisttonear(grintegrations,grrloops[[1]])
integrations[(1:nrow(integrations))%in%as.numeric(names(dtn)),d_177:=dtn]

dtn <- getdisttonear(grintegrations,grrloops[[2]])
integrations[(1:nrow(integrations))%in%as.numeric(names(dtn)),d_178:=dtn]

dtn <- getdisttonear(grintegrations,grrloops[[3]])
integrations[(1:nrow(integrations))%in%as.numeric(names(dtn)),d_179:=dtn]



#plotting
library(scales)
library(ggplot2)

pdf("distanceToRloops.pdf", height=7, width=10)

intsJurkat <- integrations[Cells=="Jurkat"&Published=="Y"]
intsJurkat <- melt(intsJurkat,measure.vars=c("d_177","d_178","d_179"))
intsJurkat$Description <- factor(intsJurkat$Description, 
                                 levels=c("Jurkat 0x IC50 (438 ng/ul)", "Jurkat 2x IC50 (426 ng/ul)",
                                          "Jurkat 5x IC50 (646 ng/ul)", "Jurkat 10x IC50"))

ggplot(intsJurkat, aes(Description,value+0.0001,fill=Description))+
    stat_boxplot(geom="errorbar")+
    geom_boxplot(outlier.color=NA)+
    facet_grid(~variable)+
    theme_light()+
#    scale_y_continuous(trans = log2_trans())+
    theme(axis.text.x = element_text(angle = 90 ))+
    coord_cartesian(ylim=c(2,400000))+
    ylab("Distance to nearest R loop")+xlab("")+ggtitle("JURKAT")

intsSupT1 <- integrations[Cells=="SupT1"&Published=="Y"&(Treatment=="")]
intsSupT1 <- melt(intsSupT1,measure.vars=c("d_177","d_178","d_179"))


intsSupT1$Description <- factor(intsSupT1$Description, 
                                 levels=c("SupT1 0x IC50 (1060 ng/ul)","SupT1 2x IC50 (844 ng/ul)",
                                          "SupT1 5x IC50 (720 ng/ul)","SupT1 10x IC50 (714 ng/ul)"))

ggplot(intsSupT1, aes(Description,value+0.0001,fill=Description))+
    stat_boxplot(geom="errorbar")+
    geom_boxplot(outlier.color=NA)+
    facet_grid(~variable)+
    theme_light()+
#    scale_y_continuous(trans = log2_trans())+
    theme(axis.text.x = element_text(angle = 90 ))+
    coord_cartesian(ylim=c(2,420000))+
    ylab("Distance to nearest R loop")+xlab("")+ggtitle("SupT1")

dev.off()

pdf("distanceToRloops_SupT1_noidea.pdf", height=12, width=8)
intsSupT1_2 <- integrations[Cells=="SupT1"&Published=="Y"&(Treatment!="")]
intsSupT1_2 <- melt(intsSupT1_2,measure.vars=c("d_177","d_178","d_179"))
intsSupT1_2[,comment:=""]
intsSupT1_2[grepl("-LEDGIN",Description),comment:="-LEDGIN"]
intsSupT1_2[grepl("2XIC50",Description),comment:="2XIC50"]
intsSupT1_2[grepl("10XIC50",Description),comment:="10XIC50"]


intsSupT1_2$comment <- factor(intsSupT1_2$comment, 
                                levels=c("-LEDGIN","2XIC50","10XIC50"))

ggplot(intsSupT1_2, aes(comment,value+0.0001,fill=comment))+
    stat_boxplot(geom="errorbar")+
    geom_boxplot(outlier.color=NA)+
    facet_grid(Treatment~variable)+
    theme_light()+
#    scale_y_continuous(trans = log2_trans())+
    theme(axis.text.x = element_text(angle = 90 ))+
    coord_cartesian(ylim=c(2,500000))+
    ylab("Distance to nearest R loop")+xlab("")+ggtitle("SupT1")


dev.off()
