args <- commandArgs(trailingOnly = TRUE)
library("MethylSeekR")
library("BSgenome.Cfamiliaris.UCSC.canFam3")
library("BSgenome.Hsapiens.UCSC.hg38")
library("BSgenome.Mmulatta.UCSC.rheMac10")
library("BSgenome.Mmusculus.UCSC.mm10")
library("BSgenome.Rnorvegicus.UCSC.rn6")


#####
## args[1]: species name - prefix
## args[2]: methylome
## args[3]: CpGislands
#####
 pref=args[1] # prefix
 meth=args[2] # methylome
 CGIs=args[3] # CpGislands
 OutDir=args[4] # otp dir

if(pref=="Human") {
    genome=Hsapiens
    sLengths=seqlengths(genome)
    sLengths=sLengths[1:23]
    names(sLengths)=sub("chr","",names(sLengths))
    seqnames(genome)=sub("chr","",seqnames(genome))
} else if(pref=="Macaque") {
    genome=Mmulatta
    sLengths=seqlengths(genome)
    sLengths=sLengths[1:21]
    names(sLengths)=sub("chr","",names(sLengths))
    seqnames(genome)=sub("chr","",seqnames(genome))
} else if(pref=="Mouse") {
    genome=Mmusculus
    sLengths=seqlengths(genome)
    sLengths=sLengths[1:20]
    names(sLengths)=sub("chr","",names(sLengths))
    seqnames(genome)=sub("chr","",seqnames(genome))
} else if(pref=="Rat") {
    genome=Rnorvegicus
    sLengths=seqlengths(genome)
    sLengths=sLengths[1:21]
    names(sLengths)=sub("chr","",names(sLengths))
    seqnames(genome)=sub("chr","",seqnames(genome))
} else if(pref=="Dog") {
    genome=Cfamiliaris
    sLengths=seqlengths(genome)
    sLengths=sLengths[1:39]
    names(sLengths)=sub("chr","",names(sLengths))
    seqnames(genome)=sub("chr","",seqnames(genome))
}


###

meth.gr<-readMethylome(FileName=meth, seqLengths=sLengths)
meth.gr=trim(meth.gr)

tab=read.table(CGIs)
tab=tab[,c(1,4,5,9)]
colnames(tab)=c("chr","start","end","id")
CpGislands.gr <-GRanges(tab)
CpGislands.gr <- suppressWarnings(resize(CpGislands.gr, 5000, fix="center"))


name=paste(c(OutDir, "/", pref, "_chr16_alphadist.pdf"), collapse="")
plotAlphaDistributionOneChr(m=meth.gr,chr.sel="16",num.cores=2, pdfFilename=name)


#FDR calculations
name=paste(c(OutDir, "/",pref, "_stats.pdf"), collapse="")
stats <- calculateFDRs(m=meth.gr, CGIs=CpGislands.gr, num.cores=2, pdfFilename=name )

FDR.cutoff<-5
m.sel<-0.5
n.sel=as.integer(names(stats$FDRs[as.character(m.sel), ][stats$FDRs[as.character(m.sel), ]<FDR.cutoff])[1])
n.sel

#save.image("prova.RData")
name=paste(c(OutDir, "/",pref, "_segmentUMRsLMRs.pdf"), collapse="")
UMRLMRsegments.gr <- segmentUMRsLMRs(m=meth.gr, meth.cutoff=m.sel, nCpG.cutoff=n.sel, num.cores=1, myGenomeSeq=genome, seqLengths=sLengths, pdfFilename=name)
#head(UMRLMRsegments.gr)

name=paste(c(OutDir, "/",pref, "_finalSegmentation.pdf"), collapse="")
plotFinalSegmentation(m=meth.gr, segs=UMRLMRsegments.gr,meth.cutoff=m.sel,pdfFilename=name )

nameGR=paste(c(OutDir, "/",pref, "_UMRsLMRs.gr.rds"), collapse="")
nameTab=paste(c(OutDir, "/",pref, "_UMRsLMRs.tab"), collapse="")
saveUMRLMRSegments(segs=UMRLMRsegments.gr, GRangesFilename=nameGR, TableFilename=nameTab)
