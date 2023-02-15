library("BPRMeth")
library("ggplot2")
library("tidyverse")
args2 <- commandArgs(trailingOnly = TRUE)

Robj=args2[1]
load(Robj)
mparam=args2[2]
print(mparam)
kparam=args2[3]
outdirectory=args2[4]

############
# extract basename
name=tail(str_split(Robj,"/")[[1]],1)
name=str_split(name, "_")[[1]][1:3]
base=paste(name[1],name[2],name[3], sep="_")
otp.name=paste(outdirectory,"/",base,"_clusters_tmp.txt", sep="")

SelectedObj=paste("m",mparam,"k",kparam, sep="")
plotname=paste0(base,"_",SelectedObj)
# get selected model
SelectedObj=get(SelectedObj)

# make table
labels=cbind(values(met_region$anno), SelectedObj$labels)
write.table(x=labels, file=otp.name, quote=F,sep="\t")

plotname.pdf=paste0(outdirectory,"/",plotname,".pdf")
plotname.png=paste0(outdirectory,"/",plotname,".png")
p=plot_cluster_profiles(SelectedObj, title=plotname)
ggsave(plot=p,filename=plotname.pdf)
ggsave(plot=p,filename=plotname.png)
