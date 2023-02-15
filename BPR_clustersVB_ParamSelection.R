library("BPRMeth")
library("ggplot2")
args <- commandArgs(trailingOnly = TRUE)

set.seed(147)
##### ARGS
out=args[1]
meth=args[2]
peak=args[3]
plotname=args[4]
imagename=args[5]
n=as.integer(args[6])
#plotname2=args[7]

imagename=paste(out,"/",imagename, sep="")

#######

met_dt <- read_met(file=meth, type = "bulk_seq")
anno_dt <- read_anno(file=peak, is_centre = T, is_window=T, upstream=n, downstream=n)

met_region <- create_region_object(met_dt = met_dt, anno_dt = anno_dt,sd_thresh=-1, cov=4)
save.image(imagename)

######
tab=c()
pk=c()

for (K in 1:10){
    K=as.integer(K)
    for (M in 1:10) {
        M=as.integer(M)
        basis_obj <- create_rbf_object(M = M)
        m=cluster_profiles_vb( X = met_region$met,
                                K = K,
                                model = "binomial",
                                basis = basis_obj,
                                vb_max_iter = 400)
        t=c(M,K, tail(m$lb,n=1))
        tab=rbind(tab,t)
        nam=paste("m", M,"k",K, sep="")
        assign(nam, m)

        pktmp=as.data.frame(as.matrix(m$pi_k))
        names=rownames(pktmp)
        rownames(pktmp)= NULL
        pktmp=cbind(names, pktmp)
        colnames(pktmp)=c("clusters","pi_k")
        pktmp$M=M
        pktmp$K=K
        pk=rbind(pk, pktmp)
        save.image(imagename)
    }
}


tab=as.data.frame(tab)
colnames(tab)=c("M", "K", "LB")
tab2=tab[which(tab$M<8),]

tab$K=as.factor(tab$K)
tab2$K=as.factor(tab2$K)
tab$M=as.factor(tab$M)
tab2$M=as.factor(tab2$M)

ms=ggplot(tab, aes(x=M, y= LB, group=K, colour=K)) +
  geom_line() +
  theme_bw()

pk=as.data.frame(pk)
colnames(pk)=c("clusters","pi_k","M", "K")
pk$K=as.factor(pk$K)
pk$M=as.factor(pk$M)

pkplot=ggplot(pk, aes(x=K, y=pi_k))+
  geom_bar(aes(fill = clusters), position = "dodge", stat="identity") +
  geom_hline(yintercept=0.1, linetype="dashed", color = "red") +
  theme_bw()+
  theme(axis.text = element_text(size = 13),
        legend.text = element_text(size=14),
        legend.title = element_text(size=15)) +
  facet_wrap(~M, ncol=3)


ks=ggplot(tab, aes(x=K, y= LB, group=M, colour=M)) +
  geom_line() +
  theme_bw()+
  theme(legend.position = "none",
        axis.text = element_text(size = 12))+
  facet_wrap(~M, ncol=3)

save.image(imagename)



plotname=paste(out,"/",plotname, sep="")
plotnameM=paste(plotname,"_",n,"_Mevolution_VB.png", sep="")
ggsave(plot=ms, filename=plotnameM, device="png", units = "in", width = 10, height = 10)

plotnamePK=paste(plotname,"_",n,"_piK_VB.png", sep="")
ggsave(plot=pkplot, filename=plotnamePK, device="png" , units = "in", width = 15, height = 10)

plotnameK=paste(plotname,"_",n,"_Kevolution_VB.png", sep="")
ggsave(plot=ks, filename=plotnameK, device="png", units = "in", width = 12, height = 7 )
