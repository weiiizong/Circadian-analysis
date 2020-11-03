#library(devtools)
#install_github("Caleb-Huo/AWFisher") 
rm(list=ls())
library(AWFisher)
observed1 = read.csv("C:/Users/zongw/Box Sync/collaboration/Kyle/NAc 02032020/github_demo/data/Example_observed_NAc.csv",row.names = 1)
observed2 = read.csv("C:/Users/zongw/Box Sync/collaboration/Kyle/NAc 02032020/github_demo/data/Example_observed_Caudate.csv",row.names = 1)

genes = intersect(observed1$genes,observed2$genes)
index1 = match(genes,observed1$genes)
observed.m1 = observed1[index1,]
index2 = match(genes,observed2$genes)
observed.m2 = observed2[index2,]
all(as.character(observed.m1$genes) == as.character(observed.m2$genes))
p.value1= observed.m1$pvalue
p.value2= observed.m2$pvalue
pmatrix = cbind(p.value1,p.value2)
row.names(pmatrix) = genes

#AW-fisher to select genes that are rhythmic in both regions
library(AWFisher)
res = AWFisher_pvalue(pmatrix)
qvalue <- p.adjust(res$pvalue, "BH") 
meta_p_q = data.frame(pvalue = res$pvalues, qvalue = qvalue, row.names = row.names(pmatrix))
sigIndex = which(res$weights[,1]==1&res$weights[,2]==1)
meta_p_q = meta_p_q[sigIndex,]
meta_p_q_sort = meta_p_q[order(meta_p_q$qvalue,decreasing = F),]

length(which(meta_p_q_sort$qvalue<0.05))


dir="C:/Users/zongw/Box Sync/collaboration/Kyle/NAc 02032020/github_demo/output/Phaseshift"
dir.create(dir,recursive = T)
setwd(dir)
library(ggplot2)
cutoff = c(0.01,0.05)
for(i in 1:length(cutoff)){
  acutoff <- cutoff[i]
  sig_genes = rownames(meta_p_q_sort[which(meta_p_q_sort$qvalue<acutoff),])
  selectedGenes1 <- match(sig_genes, observed1$genes)
  selectedGenes2 <- match(sig_genes, observed2$genes)
  phase1 = observed1$phase[selectedGenes1]
  phase2 = observed2$phase[selectedGenes2]
  
  phaseDiff = phase1 - phase2
  
  phase1[which(phaseDiff>12)] = -(24 - phase1[which(phaseDiff>12)])
  phase2[which(phaseDiff< -12)] = -(24 - phase2[which(phaseDiff< -12)])
  phaseDiff2 = phase1 - phase2
  #phase & phase-24 are equivalent where phaseDiff should be 0. Increase gap by 24 if gap is small, and decrease gap by 24 if gap is large
  #After adjustment, maximum phase differenc could be 12,nminimum could be -12 (before is from -24 to 24)
  
  color = rep("black",length(sig_genes))
  p.df = data.frame(phase1, phase2, color)
  fileName <- paste("NAc_Caudate Phase concordance","METAq",acutoff,".pdf",sep="")
  pdf(fileName, width = 7, height = 7)
  amain <- paste('n = ', length(sig_genes),'. Phase Concordance=',
                 round(100*(length(which(abs(phaseDiff2)<=4))/length(sig_genes))
                       ,0),'% \nin (+/- 4hrs)',sep='')

  pp = ggplot(data = p.df,aes(x = phase1, y = phase2))  +
    geom_point()+
    theme(axis.ticks.x=element_line(colour="black",size=1,linetype=1),
          axis.text.x=element_text(colour="black",size=10,vjust=0),
          axis.ticks.y=element_line(colour="black",size=1,linetype=1),
          axis.text.y=element_text(colour="black",size=10,hjust=0)) +
    theme_gray(base_size = 22)+
    geom_abline(intercept=0,slope=1,colour='blue',size=1)+
    geom_abline(intercept=4,slope=1,colour='darkgreen',size=1,linetype=2)+
    geom_abline(intercept=-4,slope=1,colour='darkgreen',size=1,linetype=2)+
    labs(title = amain, x = "NAc", y = "Caudate")+
    xlim(-9.1,24) + ylim(-9.1,24)
  print(pp)
  dev.off()
}

