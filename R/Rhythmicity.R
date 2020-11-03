library(limma)
library(edgeR)
library(doParallel)
rm(list=ls())
setwd("C:/Users/zongw/Box Sync/collaboration/Kyle/NAc 02032020/data")
data = read.csv("C:/Users/zongw/Box Sync/collaboration/Kyle/NAc 02032020/github_demo/data/Example_data.csv",row.names = 1)
clinical = read.csv("C:/Users/zongw/Box Sync/collaboration/Kyle/NAc 02032020/github_demo/data/Example_clinical.csv",row.names = 1)
tod = clinical$CorrectedTOD
names(tod) = clinical$pair
all(names(tod) == colnames(data))

#CIRCADIAN
library('minpack.lm')
source('C:/Users/zongw/Box Sync/collaboration/Kyle/NAc 02032020/github_demo/R/fitSinCurve.R')

n <- nrow(data)
genes<-row.names(data)
observed_para <- data.frame(genes=genes,A=numeric(n), phase=numeric(n), offset=numeric(n), peak=numeric(n), R2=numeric(n))

for (i in 1:n) {
  out <- fitSinCurve(xx=tod,observed=as.numeric(data[i,]))
  observed_para[i,-1] <- c(out$A, out$phase, out$offset, out$peak, out$R2)
  if(i%%1000==0) print(i)
}

obs_sorted<-observed_para[order(observed_para$R2, decreasing = TRUE),]

#PLOT
source("C:/Users/zongw/Box Sync/collaboration/Kyle/NAc 02032020/github_demo/R/circadianDrawing_axis.R")
WD<-"C:/Users/zongw/Box Sync/collaboration/Kyle/NAc 02032020/github_demo/output/Rhythmicity/TopPlots"
dir.create(WD,recursive = T)
setwd(WD)
tod <- tod
expr <- data
specInfo <- 'CONTROL'
labels <-rep(1, ncol(expr))
for(i in 1:10){
  agene <- as.character(obs_sorted$genes[i])
  index<-as.numeric(row.names(obs_sorted)[i])
  fileName <- paste0(specInfo,agene,'.pdf')
  pdf(fileName)
  circadianDrawing(tod=tod, expr=unlist(expr[index,]), apar=obs_sorted[i,],labels=labels, specInfo=specInfo)
  dev.off()
}

WD<-"C:/Users/zongw/Box Sync/collaboration/Kyle/NAc 02032020/github_demo/output/Rhythmicity/CorePlots"
dir.create(WD,recursive = T)
setwd(WD)
core_genes = c("ENSG00000133794", "ENSG00000134852", "ENSG00000008405",
               "ENSG00000121671", "ENSG00000126368", "ENSG00000174738",
               "ENSG00000179094", "ENSG00000132326", "ENSG00000049246")
#mapped symbols: c("ARNTL","CLOCK","CRY1", "CRY2","NR1D1","NR1D2","PER1", "PER2","PER3")
expr <- data[core_genes,]
row.names(observed_para)<-make.names(observed_para$genes, unique = TRUE)
obs_core<-observed_para[core_genes,]
specInfo <- 'CONTROL'
labels <-rep(1, ncol(expr))
for(i in 1:length(core_genes)){
  agene <- obs_core[i,"genes"]
  fileName <- paste0(specInfo,agene,'.pdf')
  pdf(fileName)
  circadianDrawing(tod=tod, expr=unlist(expr[i,]), apar=obs_core[i,],labels=labels, specInfo=specInfo)
  dev.off()
}

#Permutation to derive R^2 p-values
nullFolder <- "C:/Users/zongw/Box Sync/collaboration/Kyle/NAc 02032020/github_demo/output/Rhythmicity/null_control"
dir.create(nullFolder,recursive = T)
setwd(nullFolder)
groupName <- 'control'
thisData <- data
registerDoParallel()
B<-10#10 for demostration, 1000 was used in paper
#result <- foreach(b = 1:B) %dopar% {
for(b in 1:B){
  print(b)	
  library(minpack.lm)
  source('C:/Users/zongw/Box Sync/collaboration/Kyle/NAc 02032020/github_demo/R/fitSinCurve.R')
  
  null_pare <- data.frame(A=numeric(n), phase=numeric(n), offset=numeric(n), peak=numeric(n), R2=numeric(n))
  null_para_file <- paste('null_',groupName,'_',b,'.rdata',sep='')
  
  set.seed(b)
  shuffleTOD <- sample (tod)
  
  for (i in 1:n) {
    out <- fitSinCurve(xx=shuffleTOD,observed=unlist(thisData[i,]))
    null_pare[i,] <- c(out$A, out$phase, out$offset, out$peak, out$R2)
  }		
  save(null_pare,file=null_para_file)	
}
null_pare_A <- matrix(0,n,B)
null_pare_phase <- matrix(0,n,B)
null_pare_offset <- matrix(0,n,B)
null_pare_peak <- matrix(0,n,B)
null_pare_R2 <- matrix(0,n,B)

for(b in 1:B){
  print(b)
  file11 <- paste('null_',groupName,'_',b,'.rdata',sep='')
  
  load(file11)
  null_pare_A[,b] <- null_pare$A
  null_pare_phase[,b] <- null_pare$phase
  null_pare_offset[,b] <- null_pare$offset
  null_pare_peak[,b] <- null_pare$peak
  null_pare_R2[,b] <- null_pare$R2		
}

null_para <- list(null_para_A=null_pare_A, null_para_phase=null_pare_phase, null_para_offset=null_pare_offset, 
                  null_para_peak=null_pare_peak, null_para_R2=null_pare_R2)
null_para_file <- paste('null_',groupName,'.rdata',sep='')
save(null_para,file=null_para_file)
para_R2_pool <- c(observed_para$R2,null_para$null_para_R2)
R2Rank_para <- 1 - (rank(para_R2_pool)[1:length(observed_para$R2)] - 0.5)/length(para_R2_pool)
observed_para$pvalue <- R2Rank_para
observed_para$qvalue <- p.adjust(observed_para$pvalue, 'BH')
write.csv(observed_para,file = "C:/Users/zongw/Box Sync/collaboration/Kyle/NAc 02032020/github_demo/output/Rhythmicity/observed_para_demo.csv")

