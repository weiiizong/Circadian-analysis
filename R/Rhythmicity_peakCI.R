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
n = nrow(data)

###Bootstrap
nullFolder <- "C:/Users/zongw/Box Sync/collaboration/Kyle/NAc 02032020/github_demo/output/Rhythmicity/bootstrap_control"
dir.create(nullFolder)
setwd(nullFolder)
groupName <- 'control'
thisData <- data
registerDoParallel()
B<-10
#result <- foreach(b = 1:B) %dopar% {
for(b in 1:B){
  print(b)	
  library(minpack.lm)
  source('C:/Users/zongw/Box Sync/collaboration/Kyle/NAc 02032020/github_demo/R/fitSinCurve.R')
  
  null_pare <- data.frame(A=numeric(n), phase=numeric(n), offset=numeric(n), peak=numeric(n), R2=numeric(n))
  null_para_file <- paste('boots_',groupName,'_',b,'.rdata',sep='')
  
  set.seed(b)
  resample.index = sample(1:length(tod),replace = T)
  
  for (i in 1:n) {
    out <- fitSinCurve(xx=tod[resample.index],observed=unlist(thisData[i,resample.index]))
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
  file11 <- paste('boots_',groupName,'_',b,'.rdata',sep='')
  
  load(file11)
  null_pare_A[,b] <- null_pare$A
  null_pare_phase[,b] <- null_pare$phase
  null_pare_offset[,b] <- null_pare$offset
  null_pare_peak[,b] <- null_pare$peak
  null_pare_R2[,b] <- null_pare$R2		
}

null_para <- list(null_para_A=null_pare_A, null_para_phase=null_pare_phase, null_para_offset=null_pare_offset, 
                  null_para_peak=null_pare_peak, null_para_R2=null_pare_R2)
null_para_file <- paste('bootstrap_',groupName,'.rdata',sep='')
save(null_para,file=null_para_file)

###Calculate peak CI
rm(list=ls())
load("C:/Users/zongw/Box Sync/collaboration/Kyle/NAc 02032020/github_demo/output/Rhythmicity/bootstrap_control/bootstrap_control.rdata")
obs = read.csv("C:/Users/zongw/Box Sync/collaboration/Kyle/NAc 02032020/github_demo/output/Rhythmicity/observed_para_demo.csv",row.names = 1)

Boots_peak_mat = null_para$null_para_peak
obs_peak = obs$peak
#for each gene, adjust raw bootstrap peaks to ones that is closest to the observed peak by +/- 24
peakDiff = sweep(x = Boots_peak_mat, MARGIN = 1, STATS = obs_peak, FUN = "-")
index1 = which(peakDiff>12,arr.ind = T)
index2 = which(peakDiff<(-12),arr.ind = T)
Boots_peak_mat1 = Boots_peak_mat
Boots_peak_mat1[index1] = Boots_peak_mat[index1]-24
Boots_peak_mat1[index2] = Boots_peak_mat[index2]+24
#5% to 95% quantile
Boots_peak_leftCI = apply(Boots_peak_mat1, 1, function(x) quantile(x,0.05))
Boots_peak_rightCI = apply(Boots_peak_mat1, 1, function(x) quantile(x,0.95))
#adjust CI boundaries to [-6,18]
transTtoInterval = function(t){
  t1 = ifelse(t<(-6),t+24,t)
  t2 = ifelse(t1>18,t1-24,t1)
  return(t2)
}
obs$peak_CI_left = transTtoInterval(Boots_peak_leftCI)
obs$peak_CI_right = transTtoInterval(Boots_peak_rightCI)
write.csv(obs,file = "C:/Users/zongw/Box Sync/collaboration/Kyle/NAc 02032020/github_demo/output/Rhythmicity/observed_para_peakCI_demo.csv")


