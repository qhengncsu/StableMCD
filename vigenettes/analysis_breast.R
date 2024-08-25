source("MCD_spectral.R")
library(ggplot2)
library(gridExtra)
library(dplyr)
brca <- read.table("./data_mrna_seq_v2_rsem.txt")
brca <- t(brca)
brca <- data.frame(brca[3:1102,2:20532],row.names=brca[3:1102,1])
logtransform <- function(x){log2(x[1,1]+1)}
brca <- brca %>% mutate_all(list(~as.numeric(.)))
brca <- log2(brca+1)
sds <- apply(brca,2,sd)
order <- order(sds, decreasing = TRUE)
index <- order[1:2000]
brca <- brca[,index]
clinical <- read.table("./brca_tcga_clinical_data.tsv",sep = "\t", header=TRUE)
samples <- data.frame("Sample ID" = rownames(brca))
ER_status <- left_join(samples,clinical,by="Sample.ID")[["ER.Status.By.IHC"]]
brca["ER_status"] <- ER_status
brca <- brca[brca$ER_status %in% c("Positive","Negative"),]
x <- data.matrix(brca[,1:2000])
x_center = apply(x,2,mean)
x_centered = sweep(x,2,x_center)
svd_result = svd(x_centered)
PCs <- data.frame(x_centered%*%svd_result$v)
ptm <- proc.time()
result = bootstrap_insta(x,alpha_list=seq(0.5,0.95,by=0.05),q_list=c(2,10,100),B=50)
time <- proc.time() - ptm


spectral <- MCD_spectral(x,q=2,alpha=0.8,concentration="regular")
is_outlier_spectral = rep("outlier",dim(x)[1])
is_outlier_spectral[spectral$best] = "inlier"
data1 = data.frame(X1 = PCs[,1],X2 = PCs[,2], Class=is_outlier_spectral, ER = brca[["ER_status"]])
plot1 = ggplot()+ geom_point(data=data1, aes(x=X1,y=X2,color=ER),size=0.5) + 
  labs(y = "PC2", x = "PC1", title="Estrogen Receptor") + theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none", text = element_text(size=8)) +
  scale_color_manual(breaks=c("Positive","Negative"),
                     values=c("blue","red"))

data3 = data.frame(h=seq(0.5,0.95,by=0.05),mean_q2 = result$means[1,],sd_q2 = result$sds[1,],
                                          mean_q10 = result$means[2,],sd_q10 = result$sds[2,],
                                          mean_q100 = result$means[3,],sd_q100 = result$sds[3,])
plot2  = ggplot(data3,aes(x=h))+
  geom_point(aes(y=mean_q2,color="q=2"))+
  geom_point(aes(y=mean_q10,color="q=10"))+
  geom_point(aes(y=mean_q100,color="q=100"))+
  geom_line(aes(y=mean_q2,color="q=2"))+
  geom_line(aes(y=mean_q10,color="q=10"))+
  geom_line(aes(y=mean_q100,color="q=100"))+
  #geom_errorbar(aes(ymin=mean_q2-sd_q2,ymax=mean_q2+sd_q2,color="q=2"),width=0.01)+
  #geom_errorbar(aes(ymin=mean_q10-sd_q10,ymax=mean_q10+sd_q10,color="q=10"),width=0.01)+
  #geom_errorbar(aes(ymin=mean_q100-sd_q100,ymax=mean_q100+sd_q100,color="q=100"),width=0.01)+
  scale_color_manual("Number of PCs",breaks=c("q=2","q=10","q=100"),
                     values=c("navyblue","darkgreen","darkred"))+
  labs(y = "Instability", x = "h/n", title="Instability")+theme_bw()+
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none", text = element_text(size=8))

plot3 = ggplot()+ geom_point(data=data1, aes(x=X1,y=X2,color=Class),size=0.5) + 
  labs(y = "PC2", x = "PC1", title="SpectralMCD") + theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none", text = element_text(size=8)) +
  scale_fill_discrete(name = "Class") +
  scale_color_manual(breaks=c("inlier","outlier"),
                     values=c("blue","red"))

grid.arrange(plot1, plot2, plot3, nrow=1, ncol=3)
