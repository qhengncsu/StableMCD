library(StableMCD)
library(ggplot2)
library(gridExtra)
library(rrcov)
source("bootstrap_robpca.R")
data = read.csv("breast.csv")
x <- data.matrix(data[,1:2000])
res = PcaHubert(x,k=2,alpha=0.5)
n = dim(x)[1]
PCs = res$scores

ptm <- proc.time()
bootstrap_result1 = bootstrap_robpca(x,seq(0.5,0.975,by=0.025),2)
bootstrap_result2 = bootstrap_robpca(x,seq(0.5,0.975,by=0.025),10)
time <- proc.time() - ptm


data1 = data.frame(X1 = PCs[,1],X2 = PCs[,2], ER = data[["ER_status"]])
plot1 = ggplot()+ geom_point(data=data1, aes(x=X1,y=X2,color=ER),size=2) + 
  labs(y = "PC2", x = "PC1", title="Estrogen Receptor") + theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none", text = element_text(size=12)) +
  scale_color_manual(breaks=c("Positive","Negative"),
                     values=c("blue","red"))

data2 = data.frame(h=seq(0.5,0.975,by=0.025),mean_q2 = bootstrap_result1$means,sd_q2 = bootstrap_result1$sds,
                   mean_q10 = bootstrap_result2$means,sd_q10 = bootstrap_result2$sds)

plot2  = ggplot(data2,aes(x=h))+
  geom_point(aes(y=mean_q2,color="q=2"))+
  geom_point(aes(y=mean_q10,color="q=10"))+
  geom_errorbar(aes(ymin=mean_q2-sd_q2,ymax=mean_q2+sd_q2,color="q=2"),width=0.01)+
  geom_errorbar(aes(ymin=mean_q10-sd_q10,ymax=mean_q10+sd_q10,color="q=10"),width=0.01)+
  geom_line(aes(y=mean_q2,color="q=2"))+
  geom_line(aes(y=mean_q10,color="q=10"))+
  scale_color_manual("Number of PCs",breaks=c("q=2","q=10"),
                     values=c("navyblue","darkred"))+
  labs(y = "Instability", x = "h/n", title="Instability on Breast Cancer Data")+theme_bw()+
  theme(plot.title = element_text(hjust = 0.5), legend.position = "None", text = element_text(size=12))

result3 = mcd(PCs,0.775)
is_outlier = rep("High-SD",dim(x)[1])
ods = res$od>res$cutoff.od
is_outlier[result3$index] = "Inlier"
is_outlier[ods] = "High-OD"
data3 = data.frame(X1 = PCs[,1],X2 = PCs[,2], Class=is_outlier)
plot3 = ggplot()+ geom_point(data=data3, aes(x=X1,y=X2,color=Class),size=2) + 
  labs(y = "PC2", x = "PC1", title="Outlier Visualization (2PC)") + theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none", text = element_text(size=12)) +
  scale_color_manual(breaks=c("Inlier","High-SD","High-OD"),
                     values=c("blue","red","green"))

SDs = res$sd
ODs = res$od

cutoff.insta = sort(SDs)[floor(0.775*dim(x)[1])]

data4 = data.frame(X1 = SDs, X2=ODs)
plot4 = ggplot()+ geom_point(data=data4, aes(x=X1,y=X2),size=2, color='blue',shape=1) + 
  labs(y = "Orthogonal Distance", x = "Score Distance", title="Outlier Map (2PC) ") + theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none", text = element_text(size=12)) +
  geom_vline(data=data2, aes(xintercept=res$cutoff.sd),color='purple')+
  geom_hline(data=data2, aes(yintercept=res$cutoff.od),color='purple')+
  geom_vline(data=data2, aes(xintercept=cutoff.insta),color='red')

grid.arrange(plot1, plot2, plot3, plot4, nrow=2, ncol=2)

