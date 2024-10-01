library(StableMCD)
library(ggplot2)
library(gridExtra)
library(rrcov)
library(mrfDepth)
source("bootstrap_robpca.R")
data(glass)
glass = t(matrix(glass,nrow=750,ncol=180))

ptm <- proc.time()
bootstrap_result1 = bootstrap_robpca(glass,seq(0.5,0.975,by=0.025),2)
bootstrap_result2 = bootstrap_robpca(glass,seq(0.5,0.975,by=0.025),3)
bootstrap_result3 = bootstrap_robpca(glass,seq(0.5,0.975,by=0.025),10)
time <- proc.time() - ptm

data1 = data.frame(h=seq(0.5,0.975,by=0.025),mean_q2 = bootstrap_result1$means,sd_q2 = bootstrap_result1$sds,
                   mean_q3 = bootstrap_result2$means,sd_q3 = bootstrap_result2$sds,
                   mean_q10 = bootstrap_result3$means,sd_q10 = bootstrap_result3$sds)

plot1  = ggplot(data1,aes(x=h))+
  geom_point(aes(y=mean_q2,color="q=2"))+
  geom_point(aes(y=mean_q3,color="q=3"))+
  geom_point(aes(y=mean_q10,color="q=10"))+
  geom_line(aes(y=mean_q2,color="q=2"))+
  geom_line(aes(y=mean_q3,color="q=3"))+
  geom_line(aes(y=mean_q10,color="q=10"))+
  scale_color_manual("Number of PCs",breaks=c("q=2","q=3","q=10"),
                     values=c("navyblue","darkred","darkgreen"))+
  labs(y = "Instability", x = "h/n", title="Instability on Glass Data")+theme_bw()+
  theme(plot.title = element_text(hjust = 0.5), legend.position = "None", text = element_text(size=12))


robpca2 = PcaHubert(glass,k=3,alpha=0.5)
SDs = robpca2$sd
ODs = robpca2$od

cutoff.insta = sort(SDs)[floor(0.625*180)]

data2 = data.frame(X1 = SDs, X2=ODs)
plot2 = ggplot()+ geom_point(data=data2, aes(x=X1,y=X2),size=3, color='blue',shape=1) + 
  labs(y = "Orthogonal Distance", x = "Score Distance", title="Outlier Map (3PC) ") + theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none", text = element_text(size=12)) +
  geom_vline(data=data2, aes(xintercept=robpca2$cutoff.sd),color='purple')+
  geom_hline(data=data2, aes(yintercept=robpca2$cutoff.od),color='purple')+
  geom_vline(data=data2, aes(xintercept=cutoff.insta),color='red')

grid.arrange(plot1, plot2, nrow=1, ncol=2)
