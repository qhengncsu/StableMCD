library(StableMCD)
library(ggplot2)
library(gridExtra)
library(rrcov)
library(mrfDepth)
source("bootstrap_robpca.R")
set.seed(1)
data(octane)
octane = t(octane[,,1])
alphas = seq(21,38)/39 + 0.01
ptm <- proc.time()
bootstrap_result1 = bootstrap_robpca(octane,alphas,2,B=100)
bootstrap_result2 = bootstrap_robpca(octane,alphas,4,B=100)
bootstrap_result3 = bootstrap_robpca(octane,alphas,6,B=100)
time <- proc.time() - ptm

data1 = data.frame(h=seq(21,38),mean_q2 = bootstrap_result1$insta_means,sd_q2 = bootstrap_result1$insta_sds,
                   mean_q4 = bootstrap_result2$insta_means,sd_q4 = bootstrap_result2$insta_sds,
                   mean_q6 = bootstrap_result3$insta_means,sd_q6 = bootstrap_result3$insta_sds)

plot1  = ggplot(data1,aes(x=h))+
  geom_point(aes(y=mean_q2,color="q=2"))+
  geom_point(aes(y=mean_q4,color="q=4"))+
  geom_point(aes(y=mean_q6,color="q=6"))+
  geom_line(aes(y=mean_q2,color="q=2"))+
  geom_line(aes(y=mean_q4,color="q=4"))+
  geom_line(aes(y=mean_q6,color="q=6"))+
  scale_color_manual("Number of PCs",breaks=c("q=2","q=4","q=6"),
                     values=c("navyblue","darkred","darkgreen"))+
  labs(y = "Instability", x = "h", title="Instability on Octane Data")+theme_bw()+
  theme(plot.title = element_text(hjust = 0.5), legend.position = "None", text = element_text(size=12))


data2 = data.frame(h=seq(21,38),mean_q2 = bootstrap_result1$final_score,sd_q2 = bootstrap_result1$insta_sds,
                   mean_q4 = bootstrap_result2$final_score,sd_q4 = bootstrap_result2$insta_sds,
                   mean_q6 = bootstrap_result3$final_score,sd_q6 = bootstrap_result3$insta_sds)

plot2  = ggplot(data2,aes(x=h))+
  geom_point(aes(y=mean_q2,color="q=2"))+
  geom_point(aes(y=mean_q4,color="q=4"))+
  geom_point(aes(y=mean_q6,color="q=6"))+
  geom_line(aes(y=mean_q2,color="q=2"))+
  geom_line(aes(y=mean_q4,color="q=4"))+
  geom_line(aes(y=mean_q6,color="q=6"))+
  scale_color_manual("Number of PCs",breaks=c("q=2","q=4","q=6"),
                     values=c("navyblue","darkred","darkgreen"))+
  labs(y = "Instability", x = "h", title="Final Score on Octane Data")+theme_bw()+
  theme(plot.title = element_text(hjust = 0.5), legend.position = "None", text = element_text(size=12))


robpca3 = PcaHubert(octane,k=2,alpha=0.5)
SDs = robpca3$sd
ODs = robpca3$od

cutoff.insta = sort(SDs)[33]

data3 = data.frame(X1 = SDs, X2=ODs)
plot3 = ggplot()+ geom_point(data=data3, aes(x=X1,y=X2),size=3, color='blue',shape=1) + 
  labs(y = "Orthogonal Distance", x = "Score Distance", title="Outlier Map (2PC) ") + theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none", text = element_text(size=12)) +
  ylim(0,2)+
  geom_vline(data=data3, aes(xintercept=robpca3$cutoff.sd),color='purple')+
  geom_hline(data=data3, aes(yintercept=robpca3$cutoff.od),color='purple')+
  geom_vline(data=data3, aes(xintercept=cutoff.insta),color='red')

grid.arrange(plot1, plot2,plot3, nrow=1, ncol=3)

