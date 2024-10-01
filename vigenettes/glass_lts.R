library(ggplot2)
library(gridExtra)
library(dplyr)
library(StableMCD)
set.seed(1)
load('glass.rdata')
x = as.matrix(x[,c(25,106,230)])
y = y$PbO
alphas = seq(0.6,0.99,by=0.01)
bootstrap_result = bootstrap_lts(x,y,alphas,B=1000)

data = data.frame(alpha=alphas,insta_mean=bootstrap_result$means,insta_sd=bootstrap_result$sds)
ggplot(data,aes(x=alphas))+
  geom_point(aes(y=insta_mean,color="bootstrap"))+
  geom_line(aes(y=insta_mean,color="bootstrap"))+
  geom_errorbar(aes(ymin=insta_mean-insta_sd,ymax=insta_mean+insta_sd,color="bootstrap"),width=0.001)+
  scale_color_manual("Method",breaks=c("bootstrap"),
                     values=c("navyblue"))+
  labs(y = "Instability", x = "h/n", title="Instability Path")+theme_bw()+
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none", text = element_text(size=12))