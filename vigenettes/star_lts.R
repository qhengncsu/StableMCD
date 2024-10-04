library(ggplot2)
library(gridExtra)
library(dplyr)
library(robustbase)
library(StableMCD)
set.seed(1)
data("starsCYG")
alphas = seq(25,46)/47+0.001
x = matrix(starsCYG$log.Te,nrow=47)
y = starsCYG$log.light
bootstrap_result = bootstrap_lts(x,y,alphas,B=1000)
data = data.frame(h=25:46,insta_mean=bootstrap_result$means,insta_sd=bootstrap_result$sds)
ggplot(data,aes(x=h))+
  geom_point(aes(y=insta_mean,color="bootstrap"))+
  geom_line(aes(y=insta_mean,color="bootstrap"))+
  geom_errorbar(aes(ymin=insta_mean-insta_sd,ymax=insta_mean+insta_sd,color="bootstrap"),width=0.01)+
  scale_color_manual("Method",breaks=c("bootstrap"),
                     values=c("navyblue"))+
  labs(y = "Instability", x = "h", title="Instability Path")+theme_bw()+
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none", text = element_text(size=12))

