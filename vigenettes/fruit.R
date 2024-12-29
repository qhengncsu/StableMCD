library(ggplot2)
library(gridExtra)
library(StableMCD)
library(rrcov)

data(fruit)
x = data.matrix(fruit[,2:257])

ptm <- proc.time()
result = bootstrap_mcd(x,seq(0.5,0.975,by=0.025),B=50,classifier="MD")
time <- proc.time() - ptm
alphas = seq(0.5,0.975,by=0.025)
data = data.frame(alpha = alphas,insta_mean = result$insta_means,
                  insta_sd =result$insta_sds,
                  wd_mean = result$wd_means,
                  wd_sd = result$wd_sds,
                  final_score = result$final_score)

plot1 = ggplot(data,aes(x=alphas))+
  geom_point(aes(y=insta_mean,color="bootstrap"))+
  geom_line(aes(y=insta_mean,color="bootstrap"))+
  geom_errorbar(aes(ymin=insta_mean-insta_sd,ymax=insta_mean+insta_sd,color="bootstrap"),width=0.01)+
  scale_color_manual("Method",breaks=c("bootstrap"),
                     values=c("navyblue"))+
  labs(y = "Instability", x = "alpha", title="Instability on Fruit Data")+theme_bw()+
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none", text = element_text(size=8))

plot2 = ggplot(data,aes(x=alphas))+
  geom_point(aes(y=final_score,color="bootstrap"))+
  geom_line(aes(y=final_score,color="bootstrap"))+
  #geom_errorbar(aes(ymin=insta_mean-insta_sd,ymax=insta_mean+insta_sd,color="bootstrap"),width=0.01)+
  scale_color_manual("Method",breaks=c("bootstrap"),
                     values=c("navyblue"))+
  labs(y = "Final Score", x = "alpha", title="Final Score on Fruit Data")+theme_bw()+
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none", text = element_text(size=8))

grid.arrange(plot1, plot2, nrow=1, ncol=2)
