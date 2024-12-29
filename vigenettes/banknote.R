library(DetMCD)
library(mclust)
library(StableMCD)
library(ddalpha)
library(ggplot2)
library(gridExtra)

set.seed(12)
data(banknote)
x = as.matrix(banknote[101:200,2:7])
alphas = seq(0.5,0.99,0.01)
h = alphas*100
ptm <- proc.time()
bootstrap_result = bootstrap_mcd(x,seq(0.5,0.99,0.01),B=50,classifier='MD')
time <- proc.time() - ptm

data = data.frame(alpha = alphas,insta_mean = bootstrap_result$insta_means,
                  insta_sd =bootstrap_result$insta_sds,
                  wd_mean = bootstrap_result$wd_means,
                  wd_sd = bootstrap_result$wd_sds,
                  final_score = bootstrap_result$final_score)

plot1 = ggplot(data,aes(x=h))+
  geom_point(aes(y=insta_mean,color="bootstrap"))+
  geom_line(aes(y=insta_mean,color="bootstrap"))+
  geom_errorbar(aes(ymin=insta_mean-insta_sd,ymax=insta_mean+insta_sd,color="bootstrap"),width=0.01)+
  scale_color_manual("Method",breaks=c("bootstrap"),
                     values=c("navyblue"))+
  labs(y = "Instability", x = "h", title="Instability on Banknote Data")+theme_bw()+
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none", text = element_text(size=8))

plot2 = ggplot(data,aes(x=h))+
  geom_point(aes(y=final_score,color="bootstrap"))+
  geom_line(aes(y=final_score,color="bootstrap"))+
  #geom_errorbar(aes(ymin=insta_mean-insta_sd,ymax=insta_mean+insta_sd,color="bootstrap"),width=0.01)+
  scale_color_manual("Method",breaks=c("bootstrap"),
                     values=c("navyblue"))+
  labs(y = "Final Score", x = "h", title="Final Score on Banknote Data")+theme_bw()+
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none", text = element_text(size=8))


# subset = mcd(x,0.75)
# mds = mahalanobis(x, subset$muhat, subset$Sigmahat)
# mds = sort(mds)
# data2 = data.frame(index = 51:dim(x)[1],mds=mds[51:100])
# plot2 = ggplot(data2, aes(x=index, y=mds)) + 
#   geom_bar(stat = "identity",color='lightblue',fill='lightblue')+
#   labs(title = "Mahalanobis Distances", x = "Observation Index", y = "MD") +theme_bw()+
#   geom_vline(data=data2, aes(xintercept=bootstrap_result$best_alpha*dim(x)[1]),color='red')+
#   theme(plot.title = element_text(hjust = 0.5),text = element_text(size=12))

grid.arrange(plot1, plot2, nrow=1, ncol=2)
