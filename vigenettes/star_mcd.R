library(StableMCD)
library(DetMCD)
library(ggplot2)
library(gridExtra)
library(robustbase)
data(starsCYG)
x = as.matrix(starsCYG)
alphas = seq(25,46)/47+0.001

set.seed(1)
ptm <- proc.time()
bootstrap_result = bootstrap_mcd(x,alphas,B=100,classifier="MD")
time <- proc.time() - ptm

data = data.frame(alpha = alphas,insta_mean = bootstrap_result$insta_means,
                  insta_sd =bootstrap_result$insta_sds,
                  wd_mean = bootstrap_result$wd_means,
                  wd_sd = bootstrap_result$wd_sds,
                  iim = bootstrap_result$iim, h = seq(25,46))

plot1 = ggplot(data,aes(x=h))+
        geom_point(aes(y=insta_mean,color="bootstrap"))+
        geom_line(aes(y=insta_mean,color="bootstrap"))+
        geom_errorbar(aes(ymin=insta_mean-insta_sd,ymax=insta_mean+insta_sd,color="bootstrap"),width=0.01)+
        scale_color_manual("Method",breaks=c("bootstrap"),
                     values=c("navyblue"))+
        labs(y = "Instability", x = "h", title="Clustering Instability on Star Data")+theme_bw()+
        theme(plot.title = element_text(hjust = 0.5), legend.position = "none", text = element_text(size=8))

plot2 = ggplot(data,aes(x=h))+
  geom_point(aes(y=iim,color="bootstrap"))+
  geom_line(aes(y=iim,color="bootstrap"))+
  #geom_errorbar(aes(ymin=insta_mean-insta_sd,ymax=insta_mean+insta_sd,color="bootstrap"),width=0.01)+
  scale_color_manual("Method",breaks=c("bootstrap"),
                     values=c("navyblue"))+
  labs(y = "IIM", x = "h", title="Integrated Instability Metric on Star Data")+theme_bw()+
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none", text = element_text(size=8))

result = mcd(x,40/47+0.01)
is_outlier = rep("outlier",dim(x)[1])
is_outlier[result$index] = "inlier"
data3 = data.frame(X1 = x[,1],X2=x[,2])
data3["Class"] = is_outlier
plot3 = ggplot()+ geom_point(data=data3, aes(x=X1,y=X2,color=Class),size=3) + 
  labs(y = "Log Light", x = "Log Temperature", title="Inlier/Outlier Map for Star Data") + theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none", text = element_text(size=8)) +
  scale_color_manual("inlier/outlier",breaks=c("inlier","outlier"),
                     values=c("blue","red"))

grid.arrange(plot1, plot2, plot3, nrow=1, ncol=3)
