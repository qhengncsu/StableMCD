library(MASS)
library(StableMCD)
library(ggplot2)
library(gridExtra)

set.seed(1234)

alphas = seq(0.5,0.975,by=0.025)

x1 = rnorm(1000)
x2 = rnorm(1000)
x1[1:100] = rnorm(100,mean=5)
x2[1:100] = rnorm(100,mean=5)
#x1[101:150] = rnorm(50,mean=1000)
#x2[101:150] = rnorm(50,mean=1000)
x = cbind(x1,x2)

bootstrap_result = bootstrap_mcd(x,seq(0.5,0.975,by=0.025),50,classifier = "MD")
data1 = data.frame(alpha = alphas,insta_mean = bootstrap_result$insta_means,
                   insta_sd =bootstrap_result$insta_sds,
                   wd_mean = bootstrap_result$wd_means,
                   wd_sd = bootstrap_result$wd_sds,
                   iim = bootstrap_result$iim)

plot1 = ggplot(data1,aes(x=alphas))+
  geom_point(aes(y=insta_mean,color="bootstrap"))+
  geom_line(aes(y=insta_mean,color="bootstrap"))+
  geom_errorbar(aes(ymin=insta_mean-insta_sd,ymax=insta_mean+insta_sd,color="bootstrap"),width=0.001)+
  geom_vline(data=data1, aes(xintercept=0.9),color='red')+
  scale_color_manual("Method",breaks=c("bootstrap"),
                     values=c("navyblue"))+
  labs(y = "Instability", x = "h/n", title="Setting 1 / Clustering Instability")+theme_bw()+
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none", text = element_text(size=8))

plot2 = ggplot(data1,aes(x=alphas))+
  geom_point(aes(y=wd_mean,color="bootstrap"))+
  geom_line(aes(y=wd_mean,color="bootstrap"))+
  geom_errorbar(aes(ymin=wd_mean-wd_sd,ymax=wd_mean+wd_sd,color="bootstrap"),width=0.001)+
  geom_vline(data=data1, aes(xintercept=0.9),color='red')+
  scale_color_manual("Method",breaks=c("bootstrap"),
                     values=c("navyblue"))+
  labs(y = "WD", x = "h/n", title="Setting 1 / Wasserstein Distance")+theme_bw()+
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none", text = element_text(size=8))

plot3 = ggplot(data1,aes(x=alphas))+
  geom_point(aes(y=iim,color="bootstrap"))+
  geom_line(aes(y=iim,color="bootstrap"))+
  #geom_errorbar(aes(ymin=insta_mean-insta_sd,ymax=insta_mean+insta_sd,color="bootstrap"),width=0.001)+
  geom_vline(data=data1, aes(xintercept=0.9),color='red')+
  scale_color_manual("Method",breaks=c("bootstrap"),
                     values=c("navyblue"))+
  labs(y = "IIM", x = "h/n", title="Setting 1 / Integrated Instability Metric")+theme_bw()+
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none", text = element_text(size=8))

x1 = rnorm(1000)
x2 = rnorm(1000)
x1[1:100] = rnorm(100,mean=5)
x2[1:100] = rnorm(100,mean=5)
x1[101:150] = rnorm(50,mean=1000)
x2[101:150] = rnorm(50,mean=1000)
x = cbind(x1,x2)

bootstrap_result = bootstrap_mcd(x,seq(0.5,0.975,by=0.025),50,classifier = "MD")
data2 = data.frame(alpha = alphas,insta_mean = bootstrap_result$insta_means,
                   insta_sd =bootstrap_result$insta_sds,
                   wd_mean = bootstrap_result$wd_means,
                   wd_sd = bootstrap_result$wd_sds,
                   iim = bootstrap_result$iim)

plot4 = ggplot(data2,aes(x=alphas))+
  geom_point(aes(y=insta_mean,color="bootstrap"))+
  geom_line(aes(y=insta_mean,color="bootstrap"))+
  geom_errorbar(aes(ymin=insta_mean-insta_sd,ymax=insta_mean+insta_sd,color="bootstrap"),width=0.001)+
  geom_vline(data=data1, aes(xintercept=0.85),color='red')+
  scale_color_manual("Method",breaks=c("bootstrap"),
                     values=c("navyblue"))+
  labs(y = "Instability", x = "h/n", title="Setting 2 / Clustering Instability")+theme_bw()+
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none", text = element_text(size=8))

plot5 = ggplot(data2,aes(x=alphas))+
  geom_point(aes(y=wd_mean,color="bootstrap"))+
  geom_line(aes(y=wd_mean,color="bootstrap"))+
  geom_errorbar(aes(ymin=wd_mean-wd_sd,ymax=wd_mean+wd_sd,color="bootstrap"),width=0.001)+
  geom_vline(data=data1, aes(xintercept=0.85),color='red')+
  scale_color_manual("Method",breaks=c("bootstrap"),
                     values=c("navyblue"))+
  labs(y = "WD", x = "h/n", title="Setting 2 / Wasserstein Distance")+theme_bw()+
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none", text = element_text(size=8))

plot6 = ggplot(data2,aes(x=alphas))+
  geom_point(aes(y=iim,color="bootstrap"))+
  geom_line(aes(y=iim,color="bootstrap"))+
  #geom_errorbar(aes(ymin=insta_mean-insta_sd,ymax=insta_mean+insta_sd,color="bootstrap"),width=0.001)+
  geom_vline(data=data1, aes(xintercept=0.85),color='red')+
  scale_color_manual("Method",breaks=c("bootstrap"),
                     values=c("navyblue"))+
  labs(y = "IIM", x = "h/n", title="Setting 2 / Integrated Instability Metric")+theme_bw()+
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none", text = element_text(size=8))

x1 = rnorm(1000)
x2 = rnorm(1000)
x1[1:100] = rnorm(100,mean=5)
x2[1:100] = rnorm(100,mean=5)
x1[101:300] = seq(10,2000,by=10)
x2[101:300] = seq(10,2000,by=10)
x = cbind(x1,x2)

bootstrap_result = bootstrap_mcd(x,seq(0.5,0.975,by=0.025),50,classifier = "MD")
data3 = data.frame(alpha = alphas,insta_mean = bootstrap_result$insta_means,
                   insta_sd =bootstrap_result$insta_sds,
                   wd_mean = bootstrap_result$wd_means,
                   wd_sd = bootstrap_result$wd_sds,
                   iim = bootstrap_result$iim)

plot7 = ggplot(data3,aes(x=alphas))+
  geom_point(aes(y=insta_mean,color="bootstrap"))+
  geom_line(aes(y=insta_mean,color="bootstrap"))+
  geom_errorbar(aes(ymin=insta_mean-insta_sd,ymax=insta_mean+insta_sd,color="bootstrap"),width=0.001)+
  geom_vline(data=data1, aes(xintercept=0.7),color='red')+
  scale_color_manual("Method",breaks=c("bootstrap"),
                     values=c("navyblue"))+
  labs(y = "Instability", x = "h/n", title="Setting 3 / Clustering Instability")+theme_bw()+
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none", text = element_text(size=8))

plot8 = ggplot(data3,aes(x=alphas))+
  geom_point(aes(y=wd_mean,color="bootstrap"))+
  geom_line(aes(y=wd_mean,color="bootstrap"))+
  geom_errorbar(aes(ymin=wd_mean-wd_sd,ymax=wd_mean+wd_sd,color="bootstrap"),width=0.001)+
  geom_vline(data=data1, aes(xintercept=0.7),color='red')+
  scale_color_manual("Method",breaks=c("bootstrap"),
                     values=c("navyblue"))+
  labs(y = "WD", x = "h/n", title="Setting 3 / Wasserstein Distance")+theme_bw()+
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none", text = element_text(size=8))

plot9 = ggplot(data3,aes(x=alphas))+
  geom_point(aes(y=iim,color="bootstrap"))+
  geom_line(aes(y=iim,color="bootstrap"))+
  #geom_errorbar(aes(ymin=insta_mean-insta_sd,ymax=insta_mean+insta_sd,color="bootstrap"),width=0.001)+
  geom_vline(data=data1, aes(xintercept=0.7),color='red')+
  scale_color_manual("Method",breaks=c("bootstrap"),
                     values=c("navyblue"))+
  labs(y = "IIM", x = "h/n", title="Setting 3 / Integrated Instability Metric")+theme_bw()+
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none", text = element_text(size=8))



x1 = rnorm(1000)
x2 = rnorm(1000)
x = cbind(x1,x2)

bootstrap_result = bootstrap_mcd(x,seq(0.5,0.975,by=0.025),50,classifier = "MD")
data4 = data.frame(alpha = alphas,insta_mean = bootstrap_result$insta_means,
                   insta_sd =bootstrap_result$insta_sds,
                   wd_mean = bootstrap_result$wd_means,
                   wd_sd = bootstrap_result$wd_sds,
                   iim = bootstrap_result$iim)

plot10 = ggplot(data4,aes(x=alphas))+
  geom_point(aes(y=insta_mean,color="bootstrap"))+
  geom_line(aes(y=insta_mean,color="bootstrap"))+
  geom_errorbar(aes(ymin=insta_mean-insta_sd,ymax=insta_mean+insta_sd,color="bootstrap"),width=0.001)+
  #geom_vline(data=data1, aes(xintercept=0.7),color='red')+
  scale_color_manual("Method",breaks=c("bootstrap"),
                     values=c("navyblue"))+
  labs(y = "Instability", x = "h/n", title="Setting 4 / Clustering Instability")+theme_bw()+
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none", text = element_text(size=8))

plot11 = ggplot(data4,aes(x=alphas))+
  geom_point(aes(y=wd_mean,color="bootstrap"))+
  geom_line(aes(y=wd_mean,color="bootstrap"))+
  geom_errorbar(aes(ymin=wd_mean-wd_sd,ymax=wd_mean+wd_sd,color="bootstrap"),width=0.001)+
  #geom_vline(data=data1, aes(xintercept=0.7),color='red')+
  scale_color_manual("Method",breaks=c("bootstrap"),
                     values=c("navyblue"))+
  labs(y = "WD", x = "h/n", title="Setting 4 / Wasserstein Distance")+theme_bw()+
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none", text = element_text(size=8))

plot12 = ggplot(data4,aes(x=alphas))+
  geom_point(aes(y=iim,color="bootstrap"))+
  geom_line(aes(y=iim,color="bootstrap"))+
  #geom_errorbar(aes(ymin=insta_mean-insta_sd,ymax=insta_mean+insta_sd,color="bootstrap"),width=0.001)+
  #geom_vline(data=data1, aes(xintercept=0.7),color='red')+
  scale_color_manual("Method",breaks=c("bootstrap"),
                     values=c("navyblue"))+
  labs(y = "IIM", x = "h/n", title="Setting 4 / Integrated Instability Metric")+theme_bw()+
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none", text = element_text(size=8))

grid.arrange(plot1, plot2, plot3, plot4, plot5, plot6, plot7, plot8, plot9, plot10, plot11, plot12, nrow=4, ncol=3)
