source("MCD_spectral.R")
library(mrfDepth)
library(ggplot2)
library(gridExtra)
data(glass)
glass = t(matrix(glass,nrow=750,ncol=180))
ptm <- proc.time()
result = bootstrap_insta(glass,alpha_list=seq(0.5,0.975,by=0.025),q_list=c(2,3,10),B=100)
time <- proc.time() - ptm
x = glass
x_center = apply(x,2,mean)
x_centered = sweep(x,2,x_center)
svd_result = svd(x_centered)
PCs <- x_centered%*%svd_result$v

data1 = data.frame(h=seq(0.5,0.975,by=0.025),mean_q2 = result$means[1,],sd_q2 = result$sds[1,],
                   mean_q3 = result$means[2,],sd_q3 = result$sds[2,],
                   mean_q10 = result$means[3,],sd_q10 = result$sds[3,])
plot1  = ggplot(data1,aes(x=h))+
  geom_point(aes(y=mean_q2,color="q=2"))+
  geom_point(aes(y=mean_q3,color="q=3"))+
  geom_point(aes(y=mean_q10,color="q=10"))+
  geom_line(aes(y=mean_q2,color="q=2"))+
  geom_line(aes(y=mean_q3,color="q=3"))+
  geom_line(aes(y=mean_q10,color="q=10"))+
  scale_color_manual("Number of PCs",breaks=c("q=2","q=3","q=10"),
                     values=c("navyblue","darkgreen","darkred"))+
  labs(y = "Instability", x = "h/n", title="Instability")+theme_bw()+
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none", text = element_text(size=8))

spectral <- MCD_spectral(x,q=2,alpha=0.775)
is_outlier_spectral = rep("outlier",dim(x)[1])
is_outlier_spectral[spectral$best] = "inlier"
data2 = data.frame(X1 = PCs[,1],X2=PCs[,2],X3=PCs[,3])
data2["Class"] = is_outlier_spectral
plot2 = ggplot()+ geom_point(data=data2, aes(x=X1,y=X2,color=Class),size=1) + 
  labs(y = "PC2", x = "PC1", title="SpectralMCD") + theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none", text = element_text(size=8)) +
  scale_color_manual("inlier/outlier",breaks=c("inlier","outlier"),
                     values=c("blue","red"))

spectral <- MCD_spectral(x,q=3,alpha=0.625)
is_outlier_spectral = rep("outlier",dim(x)[1])
is_outlier_spectral[spectral$best] = "inlier"
data3 = data.frame(X1 = PCs[,1],X2=PCs[,2],X3=PCs[,3])
data3["Class"] = is_outlier_spectral
plot3 = ggplot()+ geom_point(data=data3, aes(x=X1,y=X2,color=Class),size=1) + 
  labs(y = "PC2", x = "PC1", title="SpectralMCD") + theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none", text = element_text(size=8)) +
  scale_color_manual("inlier/outlier",breaks=c("inlier","outlier"),
                     values=c("blue","red"))

grid.arrange(plot1, plot2,plot3, nrow=1, ncol=3)
