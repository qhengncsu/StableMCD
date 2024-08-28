library(robustbase)
library(ggplot2)
library(gridExtra)
library(dplyr)
library(fsdaR)
library(StableMCD)
load('glass.rdata')
ltsReg_find_h	<- function(y,X,alphas,do.set.seed=TRUE){
  n	<- length(y)
  p	<- ncol(X)
  if(is.null(p))	p	<- 1
  p	<- p+1		# count intercept
  #h	<- floor((n+p+1)/2):n
  #names	<- list(h,c("alpha","h",as.character(1:p),
  #"RSS","sigma","cum3","cum4","|cum4|","T","IC",
  #"sigma.df","cum3.df","cum4.df","|cum4|.df","T.df"))
  #m.results	<- matrix(nrow=length(h),ncol=length(names[[2]]),dimnames=names)
  ###############
  #	loop over h
  #	compute LTS using ltsReg
  #	compute derived statistics for each h
  thieles <- rep(0,length(alphas))
  for(i in 1:length(alphas))
  {
    alpha	<- alphas[i]
    h <- floor(n*alpha)
    if(do.set.seed)	set.seed(0)
    lts_result <- lts(X,y,alpha)
    beta = lts_result$betas[[1]]
    if(p>2)
    {	residuals 	<- y- beta[1] - X %*% beta[2:p]	}
    else
    {	residuals 	<- y- beta[1] - X  *  beta[2]		}
    #if(h==n)	lts$best	<- 1:n		# ltsReg returns null	
    good.residuals	<- residuals[lts_result$indexes[[1]]+1]
    RSS <- sum(good.residuals^2)
    #	statistics normalised by sample size h
    sigma2	<- RSS/h					
    k3	<- h*sum(good.residuals^3)/RSS^(3/2)/sqrt(6)
    k4	<- sqrt(h)*(h*sum(good.residuals^4)/RSS^2 - 3)/sqrt(24)	
    Thiele	<- k3^2+k4^2
    #	statistics normalised by degrees of freedom h-p not counting intercept
    sigma2.df	<- RSS/(h-p)
    k3.df	<- (h-p)*sum(good.residuals^3)/RSS^(3/2)/sqrt(6)
    k4.df	<- sqrt(h-p)*((h-p)*sum(good.residuals^4)/RSS^2 - 3)/sqrt(24)	
    Thiele.df	<- k3.df^2+k4.df^2
    #	information criteria
    IC	<- log(sigma2)+2*log(log(h))*(1-h/n)
    #m.results[i,]	<- c(alpha,lts$quan,lts$raw.coefficients,
    #RSS,sqrt(sigma2),k3,k4,abs(k4),Thiele,IC,
    #sqrt(sigma2.df),k3.df,k4.df,abs(k4.df),Thiele.df)
    thieles[i] <- Thiele
    
  }
  best_alpha = alphas[which.min(thieles)]
  return(list(best_alpha = best_alpha, normality = thieles))
}
x_original <- x
x = as.matrix(x[,c(25,106,230)])
y = y$PbO
init = lts(x,y,0.5)$indexes[[1]]+1
fs_result <- fsreg(y~x,control=FSR_control(plot=F),trace=FALSE,bsb=init)
alphas = seq(0.6,0.99,by=0.01)
normality_result = ltsReg_find_h(y,x,alphas)
bootstrap_result = bootstrap_lts(x,y,alphas,100)
xy = cbind(x,y)
xy_center = apply(xy,2,mean)
xy_centered = sweep(xy,2,xy_center)
svd_result = svd(xy_centered)
PCs <- data.frame(xy_centered%*%svd_result$v[,1:2])

lts_normality = lts(x,y,normality_result$best_alpha)
is_outlier_normality = rep("outlier",dim(x)[1])
is_outlier_normality[lts_normality$indexes[[1]]+1] = "inlier"
residual_normality = y - lts_normality$betas[[1]][1] - x%*%lts_normality$betas[[1]][2:4]
highlight = rep(4,180)
highlight[lts_normality$best] = 1
data1 = data.frame(obs=seq(1,180,by=1),residual = residual_normality)
plot1 <- ggplot(data = data1, mapping = aes(x = obs, y = residual)) +
  geom_hline(aes(yintercept = 0)) +
  geom_segment(mapping = aes(xend = obs, yend = 0),fill=highlight,linewidth=0.75)+
  labs(y = "Residual", x = "Observation", title="Normality Test")+theme_bw()+
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none", text = element_text(size=12))

# data1 = data.frame(PC1=PCs$X1,PC2=PCs$X2, Class=is_outlier_normality,X1=x[,1],X2=x[,2],X3=x[,3],y=y)
# plot1 = ggplot()+ geom_point(data=data1[is_outlier_normality=="inlier",], aes(x=X3,y=y,color="inlier"),size=1) + 
#   geom_point(data=data1[is_outlier_normality=="outlier",], aes(x=X3,y=y,color="outlier"),size=1)+
#   labs(y = "PC2", x = "PC1", title="Normality Test") + theme_bw() + 
#   theme(plot.title = element_text(hjust = 0.5), legend.position = "none", text = element_text(size=12)) +
#   scale_color_manual(breaks=c("outlier","inlier"),
#                      values=c("red","blue"))

is_outlier_fs = rep("inlier",dim(x)[1])
is_outlier_fs[fs_result$outliers] = "outlier"
lts_fs = lts(x,y,(180-length(fs_result$outliers))/180)
is_outlier_fs = rep("outlier",dim(x)[1])
residual_fs = y - lts_fs$betas[[1]][1] - x%*%lts_fs$betas[[1]][2:4]

highlight = rep(4,180)
highlight[lts_fs$indexes[[1]]+1] = 1
highlight[c(93,163,172)] = 2
data2 = data.frame(obs=seq(1,180,by=1),residual = residual_fs,highlight=highlight)
plot2 <- ggplot(data = data2, mapping = aes(x = obs, y = residual)) +
  geom_hline(aes(yintercept = 0)) +
  geom_segment(mapping = aes(xend = obs, yend = 0),color=highlight,linewidth=0.75)+
  labs(y = "Residual", x = "Observation", title="Forward Search")+theme_bw()+
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none", text = element_text(size=12))
  
# data2 = data.frame(PC1=PCs$X1,PC2=PCs$X2, Class=is_outlier_fs,X1=x[,1],X2=x[,2],X3=x[,3],y=y)
# plot2 = ggplot()+ geom_point(data=data1[is_outlier_fs=="inlier",], aes(x=X3,y=y,color="inlier"),size=1) + 
#   geom_point(data=data1[is_outlier_fs=="outlier",], aes(x=X3,y=y,color="outlier"),size=1)+
#   labs(y = "PC2", x = "PC1", title="Forward Search") + theme_bw() + 
#   theme(plot.title = element_text(hjust = 0.5), legend.position = "none", text = element_text(size=12)) +
#   scale_color_manual(breaks=c("outlier","inlier"),
#                      values=c("red","blue"))

data3 = data.frame(alpha=alphas,insta_mean=bootstrap_result$insta_means,insta_sd=bootstrap_result$insta_sds)
plot3 = ggplot(data3,aes(x=alphas))+
  geom_point(aes(y=insta_mean,color="bootstrap"))+
  geom_line(aes(y=insta_mean,color="bootstrap"))+
  geom_errorbar(aes(ymin=insta_mean-insta_sd,ymax=insta_mean+insta_sd,color="bootstrap"),width=0.001)+
  scale_color_manual("Method",breaks=c("bootstrap"),
                     values=c("navyblue"))+
  labs(y = "Instability", x = "h/n", title="Instability Path")+theme_bw()+
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none", text = element_text(size=12))

lts_bootstrap = lts(x,y,bootstrap_result$best_alpha)
is_outlier_bootstrap = rep("outlier",dim(x)[1])
is_outlier_bootstrap[lts_bootstrap$indexes[[1]]+1] = "inlier"
residual_bootstrap = y - lts_bootstrap$betas[[1]][1] - x%*%lts_bootstrap$betas[[1]][2:4]
highlight = rep(4,180)
highlight[lts_bootstrap$indexes[[1]]+1] = 1
data4 = data.frame(obs=seq(1,180,by=1),residual = residual_bootstrap,highlight=highlight)

plot4 <- ggplot(data = data4, mapping = aes(x = obs, y = residual)) +
  geom_hline(aes(yintercept = 0)) +
  geom_segment(mapping = aes(xend = obs, yend = 0),color=highlight,linewidth=0.75)+
  labs(y = "Residual", x = "Observation", title="Bootstrap")+theme_bw()+
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none", text = element_text(size=12))
# data4 = data.frame(PC1=PCs$X1,PC2=PCs$X2, Class=is_outlier_bootstrap)
# plot4 = ggplot()+ geom_point(data=data1[is_outlier_bootstrap=="inlier",], aes(x=PC1,y=PC2,color="inlier"),size=1) + 
#   geom_point(data=data1[is_outlier_bootstrap=="outlier",], aes(x=PC1,y=PC2,color="outlier"),size=1)+
#   labs(y = "PC2", x = "PC1", title="Bootstrap") + theme_bw() + 
#   theme(plot.title = element_text(hjust = 0.5), legend.position = "none", text = element_text(size=12)) +
#   scale_color_manual(breaks=c("outlier","inlier"),
#                      values=c("red","blue"))



grid.arrange(plot1, plot2, plot3, plot4, nrow=2, ncol=2)
