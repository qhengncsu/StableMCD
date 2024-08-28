library(robustbase)
library(ggplot2)
library(gridExtra)
library(dplyr)
library(fsdaR)
library(StableMCD)
data("starsCYG")
alphas = rep(0,22)

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
alphas = seq(25,46)/47+0.001

x = matrix(starsCYG$log.Te,nrow=47)
y = starsCYG$log.light
bootstrap_result = bootstrap_lts(x,y,alphas,500)
normality_result = ltsReg_find_h(y,x,alphas)
init = lts(x,y,0.5)$indexes[[1]]+1
fs_result <- fsreg(y~x,control=FSR_control(plot=F),trace=FALSE,bsb=init)



lts_normality = lts(x,y,normality_result$best_alpha)
is_outlier_normality = rep("outlier",dim(x)[1])
is_outlier_normality[lts_normality$indexes[[1]]+1] = "inlier"
data1 = data.frame(logTemp=starsCYG$log.Te,logLight=y, Class=is_outlier_normality)
plot1 = ggplot()+ geom_point(data=data1, aes(x=logTemp,y=y,color=is_outlier_normality)) + 
  labs(y = "Log Light", x = "Log Temperature", title="Normality Test") + theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none", text = element_text(size=12)) +
  scale_color_manual(breaks=c("inlier","outlier"),
                     values=c("blue","red"))

is_outlier_fs = rep("inlier",dim(x)[1])
is_outlier_fs[fs_result$outliers] = "outlier"
data2 = data.frame(logTemp=starsCYG$log.Te,logLight=y, Class=is_outlier_fs)
plot2 = ggplot()+ geom_point(data=data2, aes(x=logTemp,y=y,color=is_outlier_fs)) + 
  labs(y = "Log Light", x = "Log Temperature", title="Forward Search") + theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none", text = element_text(size=12)) +
  scale_color_manual(breaks=c("inlier","outlier"),
                     values=c("blue","red"))

data3 = data.frame(h=25:46,insta_mean=bootstrap_result$insta_means,insta_sd=bootstrap_result$insta_sds)
plot3 = ggplot(data3,aes(x=h))+
  geom_point(aes(y=insta_mean,color="bootstrap"))+
  geom_line(aes(y=insta_mean,color="bootstrap"))+
  geom_errorbar(aes(ymin=insta_mean-insta_sd,ymax=insta_mean+insta_sd,color="bootstrap"),width=0.01)+
  scale_color_manual("Method",breaks=c("bootstrap"),
                     values=c("navyblue"))+
  labs(y = "Instability", x = "h", title="Instability Path")+theme_bw()+
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none", text = element_text(size=12))

lts_bootstrap = lts(x,y,alpha=bootstrap_result$best_alpha)
is_outlier_bootstrap = rep("outlier",dim(x)[1])
is_outlier_bootstrap[lts_bootstrap$indexes[[1]]+1] = "inlier"
data4 = data.frame(logTemp=starsCYG$log.Te,logLight=y, Class=is_outlier_bootstrap)
plot4 = ggplot()+ geom_point(data=data4, aes(x=logTemp,y=y,color=is_outlier_bootstrap)) + 
  labs(y = "Log Light", x = "Log Temperature", title="Bootstrap") + theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none", text = element_text(size=12)) +
  scale_color_manual(breaks=c("inlier","outlier"),
                     values=c("blue","red"))



grid.arrange(plot1, plot2, plot3, plot4, nrow=2, ncol=2)
