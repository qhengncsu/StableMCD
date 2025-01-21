library(MASS)
library(StableMCD)
library(ggplot2)
library(gridExtra)

n = 400
p = 40

FDB <- function(x, alpha = 0.75, depth = "pro", k = 1000) 
{
  na.x <- complete.cases(x)
  if (sum(na.x) != nrow(x)) {
    x <- x[na.x, ]
    warning(paste("Observetions #:", which(na.x == 0), "where removed"))
  }
  Data <- data.matrix(x)
  n <- nrow(Data)
  p <- ncol(Data)
  if (depth == "pro") {
    pro <- depth.projection(x, x, num.directions = k)
    index11 <- order(pro, decreasing = TRUE)[1:(alpha * n)]
  }
  subset <- x[index11, ]
  hat_mu <- colMeans(subset)
  hat_sigma <- cov(subset)
  MD_c0 <- mahalanobis(x, hat_mu, hat_sigma)
  c_s <- median(MD_c0)/qchisq(0.5, p)
  sigma_raw <- c_s * hat_sigma
  MD <- mahalanobis(x, hat_mu, sigma_raw)
  trunc <- which(MD >= qchisq(0.975, p))
  x_trunc <- x[-trunc, ]
  center <- colMeans(x_trunc)
  cov <- cov(x_trunc)
  return(list(center = center, cov = cov, best = index11, raw.center = hat_mu, 
              raw.cov = sigma_raw, rew.md = MD))
}

Schmid_orthogonalization<-function(vec_u){
  n<-dim(vec_u)[1]
  s<-dim(vec_u)[2]
  beta1<-vec_u[,1]
  result<-beta1
  if(s>1){
    for(i in 2:s){
      vec<-rep(0,n)
      for(j in 1:(i-1)){
        beta0<-get(paste0("beta",j))
        vec<-vec-(sum(vec_u[,i]*beta0)/sum(beta0^2))*beta0
      }
      assign(paste0("beta",i),vec_u[,i]+vec)
      result<-cbind(result,get(paste0("beta",i)))
    }
  }
  return(result)
}

pvector1 <- function(p,dis){
  library(MASS)
  Sigma <- diag(length(p))
  p1 <- mvrnorm(n=1, p, Sigma)
  dis1 <- sqrt(sum((p1-p)^2))
  p11 <- (dis/dis1)*(p1-p)+p
  return(p11)
}

set.seed(1)
G <- matrix(0.75,p,p)
diag(G) = 1


y<-mvrnorm(n,rep(0,p),diag(p))
index_opt<-sample(1:n,0.05*n)
dis <- 5*(p^(1/4))
p11 <- rep(dis/sqrt(p),p)
y[index_opt,]<-mvrnorm(0.05*n,p11,diag(p))
x<-y%*%G

alphas = seq(0.5,0.975,by=0.025)
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
  geom_vline(data=data1, aes(xintercept=0.95),color='red')+
  scale_color_manual("Method",breaks=c("bootstrap"),
                     values=c("navyblue"))+
  labs(y = "Instability", x = "h/n", title="Setting 5 / Clustering Instability")+theme_bw()+
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none", text = element_text(size=8))

plot2 = ggplot(data1,aes(x=alphas))+
  geom_point(aes(y=wd_mean,color="bootstrap"))+
  geom_line(aes(y=wd_mean,color="bootstrap"))+
  geom_errorbar(aes(ymin=wd_mean-wd_sd,ymax=wd_mean+wd_sd,color="bootstrap"),width=0.001)+
  geom_vline(data=data1, aes(xintercept=0.95),color='red')+
  scale_color_manual("Method",breaks=c("bootstrap"),
                     values=c("navyblue"))+
  labs(y = "WD", x = "h/n", title="Setting 5 / Wasserstein Distance")+theme_bw()+
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none", text = element_text(size=8))

plot3 = ggplot(data1,aes(x=alphas))+
  geom_point(aes(y=iim,color="bootstrap"))+
  geom_line(aes(y=iim,color="bootstrap"))+
  #geom_errorbar(aes(ymin=insta_mean-insta_sd,ymax=insta_mean+insta_sd,color="bootstrap"),width=0.001)+
  geom_vline(data=data1, aes(xintercept=0.95),color='red')+
  scale_color_manual("Method",breaks=c("bootstrap"),
                     values=c("navyblue"))+
  labs(y = "IIM", x = "h/n", title="Setting 5 / Integrated Instability Metric")+theme_bw()+
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none", text = element_text(size=8))

y<-mvrnorm(n,rep(0,p),diag(p))
index_opt<-sample(1:n,0.25*n)
dis <- 50*(p^(1/4))
p11 <- rep(dis/sqrt(p),p)
y[index_opt[1:80],]<-mvrnorm(0.2*n,p11,diag(p))
vec_u<-matrix(c(rep(1,p),sample(0:(10*p),p,replace = TRUE)),byrow=F,ncol=2)
a<-Schmid_orthogonalization(vec_u)[,2]
unit_a<-a/sqrt(sum(a^2))
y[index_opt[81:100],]<-mvrnorm(0.05*n,5*sqrt(p)*unit_a,0.01^2*diag(p))
x<-y%*%G

alphas = seq(0.5,0.975,by=0.025)
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
  geom_vline(data=data1, aes(xintercept=0.75),color='red')+
  scale_color_manual("Method",breaks=c("bootstrap"),
                     values=c("navyblue"))+
  labs(y = "Instability", x = "h/n", title="Setting 6 / Clustering Instability")+theme_bw()+
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none", text = element_text(size=8))

plot5 = ggplot(data2,aes(x=alphas))+
  geom_point(aes(y=wd_mean,color="bootstrap"))+
  geom_line(aes(y=wd_mean,color="bootstrap"))+
  geom_errorbar(aes(ymin=wd_mean-wd_sd,ymax=wd_mean+wd_sd,color="bootstrap"),width=0.001)+
  geom_vline(data=data1, aes(xintercept=0.75),color='red')+
  scale_color_manual("Method",breaks=c("bootstrap"),
                     values=c("navyblue"))+
  labs(y = "WD", x = "h/n", title="Setting 6 / Wasserstein Distance")+theme_bw()+
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none", text = element_text(size=8))

plot6 = ggplot(data2,aes(x=alphas))+
  geom_point(aes(y=iim,color="bootstrap"))+
  geom_line(aes(y=iim,color="bootstrap"))+
  #geom_errorbar(aes(ymin=insta_mean-insta_sd,ymax=insta_mean+insta_sd,color="bootstrap"),width=0.001)+
  geom_vline(data=data1, aes(xintercept=0.75),color='red')+
  scale_color_manual("Method",breaks=c("bootstrap"),
                     values=c("navyblue"))+
  labs(y = "IIM", x = "h/n", title="Setting 6 / Integrated Instability Metric")+theme_bw()+
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none", text = element_text(size=8))

y<-mvrnorm(n,rep(0,p),diag(p))
index_opt<-sample(1:n,0.2*n)
dis <- 5*(p^(1/4))
for (l in index_opt[1:20]){
  p11 = pvector1(rep(0,p),dis)
  y[l,] <- mvrnorm(1,p11,diag(p))}
y[index_opt[21:80],]<-mvrnorm(0.15*n,rep(0,p),50*diag(p))
x <- y%*%G

alphas = seq(0.5,0.975,by=0.025)
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
  geom_vline(data=data1, aes(xintercept=0.8),color='red')+
  scale_color_manual("Method",breaks=c("bootstrap"),
                     values=c("navyblue"))+
  labs(y = "Instability", x = "h/n", title="Setting 7 / Clustering Instability")+theme_bw()+
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none", text = element_text(size=8))

plot8 = ggplot(data3,aes(x=alphas))+
  geom_point(aes(y=wd_mean,color="bootstrap"))+
  geom_line(aes(y=wd_mean,color="bootstrap"))+
  geom_errorbar(aes(ymin=wd_mean-wd_sd,ymax=wd_mean+wd_sd,color="bootstrap"),width=0.001)+
  geom_vline(data=data1, aes(xintercept=0.8),color='red')+
  scale_color_manual("Method",breaks=c("bootstrap"),
                     values=c("navyblue"))+
  labs(y = "WD", x = "h/n", title="Setting 7 / Wasserstein Distance")+theme_bw()+
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none", text = element_text(size=8))

plot9 = ggplot(data3,aes(x=alphas))+
  geom_point(aes(y=iim,color="bootstrap"))+
  geom_line(aes(y=iim,color="bootstrap"))+
  #geom_errorbar(aes(ymin=insta_mean-insta_sd,ymax=insta_mean+insta_sd,color="bootstrap"),width=0.001)+
  geom_vline(data=data1, aes(xintercept=0.8),color='red')+
  scale_color_manual("Method",breaks=c("bootstrap"),
                     values=c("navyblue"))+
  labs(y = "IIM", x = "h/n", title="Setting 7 / Integrated Instability Metric")+theme_bw()+
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none", text = element_text(size=8))


y<-mvrnorm(n,rep(0,p),diag(p))
index_opt<-sample(1:n,0.35*n)
dis <- 5*(p^(1/4))
for (l in index_opt[1:70]){
  p11 = pvector1(rep(0,p),dis)
  y[l,] <- mvrnorm(1,p11,diag(p))}
dis <- 50*(p^(1/4))
p11 <- rep(dis/sqrt(p),p)
y[index_opt[71:140],]<-mvrnorm(0.175*n,p11,diag(p))
x <- y%*%G

alphas = seq(0.5,0.975,by=0.025)
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
  geom_vline(data=data1, aes(xintercept=0.65),color='red')+
  scale_color_manual("Method",breaks=c("bootstrap"),
                     values=c("navyblue"))+
  labs(y = "Instability", x = "h/n", title="Setting 8 / Clustering Instability")+theme_bw()+
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none", text = element_text(size=8))

plot11 = ggplot(data4,aes(x=alphas))+
  geom_point(aes(y=wd_mean,color="bootstrap"))+
  geom_line(aes(y=wd_mean,color="bootstrap"))+
  geom_errorbar(aes(ymin=wd_mean-wd_sd,ymax=wd_mean+wd_sd,color="bootstrap"),width=0.001)+
  geom_vline(data=data1, aes(xintercept=0.65),color='red')+
  scale_color_manual("Method",breaks=c("bootstrap"),
                     values=c("navyblue"))+
  labs(y = "WD", x = "h/n", title="Setting 8 / Wasserstein Distance")+theme_bw()+
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none", text = element_text(size=8))

plot12 = ggplot(data4,aes(x=alphas))+
  geom_point(aes(y=iim,color="bootstrap"))+
  geom_line(aes(y=iim,color="bootstrap"))+
  #geom_errorbar(aes(ymin=insta_mean-insta_sd,ymax=insta_mean+insta_sd,color="bootstrap"),width=0.001)+
  geom_vline(data=data1, aes(xintercept=0.65),color='red')+
  scale_color_manual("Method",breaks=c("bootstrap"),
                     values=c("navyblue"))+
  labs(y = "IIM", x = "h/n", title="Setting 8 / Integrated Instability Metric")+theme_bw()+
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none", text = element_text(size=8))


grid.arrange(plot1, plot2, plot3, plot4, plot5, plot6, plot7, plot8, plot9, plot10, plot11, plot12, nrow=4, ncol=3)
