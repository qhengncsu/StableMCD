library(MASS)
library(StableMCD)
library(ggplot2)
library(gridExtra)

n = 400
p = 40

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

y<-mvrnorm(n,rep(0,p),diag(p))
index_opt<-sample(1:n,0.35*n)
dis <- 50*(p^(1/4))
for (l in index_opt[1:70]){
  p11 = pvector1(rep(0,p),dis)
  y[l,] <- mvrnorm(1,p11,diag(p))}
dis <- 5*(p^(1/4))
p11 <- rep(dis/sqrt(p),p)
y[index_opt[71:140],]<-mvrnorm(0.175*n,p11,diag(p))
x <- y%*%G

alphas = seq(0.5,0.975,by=0.025)
bootstrap_result = bootstrap_mcd(x,seq(0.5,0.975,by=0.025),50,classifier = "MD",sd_ratio=1.5)
data1 = data.frame(alpha = alphas,insta_mean = bootstrap_result$insta_means,
                   insta_sd =bootstrap_result$insta_sds,
                   wd_mean = bootstrap_result$wd_means,
                   wd_sd = bootstrap_result$wd_sds,
                   iim = bootstrap_result$iim)

plot1 = ggplot(data1,aes(x=alphas))+
  geom_point(aes(y=iim,color="bootstrap"))+
  geom_line(aes(y=iim,color="bootstrap"))+
  #geom_errorbar(aes(ymin=insta_mean-insta_sd,ymax=insta_mean+insta_sd,color="bootstrap"),width=0.001)+
  geom_vline(data=data1, aes(xintercept=0.65),color='red')+
  scale_color_manual("Method",breaks=c("bootstrap"),
                     values=c("navyblue"))+
  labs(y = "IIM", x = "h/n", title="Setting 8 / Integrated Instability Metric usng lambda=1.5")+theme_bw()+
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none", text = element_text(size=10))

trimmed = mcd(x,0.8)
x_trimmed = x[trimmed$index,]
alphas = seq(0.5,0.975,by=0.025)
bootstrap_result = bootstrap_mcd(x_trimmed,seq(0.5,0.975,by=0.025),50,classifier = "MD",sd_ratio=3)
data2 = data.frame(alpha = alphas,insta_mean = bootstrap_result$insta_means,
                   insta_sd =bootstrap_result$insta_sds,
                   wd_mean = bootstrap_result$wd_means,
                   wd_sd = bootstrap_result$wd_sds,
                   iim = bootstrap_result$iim)

plot2 = ggplot(data2,aes(x=alphas))+
  geom_point(aes(y=iim,color="bootstrap"))+
  geom_line(aes(y=iim,color="bootstrap"))+
  #geom_errorbar(aes(ymin=insta_mean-insta_sd,ymax=insta_mean+insta_sd,color="bootstrap"),width=0.001)+
  geom_vline(data=data1, aes(xintercept=0.8125),color='red')+
  scale_color_manual("Method",breaks=c("bootstrap"),
                     values=c("navyblue"))+
  labs(y = "IIM", x = "h/n", title="Setting 8 / Integrated Instability Metric after Trimming")+theme_bw()+
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none", text = element_text(size=10))

grid.arrange(plot1, plot2, nrow=1, ncol=2)
