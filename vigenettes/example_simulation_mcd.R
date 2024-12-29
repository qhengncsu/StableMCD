library(MASS)
library(ggplot2)
library(StableMCD)
n = 400
p = 40
outlier_ratio = 0.25
outlier_type = "Random"

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
r = 5
G <- matrix(0.75,p,p)
diag(G) = 1
reps = 1


y<-mvrnorm(n,rep(0,p),diag(p))
index_opt<-sample(1:n,outlier_ratio*n)
if(outlier_type=="Point"){
  vec_u<-matrix(c(rep(1,p),sample(0:(10*p),p,replace = TRUE)),byrow=F,ncol=2)
  a<-Schmid_orthogonalization(vec_u)[,2]
  unit_a<-a/sqrt(sum(a^2))
  y[index_opt,]<-mvrnorm(outlier_ratio*n,r*sqrt(p)*unit_a,0.01^2*diag(p))
}else if(outlier_type=="Cluster"){
  dis <- r*(p^(1/4))
  p11 <- rep(dis/sqrt(p),p)
  y[index_opt[1:10],]<-mvrnorm(0.1*outlier_ratio*n,p11,diag(p))
  y[index_opt[11:100],]<-mvrnorm(0.9*outlier_ratio*n,10*p11,diag(p))
}else if(outlier_type=="Random"){
  dis <- r*(p^(1/4))
  p11 <- rep(dis,p)
  for (l in index_opt[1:10]){
    p11=pvector1(rep(0,p),dis)
    y[l,]<-mvrnorm(1,p11,diag(p))}
  for (l in 11:100){
    p11=pvector1(rep(0,p),dis)
    y[index_opt[l],]<-mvrnorm(1,(l-10)*p11,diag(p))}
}else if(outlier_type=="Radial"){
  y[index_opt,]<-mvrnorm(outlier_ratio*n,rep(0,p),5*diag(p))
}
x<-y%*%G

ptm <- proc.time()
result = bootstrap_mcd(x,seq(0.7,0.8,by=0.005),B=100,classifier='MD')
print(result$best_alpha)
time <- proc.time() - ptm

#result = mcd(x,0.75)

data1 = data.frame(h = seq(0.7,0.8,by=0.005)*dim(x)[1],insta=as.vector(result$means),sd = as.vector(result$sds))
ggplot(data1, aes(x=h, y=insta)) + 
  geom_line(colour="navyblue") +
  geom_point(colour="navyblue")+
  geom_errorbar(aes(ymin=insta-sd,ymax=insta+sd),width=0.01)+
  labs(y = "Instability", x = "h", title="Instability")+theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),text = element_text(size=12))


set.seed(1234)
n = 1000
x1 = rnorm(n)
x2 = rnorm(n)

x = cbind(x1, x2)#

x[1:50, ] = rep(c(-10, 10), each = 5) + rnorm(50)
x[51:60, ] = 200 + rnorm(10)
plot(x)
plot(mahalanobis(x, rep(0, 2), diag(2)))


dvals = seq(0.9,0.99,by=0.01)
result = bootstrap_mcd(x,dvals,B=100,classifier="MD")
cbind(dvals, result$means)
plot(dvals, result$means, type = "b")