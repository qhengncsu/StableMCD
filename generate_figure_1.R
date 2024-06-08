source('MCD_spectral.R')
library(ggplot2)
outlier_ratio = 0.25
p = 500
l = 20

generate_corr <- function(CN,p){
  Lambda = rep(0,p)
  Lambda[1] = 1
  Lambda[p] = CN
  Lambda[2:(p-1)] = runif(p-2,1,CN)
  Lambda[2:(p-1)] = sort(Lambda[2:(p-1)])
  Y = matrix(rnorm(p*p),p)
  svd_result = svd(Y)
  U = svd_result$v
  Sigma0 = U%*%diag(Lambda)%*%t(U)
  sd_Sigma0 = sqrt(diag(Sigma0))
  R0 = diag(sd_Sigma0^{-1})%*%Sigma0%*%diag(sd_Sigma0^{-1})
  while(cond(R0)<CN-0.1 | cond(R0)>CN+0.1){
    eig_result = eigen(R0)
    Lambda = eig_result$values
    U = eig_result$vectors
    Lambda[1] = CN*Lambda[p]
    R0 = U%*%diag(Lambda)%*%t(U)
    sd_R0 = sqrt(diag(R0))
    R0 = diag(sd_R0^{-1})%*%R0%*%diag(sd_R0^{-1})
  }
  eig_result = eigen(R0)
  Lambda = eig_result$values
  U = eig_result$vectors
  return(list(R0=R0,U=U,Lambda=Lambda))
}

control_CN <- function(x, index, CN){
  Sigmahat = cov(x[index, ])
  d = eigen(Sigmahat)$values
  Lambda_max = max(d)
  Lambda_min = min(d)
  if(Lambda_max <= CN*Lambda_min){
    return(Sigmahat)
  }else{
    rho = (Lambda_max-CN*Lambda_min)/(CN-1+Lambda_max-CN*Lambda_min)
    print(rho)
    Sigmahat = rho*diag(dim(x)[2])+(1-rho)*Sigmahat
    return(Sigmahat)
  }
}

CN = 50
n = 300
corr_result = generate_corr(CN,p)
x <- mvrnorm(n,rep(0,p),corr_result$R0)
index_opt<-sample(1:n,outlier_ratio*n)
Sigma_inv = solve(corr_result$R0)
for(j in index_opt){
  if(l>1){
    dir_index = sample(c(seq(p-l+1,p)),1)
  }else{
    dir_index = p
  }
  direction = corr_result$U[,dir_index]
  x[j,] <- x[j,]+50*direction
}

bootstrap_result = bootstrap_insta(x,alpha_list=seq(0.5,0.95,by=0.05),q_list=c(2,10,50))

data = data.frame(h=seq(0.5,0.95,by=0.05),mean_q2 = bootstrap_result$means[1,],sd_q2 = bootstrap_result$sds[1,],
                   mean_q10 = bootstrap_result$means[2,],sd_q10 = bootstrap_result$sds[2,],
                   mean_q50 = bootstrap_result$means[3,],sd_q50 = bootstrap_result$sds[3,])

plot  = ggplot(data,aes(x=h))+
  geom_point(aes(y=mean_q2,color="q=2"))+
  geom_point(aes(y=mean_q10,color="q=10"))+
  geom_point(aes(y=mean_q50,color="q=50"))+
  geom_line(aes(y=mean_q2,color="q=2"))+
  geom_line(aes(y=mean_q10,color="q=10"))+
  geom_line(aes(y=mean_q50,color="q=50"))+
  scale_color_manual("Number of PCs",breaks=c("q=2","q=10","q=50"),
                     values=c("navyblue","darkgreen","darkred"))+
  labs(y = "Instability", x = "h/n", title="Instability")+theme_bw()+
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none", text = element_text(size=8))