source("MCD_spectral.R")
library(ggplot2)
library(gridExtra)
data(fruit)
x = data.matrix(fruit[,2:257])
x_center = apply(x,2,mean)
x_centered = sweep(x,2,x_center)
svd_result = svd(x_centered)
PCs <- x_centered%*%svd_result$v
data1 = data.frame(X1 = PCs[,1],X2=PCs[,2],cultivar=lapply(fruit["cultivar"],as.character))

ptm <- proc.time()
result = bootstrap_insta(x,alpha_list=seq(0.5,0.95,by=0.05),q_list=2,B=50)
time <- proc.time() - ptm

FDB <- function (x, alpha = 0.75, depth = "pro", k = 1000) 
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
  if (depth == "l2") {
    l2 <- l2_fast(x)
    index11 <- order(l2, decreasing = TRUE)[1:(alpha * n)]
  }
  subset <- x[index11, ]
  hat_mu <- colMeans(subset)
  hat_sigma <- cov(subset)
  MD_c0 <- mahalanobis(x, hat_mu, hat_sigma)
  c_s <- median(MD_c0)/qchisq(0.5, p)
  sigma_raw <- c_s * hat_sigma
  MD <- mahalanobis(x, hat_mu, sigma_raw)
  trunc <- which(MD <= qchisq(0.975, p))
  x_trunc <- x[trunc, ]
  center <- colMeans(x_trunc)
  cov <- cov(x_trunc)
  return(list(center = center, cov = cov, best = trunc, raw.center = hat_mu, 
              raw.cov = sigma_raw, raw.md = MD))
}


pro <- FDB(x,alpha=0.5,k=10*dim(x)[2],depth="pro")
spectral <- MCD_spectral(x,q=2,alpha=0.85)
is_outlier_pro = rep("outlier",dim(x)[1])
is_outlier_pro[pro$best] = "inlier"
is_outlier_spectral = rep("outlier",dim(x)[1])
is_outlier_spectral[spectral$best] = "inlier"

data2 = data.frame(h = seq(0.5,0.95,by=0.05),insta=as.vector(result$means),sd = as.vector(result$sds))
plot1 = ggplot()+ geom_point(data=data1, aes(x=X1,y=X2,color=cultivar),size=0.5) + 
        labs(y = "PC2", x = "PC1", title="Cultivar") + theme_bw() + 
        theme(plot.title = element_text(hjust = 0.5), legend.position = "none", text = element_text(size=8)) +
        scale_color_manual(breaks=c("D","HA","M"),
                           values=c("blue","red","green"))
plot2 = ggplot(data2, aes(x=h, y=insta)) + 
        geom_line(colour="navyblue") +
        geom_point(colour="navyblue")+
        labs(y = "Instability", x = "h/n", title="Instability")+theme_bw()+
        theme(plot.title = element_text(hjust = 0.5),text = element_text(size=8))
data1["Class"] = is_outlier_pro 
# plot3 = ggplot()+ geom_point(data=data1, aes(x=X1,y=X2,color=Class)) + 
#         labs(y = "PC2", x = "PC1", title="FDB") + theme_bw() + 
#         theme(plot.title = element_text(hjust = 0.5), legend.position = c(0.9,0.88), text = element_text(size=12)) +
#         scale_fill_discrete(name = "Class") +
#         scale_color_manual(breaks=c("inlier","outlier"),
#                            values=c("blue","red"))

data1["Class"] = is_outlier_spectral
plot3 = ggplot()+ geom_point(data=data1, aes(x=X1,y=X2,color=Class),size=0.5) + 
        labs(y = "PC2", x = "PC1", title="SpectralMCD") + theme_bw() + 
        theme(plot.title = element_text(hjust = 0.5), legend.position = "none", text = element_text(size=8)) +
        scale_color_manual("inlier/outlier",breaks=c("inlier","outlier"),
                          values=c("blue","red"))

grid.arrange(plot1, plot2,plot3, nrow=1, ncol=3)
