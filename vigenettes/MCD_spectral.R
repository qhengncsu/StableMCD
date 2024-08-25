library(ddalpha)
library(MASS)
library(DetMCD)
library(FDB)
library(pracma)
library(rrcov)
library(class)

# diag_concentration <- function(PCs,index,max_iter=100,verbose=FALSE){
#   i <- 0
#   old_sum_distances <- Inf
#   while(i<=max_iter-1){
#     i <- i + 1
#     subset <- PCs[index,]
#     muhat <- apply(subset,2,mean)
#     sigma2hat <- apply(subset,2,var)
#     PCs_minus_muhat <- sweep(PCs,2,muhat)
#     distances <- rowSums(sweep(PCs_minus_muhat,2,1/sigma2hat,'*') * PCs_minus_muhat)
#     index_new <- order(distances)[1:length(index)]
#     sum_distances <- sum(distances[index_new])
#     len_indexdiff <- length(setdiff(index_new,index))
#     if(verbose){
#       cat(sprintf("Iteration %d: %d indices updated.\n", i, len_indexdiff))
#     }
#     if(len_indexdiff ==0 | abs(sum_distances-old_sum_distances) < 1e-4){
#       break
#     }else{
#       old_sum_distances <- sum_distances
#     }
#     index <- index_new
#   }
#   return(index)
# }

regular_concentration <- function(PCs,index,max_iter=100,verbose=FALSE){
  i <- 0
  old_sum_distances <- Inf
  while(i<=max_iter-1){
    i <- i + 1
    subset <- PCs[index,]
    muhat <- apply(subset,2,mean)
    Sigmahat <- cov(subset)
    MD <- mahalanobis(PCs, muhat, Sigmahat)
    index_new <- order(MD)[1:length(index)]
    sum_distances <- sum(MD[index_new])
    len_indexdiff <- length(setdiff(index_new,index))
    if(verbose){
      cat(sprintf("Iteration %d: %d indices updated.\n", i, len_indexdiff))
    }
    if(len_indexdiff ==0 | abs(sum_distances-old_sum_distances) < 1e-4){
      break
    }else{
      old_sum_distances <- sum_distances
    }
    index <- index_new
  }
  return(list(index=index,muhat=muhat,Sigmahat=Sigmahat))
}

MCD_spectral <- function(x, alpha = 0.75, q = min(dim(x)[1],dim(x)[2]), 
                         concentration = "regular", k = NULL, verbose=TRUE) 
{
  n <- nrow(x)
  p <- ncol(x)
  x_center = apply(x,2,mean)
  x_centered = sweep(x,2,x_center)
  svd_result = svd(x_centered)
  PCs <- x_centered%*%svd_result$v[,1:q]
  if(is.null(k)){
    k <- max(1000,10*q)
  }
  pro <- depth.projection(PCs, PCs, num.directions = k)
  depth_order <- order(pro, decreasing = TRUE)
  index <- depth_order[1:(alpha * n)]
  result <- regular_concentration(PCs,index,verbose=verbose)
  return(list(muhat=result$muhat,Sigmahat=result$Sigmahat,best=result$index,
              depth_order=depth_order,
              PCs=PCs,x_center=x_center,v=svd_result$v))
}

get_instability <- function(is_outliers1,is_outliers2, h){
  n <- length(is_outliers1)
  instability <- sum(abs(is_outliers1-is_outliers2))
  p <- instability/n
  c <- (choose(h,2) + choose(n-h,2))/choose(n,2)
  return(p*(1-p)/(c*(1-c))-1)
}

bootstrap_insta <- function(x, alpha_list, q_list,k=NULL, B=100){
  n <- nrow(x)
  instabilities = list()
  for(i in 1:length(q_list)){
    instabilities[[i]] = matrix(0,B,length(alpha_list))
  }
  for(b in 1:B){
    x1 = x[sample(1:n,n,replace=TRUE),]
    x2 = x[sample(1:n,n,replace=TRUE),]
    x1_center = colMeans(x1)
    x2_center = colMeans(x2)
    x1_centered = sweep(x1,2,x1_center)
    x2_centered = sweep(x2,2,x2_center)
    v1 = svd(x1_centered)$v
    v2 = svd(x2_centered)$v
    PCs1_b = x1_centered%*%v1
    PCs2_b = x2_centered%*%v2
    PCs1 = sweep(x,2,x1_center)%*%v1
    PCs2 = sweep(x,2,x2_center)%*%v2
    for(i in 1:length(q_list)){
      q = q_list[i]
      if(is.null(k)){
        k <- max(1000,10*q)
      }
      depths1 = depth.projection(PCs1_b[,1:q], PCs1_b[,1:q], num.directions = k)
      depths2 = depth.projection(PCs2_b[,1:q], PCs2_b[,1:q], num.directions = k)
      depth_order1 = order(depths1, decreasing = TRUE)
      depth_order2 = order(depths2, decreasing = TRUE)
      for(j in 1:length(alpha_list)){
        alpha = alpha_list[j]
        alphan = floor(alpha*n)
        index1 = depth_order1[1:alphan]
        index2 = depth_order2[1:alphan]
        result1 = regular_concentration(PCs1_b[,1:q],index1)
        result2 = regular_concentration(PCs2_b[,1:q],index2)
        index1 = result1$index
        index2 = result2$index
        depths1 = depth.projection(PCs1[,1:q],PCs1_b[index1,1:q],num.directions = k)
        depths2 = depth.projection(PCs2[,1:q],PCs2_b[index2,1:q],num.directions = k)
        # order1 = order(depths1,decreasing = TRUE)
        # order2 = order(depths2,decreasing = TRUE)
        #result1 =  regular_concentration(PCs1_b[,1:q],index1)
        #result2 =  regular_concentration(PCs2_b[,1:q],index2)
        #MD1 <- mahalanobis(PCs1[,1:q], result1$muhat, result1$Sigmahat)
        #MD2 <- mahalanobis(PCs2[,1:q], result2$muhat, result2$Sigmahat)
        #depths1 = depth.projection(PCs1[,1:q],PCs1_b[index1,1:q],num.directions = k)
        #depths2 = depth.projection(PCs2[,1:q],PCs2_b[index2,1:q],num.directions = k)
        order1 = order(depths1,decreasing = TRUE)
        order2 = order(depths2,decreasing = TRUE)
        is_outliers1 = rep(1,n)
        is_outliers2 = rep(1,n)
        is_outliers1[order1[1:(n*alpha)]] = 0
        is_outliers2[order2[1:(n*alpha)]] = 0
        instabilities[[i]][b,j] = get_instability(is_outliers1,is_outliers2,alphan)
        #instabilities[[i]][b,j] = -cor(depths1,depths2,method='kendall')
      }
    }
  }
  means = matrix(0,length(q_list),length(alpha_list))
  sds = matrix(0,length(q_list),length(alpha_list))
  for(i in 1:length(q_list)){
    for(j in 1:length(alpha_list)){
      means[i,j] = mean(instabilities[[i]][,j])
      sds[i,j] = sd(instabilities[[i]][,j])
    }
  }
  best_index = which(means == min(means), arr.ind = TRUE)[1,]
  best_q = q_list[best_index[1]]
  best_alpha = alpha_list[best_index[2]]
  return(list(means=means,sds=sds,best_q=best_q,best_alpha=best_alpha))
}


