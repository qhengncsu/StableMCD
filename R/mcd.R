get_instability <- function(is_outliers1,is_outliers2, h){
  n <- length(is_outliers1)
  instability <- sum(abs(is_outliers1-is_outliers2))
  p <- instability/n
  c <- (choose(h,2) + choose(n-h,2))/choose(n,2)
  return(p*(1-p)/(c*(1-c))-1)
}

concentration <- function(x,index,h,max_iter=100,verbose=T){
  i <- 0
  old_sum_distances <- Inf
  while(i<=max_iter-1){
    i <- i + 1
    subset <- x[index,]
    muhat <- apply(subset,2,mean)
    Sigmahat <- cov(subset)
    MD <- mahalanobis(x, muhat, Sigmahat)
    index_new <- order(MD)[1:h]
    sum_distances <- sum(MD[index_new])
    len_indexdiff <- length(setdiff(index_new,index))
    if(verbose){
      cat(sprintf("Iteration %d: %d indices updated.\n", i, len_indexdiff))
    }
    if(len_indexdiff ==0 | abs(sum_distances-old_sum_distances) < 1e-8){
      break
    }else{
      old_sum_distances <- sum_distances
    }
    index <- index_new
  }
  return(list(index=index,muhat=muhat,Sigmahat=Sigmahat))
}

mcd <- function(x,alpha,verbose=T){
  n <- nrow(x)
  h = floor(alpha*n)
  depths = proj_depth(x,x,3,multiplier=100)
  depth_order = order(depths, decreasing = TRUE)
  index = depth_order[1:h]
  return(concentration(x,index,h,verbose=verbose))
}

bootstrap_mcd <- function(x, alphas, B=50, classifier='depth'){
  n <- nrow(x)
  instabilities = list()
  for(i in 1:length(alphas)){
    instabilities[[i]] = rep(0,B)
  }
  depths = proj_depth(x,x,3,multiplier=100)
  for(b in 1:B){
    index1_bootstrap = sample(1:n,n,replace=TRUE)
    index2_bootstrap = sample(1:n,n,replace=TRUE)
    x1 = x[index1_bootstrap,]
    x2 = x[index2_bootstrap,]
    depths1 = depths[index1_bootstrap]
    depths2 = depths[index2_bootstrap]
    #depths1 = proj_depth(x1,x1,1,multiplier=10)
    #depths2 = proj_depth(x2,x2,1,multiplier=10)
    depth_order1 = order(depths1, decreasing = TRUE)
    depth_order2 = order(depths2, decreasing = TRUE)
    for(i in 1:length(alphas)){
      h = floor(alphas[i]*n)
      index1 = depth_order1[1:h]
      index2 = depth_order2[1:h]
      result1 = concentration(x1,index1,h,verbose=F)
      result2 = concentration(x2,index2,h,verbose=F)
      index1 = result1$index
      index2 = result2$index
      if(classifier=="depth"){
        depths1 = proj_depth(x,x1[index1,],1,multiplier=100)
        depths2 = proj_depth(x,x2[index2,],1,multiplier=100)
        order1 = order(depths1,decreasing = TRUE)
        order2 = order(depths2,decreasing = TRUE)
      }else if(classifier=="MD"){
        MD1 <- mahalanobis(x, result1$muhat, result1$Sigmahat)
        MD2 <- mahalanobis(x, result2$muhat, result2$Sigmahat)
        order1 = order(MD1,decreasing = FALSE)
        order2 = order(MD2,decreasing = FALSE)
      }else{
        stop("Invalid Classifier!")
      }
      is_outliers1 = rep(1,n)
      is_outliers2 = rep(1,n)
      is_outliers1[order1[1:h]] = 0
      is_outliers2[order2[1:h]] = 0
      instabilities[[i]][b] = get_instability(is_outliers1,is_outliers2,h)
      #instabilities[[i]][b,j] = -cor(depths1,depths2,method='kendall')
    }
    if(b%%10==0){
      cat(sprintf("Bootstrap pair %d completed!\n", b))
    }
  }
  means = rep(0,length(alphas))
  sds = rep(0,length(alphas))
  for(i in 1:length(alphas)){
    #quartiles <- quantile(instabilities[[i]], probs = c(0.25, 0.5, 0.75))
    order = order(instabilities[[i]])
    #trimeans[i] = (quartiles[1]+2*quartiles[2]+quartiles[3])/4
    means[i] = mean(instabilities[[i]])
    h = floor(alphas[i]*n)
    #retain = 1 - (1 - pbinom(h,n,h/n))^2
    #means[i] = mean(instabilities[[i]][order[1:(0.95*B)]])
    sds[i] = sd(instabilities[[i]])
  }
  best_index = which(means == min(means))
  best_alpha = alphas[best_index]
  return(list(best_alpha=best_alpha,means=means,sds=sds,instabilities=instabilities))
}
