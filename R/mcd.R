get_instability <- function(is_outliers1,is_outliers2, h){
  n <- length(is_outliers1)
  instability <- sum(abs(is_outliers1-is_outliers2))
  p <- instability/n
  c <- (choose(h,2) + choose(n-h,2))/choose(n,2)
  #return(p*(1-p)/(c*(1-c))-1)
  return(p/(1-c))
}

concentration <- function(x,index,h,max_iter=100,verbose=T){
  i <- 0
  old_sum_distances <- Inf
  while(i<=max_iter-1){
    i <- i + 1
    subset <- x[index,]
    muhat <- apply(subset,2,mean)
    Sigmahat <- cov(subset)*nrow(x)/(nrow(x)-1)
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

mcd <- function(x,alpha,verbose=T, reweighting=F){
  n <- nrow(x)
  h = floor(alpha*n)
  depths = proj_depth(x,x,3,multiplier=100)
  depth_order = order(depths, decreasing = TRUE)
  index = depth_order[1:h]
  res = concentration(x,index,h,verbose=verbose)
  if(!reweighting){
    return(res)
  }else{
    MD_c0 <- mahalanobis(x, res$muhat, res$Sigmahat)
    c_s <- median(MD_c0)/qchisq(0.5, p)
    sigma_raw <- c_s * res$Sigmahat
    MD <- mahalanobis(x, res$muhat, sigma_raw)
    trunc <- which(MD >= qchisq(0.975, p))
    new_index = setdiff(res$index,trunc)
    subset <- x[new_index,]
    muhat <- apply(subset,2,mean)
    Sigmahat <- cov(subset)*nrow(x)/(nrow(x)-1)
    return(list(index=new_index,muhat=muhat,Sigmahat=Sigmahat))
  }
}

l1_depth <- function (x, data){
  if (!(is.matrix(data) && is.numeric(data) || is.data.frame(data) && 
        prod(sapply(data, is.numeric))) || ncol(data) < 2) {
    stop("Argument \"data\" should be a numeric matrix of at least 2-dimensional data")
  }
  if (is.data.frame(data)) 
    data = data.matrix(data)
  if (!is.matrix(x)) {
    if (is.vector(x)) 
      x <- matrix(x, nrow = 1)
    if (is.data.frame(x)) 
      x = data.matrix(x)
  }
  mean <- colMeans(data)
  cov <- cov(data)
  if (sum(is.na(cov)) == 0) {
    cov.eig <- eigen(cov)
    B <- cov.eig$vectors %*% diag(sqrt(cov.eig$values))
    lambda <- solve(B)
  }
  else {
    lambda = diag(ncol(data))
  }
  depths <- rep(-1, nrow(x))
  x_scaled <- x %*% t(lambda)
  data_scaled <- data %*% t(lambda) 
  for (i in 1:nrow(x)){
    #tmp1 <- t(lambda %*% (x[i, ] - t(data)))
    tmp1 <- t((x_scaled[i, ] - t(data_scaled)))
    tmp2 <- 1/sqrt(rowSums(tmp1^2))
    tmp2[is.infinite(tmp2)] <- 0
    depths[i] <- 1 - sqrt(sum((colSums(tmp2 * tmp1)/nrow(data))^2))
  }
  return(depths)
}

# Helper function to compute matrix square root
matrix_sqrt <- function(mat) {
  eig <- eigen(mat)
  eig$values[eig$values < 0] <- 0 # Handle numerical precision issues
  eig$vectors %*% diag(sqrt(eig$values)) %*% t(eig$vectors)
}

wasserstein_distance <- function(mu1, Sigma1, mu2, Sigma2) {
  # Ensure the inputs are numeric
  mu1 <- as.numeric(mu1)
  mu2 <- as.numeric(mu2)
  Sigma1 <- as.matrix(Sigma1)
  Sigma2 <- as.matrix(Sigma2)
  
  # Check dimensional consistency
  if (length(mu1) != length(mu2) || nrow(Sigma1) != ncol(Sigma1) || 
      nrow(Sigma2) != ncol(Sigma2) || nrow(Sigma1) != nrow(Sigma2)) {
    stop("Dimension mismatch between mean vectors and covariance matrices.")
  }
  
  # Calculate mean term (L2 norm squared)
  mean_diff <- sum((mu1 - mu2)^2)
  
  # Compute square root of Sigma1
  Sigma1_sqrt <- matrix_sqrt(Sigma1)
  
  # Middle term: Sigma1_sqrt %*% Sigma2 %*% Sigma1_sqrt
  middle <- Sigma1_sqrt %*% Sigma2 %*% Sigma1_sqrt
  middle_sqrt <- matrix_sqrt(middle)
  
  # Compute trace term
  trace_term <- sum(diag(Sigma1 + Sigma2 - 2 * middle_sqrt))
  
  # Wasserstein distance
  W2 <- mean_diff + trace_term
  return(sqrt(W2)) # Return the 2-Wasserstein distance
}


bootstrap_mcd <- function(x, alphas, B=50, classifier='MD', sd_ratio=3){
  n <- nrow(x)
  instabilities = list()
  wds = list()
  for(i in 1:length(alphas)){
    instabilities[[i]] = rep(0,B)
    wds[[i]] = rep(0,B)
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
      if(classifier=="proj_depth"){
        depths1 = proj_depth(x,x1[index1,],1,multiplier=100)
        depths2 = proj_depth(x,x2[index2,],1,multiplier=100)
        order1 = order(depths1,decreasing = TRUE)
        order2 = order(depths2,decreasing = TRUE)
      }else if(classifier=="MD"){
        MD1 <- mahalanobis(x, result1$muhat, result1$Sigmahat)
        MD2 <- mahalanobis(x, result2$muhat, result2$Sigmahat)
        order1 = order(MD1,decreasing = FALSE)
        order2 = order(MD2,decreasing = FALSE)
      }else if(classifier=="l1_depth"){
        depths1 = l1_depth(x,x1[index1,])
        depths2 = l1_depth(x,x2[index2,])
        order1 = order(depths1,decreasing = TRUE)
        order2 = order(depths2,decreasing = TRUE)
      }else{
        stop("Invalid Classifier!")
      }
      is_outliers1 = rep(1,n)
      is_outliers2 = rep(1,n)
      is_outliers1[order1[1:h]] = 0
      is_outliers2[order2[1:h]] = 0
      instabilities[[i]][b] = log(1+get_instability(is_outliers1,is_outliers2,h))
      wds[[i]][b] = log(1+wasserstein_distance(result1$muhat,result1$Sigmahat,result2$muhat,result2$Sigmahat))
    }
    if(b%%10==0){
      cat(sprintf("Bootstrap pair %d completed!\n", b))
    }
  }
  insta_means = rep(0,length(alphas))
  insta_sds = rep(0,length(alphas))
  wd_means = rep(0,length(alphas))
  wd_sds = rep(0,length(alphas))
  for(i in 1:length(alphas)){
    order = order(instabilities[[i]])
    insta_means[i] = mean(instabilities[[i]],na.rm=T)
    insta_sds[i] = sd(instabilities[[i]],na.rm=T)
    wd_means[i] = mean(wds[[i]],na.rm=T)
    wd_sds[i] = sd(wds[[i]],na.rm=T)
    #h = floor(alphas[i]*n)
  }
  #scaled_wd_means = (wd_means - min(wd_means))/(max(wd_means)-min(wd_means))
  #scaled_instas = (insta_means - min(insta_means))/(max(insta_means)-min(insta_means))
  sd_instas = sd(insta_means)
  sd_wd = sd(wd_means)
  beta = sd_instas/(sd_instas+sd_ratio*sd_wd)
  iim = (1-beta)*insta_means + beta*(wd_means - min(wd_means))
  #final_score = scaled_wd_means+scaled_instas
  best_index = which(iim == min(iim))
  best_alpha = alphas[best_index]
  return(list(best_alpha=best_alpha,iim=iim,insta_means=insta_means,insta_sds=insta_sds,
              wd_means=wd_means,wd_sds=wd_sds,instabilities=instabilities,wds = wds))
}
