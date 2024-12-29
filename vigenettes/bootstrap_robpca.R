library(rrcov)
bootstrap_robpca <- function(x, alphas, q, B=50, classifier='MD'){
  n <- nrow(x)
  instabilities = list()
  wds = list()
  for(i in 1:length(alphas)){
    instabilities[[i]] = rep(0,B)
    wds[[i]] = rep(0,B)
  }
  for(b in 1:B){
    index1_bootstrap = sample(1:n,n,replace=TRUE)
    index2_bootstrap = sample(1:n,n,replace=TRUE)
    x1 = x[index1_bootstrap,]
    x2 = x[index2_bootstrap,]
    res1 = PcaHubert(x1,k=q,alpha=0.5)
    res2 = PcaHubert(x2,k=q,alpha=0.5)
    x1_center = res1$center
    x2_center = res2$center
    v1 = res1$loadings
    v2 = res2$loadings
    PCs1_b = sweep(x1,2,x1_center)%*%v1
    PCs2_b = sweep(x2,2,x2_center)%*%v2
    PCs1 = sweep(x,2,x1_center)%*%v1
    PCs2 = sweep(x,2,x2_center)%*%v2
    depths1 = proj_depth(PCs1_b, PCs1_b, 1, multiplier=100)
    depths2 = proj_depth(PCs2_b, PCs2_b, 1, multiplier=100)
    depth_order1 = order(depths1, decreasing = TRUE)
    depth_order2 = order(depths2, decreasing = TRUE)
    for(i in 1:length(alphas)){
      alpha = alphas[i]
      h = floor(alpha*n)
      index1 = depth_order1[1:h]
      index2 = depth_order2[1:h]
      result1 = concentration(PCs1_b,index1,h,verbose=F)
      result2 = concentration(PCs2_b,index2,h,verbose=F)
      index1 = result1$index
      index2 = result2$index
      if(classifier=='depth'){
        depths3 = proj_depth(PCs1,PCs1_b[index1,],1,multiplier=100)
        depths4 = proj_depth(PCs2,PCs2_b[index2,],1,multiplier=100)
        order1 = order(depths3,decreasing=TRUE)
        order2 = order(depths4,decreasing=TRUE)
      }else{
        MD1 <- mahalanobis(PCs1, result1$muhat, result1$Sigmahat)
        MD2 <- mahalanobis(PCs2, result2$muhat, result2$Sigmahat)
        order1 = order(MD1)
        order2 = order(MD2)
      }
      is_outliers1 = rep(1,n)
      is_outliers2 = rep(1,n)
      is_outliers1[order1[1:h]] = 0
      is_outliers2[order2[1:h]] = 0
      instabilities[[i]][b] = get_instability(is_outliers1,is_outliers2,h)
      wds[[i]][b] = log(wasserstein_distance(result1$muhat,result1$Sigmahat,result2$muhat,result2$Sigmahat))
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
    insta_means[i] = mean(instabilities[[i]])
    insta_sds[i] = sd(instabilities[[i]])
    wd_means[i] = mean(wds[[i]])
    wd_sds[i] = sd(wds[[i]])
    #h = floor(alphas[i]*n)
  }
  scaled_wd_means = (wd_means - min(wd_means))/(max(wd_means)-min(wd_means))
  final_score = insta_means + 0.5*median(insta_means)*scaled_wd_means
  best_index = which(final_score == min(final_score))
  best_alpha = alphas[best_index]
  return(list(best_alpha=best_alpha,final_score=final_score,insta_means=insta_means,insta_sds=insta_sds,
              wd_means=wd_means,wd_sds=wd_sds,instabilities=instabilities,wds = wds))
}