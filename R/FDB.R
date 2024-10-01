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
    pro <- proj_depth(x, x, 1, multiplier=floor(k/dim(x)[2]))
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
              raw.cov = sigma_raw, raw.md = MD, trunc =trunc))
}