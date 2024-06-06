source("MCD_spectral.R")
library(ggplot2)
library(gridExtra)

findPsix = function(y, h,
                    x) {
  # computes psi(x) for the estimates based on 
  # the bootstrapped sample y
  
  # adapted from rrcov:::unimcd
  
  # y is bootstrapped sample
  # h is size of h-subset
  # x is original data sample
  
  quan = h
  out = list()
  ncas = length(y)
  len = ncas - quan + 1
  if (len == 1) {
    out$tmcd = mean(y)
    out$smcd = sqrt(var(y))
  }
  else {
    ay = c()
    I = order(y)
    y = y[I]
    ay[1] = sum(y[1:quan])
    for (samp in 2:len) {
      ay[samp] = ay[samp - 1] - y[samp - 1] + y[samp + 
                                                  quan - 1]
    }
    ay2 = ay^2/quan
    sq = c()
    sq[1] = sum(y[1:quan]^2) - ay2[1]
    for (samp in 2:len) {
      sq[samp] = sq[samp - 1] - y[samp - 1]^2 + y[samp + 
                                                    quan - 1]^2 - ay2[samp] + ay2[samp - 1]
    }
    sqmin = min(sq)
    Isq = order(sq)
    ndup = sum(sq == sqmin)
    ii = Isq[1:ndup]
    # 
    slutn = c()
    slutn[1:ndup] = ay[ii]
    initmean = slutn[floor((ndup + 1)/2)]/quan
    initmean = median(y[ii:(ii+h-1)])
    psix = rep(1, length(x))
    psix[order(abs(x-initmean))[1:h]] = 0
    
  }
  return(psix)
}

n = 10^3 # sample size
hvals = seq(0.5,0.95,by=0.05) * n # h-values
B = 100

x = rnorm(n)
x[1:(n*0.1)] = rnorm(n*0.1, 10);
x[(n*0.1+1):(n*0.2)] = rnorm(n*0.1, -10)

result = matrix(0, length(hvals), B)
for (j in 1:length(hvals)) {
  h = hvals[j]
  cprime = (choose(h, 2) + choose(n - h, 2)) / choose(n, 2)
  for (i in 1:B) {
    xdot  = sample(x, n, replace = TRUE)
    xddot = sample(x, n, replace = TRUE)
    dx = mean(abs(findPsix(y = xdot,h = h, x = x) - 
                    findPsix(y = xddot, h = h, x = x)))
    result[j, i] = dx*(1-dx) / (cprime *(1-cprime)) - 1
  }
}

data1 = data.frame(h = hvals, insta = rowMeans(result), sd = apply(result,1,std)) 
plot1 = ggplot(data1, aes(x=h, y=insta)) +
  geom_point() +
  geom_line() +
  labs(y = "Instability", x = "h")+theme_bw()+
  theme(text = element_text(size=8))

x = rnorm(n)
x[1:(n * 0.15)] = rnorm(n * 0.15, 5)
x[(n * 0.15 + 1):(n * 0.15 + n * 0.05)] = rnorm(n * 0.05, 10^3)


result = matrix(0, length(hvals), B)
for (j in 1:length(hvals)) {
  h = hvals[j]
  cprime = (choose(h, 2) + choose(n - h, 2)) / choose(n, 2)
  for (i in 1:B) {
    xdot  = sample(x, n, replace = TRUE)
    xddot = sample(x, n, replace = TRUE)
    dx = mean(abs(findPsix(y = xdot,h = h, x = x) - 
                    findPsix(y = xddot, h = h, x = x)))
    result[j, i] = dx*(1-dx) / (cprime *(1-cprime)) - 1
  }
}

data2 = data.frame(h = hvals, insta = rowMeans(result), sd = apply(result,1,std)) 
plot2 = ggplot(data2, aes(x=h, y=insta)) +
  geom_point() +
  geom_line() +
  labs(y = "Instability", x = "h")+theme_bw()+
  theme(text = element_text(size=8))

x1 = rnorm(n)
x1[1:(n*0.1)] = rnorm(n*0.1, 10);
x1[(n*0.1+1):(n*0.2)] = rnorm(n*0.1, -10)
x2 = rnorm(n)
x2[1:(n*0.1)] = rnorm(n*0.1, 10);
x2[(n*0.1+1):(n*0.2)] = rnorm(n*0.1, -10)
x = cbind(x1,x2)

bootstrap_result = bootstrap_insta(x,alpha_list = seq(0.5,0.95,by=0.05),q_list=2)
data3 = data.frame(h = hvals, insta = as.vector(bootstrap_result$means), sd = as.vector(bootstrap_result$sds))
plot3 = ggplot(data3, aes(x=h, y=insta)) +
  geom_point() +
  geom_line() +
  labs(y = "Instability", x = "h")+ theme_bw()+
  theme(text = element_text(size=8))

x1 = rnorm(n)
x1[1:(n * 0.15)] = rnorm(n * 0.15, 5)
x1[(n * 0.15 + 1):(n * 0.15 + n * 0.05)] = rnorm(n * 0.05, 10^3)
x2 = rnorm(n)
x2[1:(n * 0.15)] = rnorm(n * 0.15, 5)
x2[(n * 0.15 + 1):(n * 0.15 + n * 0.05)] = rnorm(n * 0.05, 10^3)
x = cbind(x1,x2)

bootstrap_result = bootstrap_insta(x,alpha_list = seq(0.5,0.95,by=0.05),q_list=2)
data4 = data.frame(h = hvals, insta = as.vector(bootstrap_result$means), sd = as.vector(bootstrap_result$sds))
plot4 = ggplot(data4, aes(x=h, y=insta)) +
  geom_point() +
  geom_line() +
  labs(y = "Instability", x = "h")+ theme_bw()+
  theme(text = element_text(size=8))


grid.arrange(plot1, plot2,plot3, plot4, nrow=2, ncol=2)
