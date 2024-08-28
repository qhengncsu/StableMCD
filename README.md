This repo provides code for two publications

1. High-dimensional Outlier Detection via Stability

Qiang Heng, Hui Shen, Kenneth Lange (2024+), arXiv: [2401.14359](https://arxiv.org/abs/2401.14359v3)

2. Bootstrap Estimation of the Proportion of Outliers In Robust Regression

Qiang Heng, Kenneth Lange (2024+).

## Installation

```R
library(devtools)
install_github("qhengncsu/StableMCD")
```

## Quick Start for "Bootstrap Estimation of the Proportion of Outliers In Robust Regression"

```R
library(StableMCD)
# generate data
X = matrix(rnorm(400*5),nrow=400)
beta = rep(1,5)
y = X%*%beta + rnorm(100)

#add outlier corruption to y
y[1:40] = y[1:40] + 10

#add outlier corruption to X
X[41:80,1] =  X[41:80,1] + 10

# Bootstrap!
result = bootstra_lts(X,y,seq(0.5,0.975,by=0.025))
plot(seq(0.5,0.975,by=0.025), result$insta_means)
```

