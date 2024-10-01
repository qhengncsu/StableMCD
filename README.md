This repo implements the methods described in the following two publications:

1. A Stability Framework for Parameter Selection in the Minimum Covariance Determinant Problem

Qiang Heng, Hui Shen, Kenneth Lange (2024+)

2. Bootstrap Estimation of the Proportion of Outliers In Robust Regression

Qiang Heng, Kenneth Lange (2024+).

## Installation

```R
library(devtools)
install_github("qhengncsu/StableMCD")
```
## Quick Start for "A Stability Framework for Parameter Selection in the Minimum Covariance Determinant Problem"

### Simulated Data
```R
library(StableMCD)

# Regular Outliers
n = 1000
x1 = rnorm(n)
x1[1:(n * 0.2)] = rnorm(n * 0.2, 10)
x2 = rnorm(n)
x = cbind(x1,x2)

result = bootstrap_mcd(x,seq(0.5,0.975,by=0.025),B=200,classifier="MD")
plot(seq(0.5,0.975,by=0.025), result$means, type = "b")

# Masking outliers
n = 1000
x1 = rnorm(n)
x1[1:(n * 0.15)] = rnorm(n * 0.15, 7)
x1[(n * 0.15 + 1):(n * 0.15 + n * 0.05)] = rnorm(n * 0.05, 10^3)
x2 = rnorm(n)
x = cbind(x1,x2)

result = bootstrap_mcd(x,seq(0.5,0.975,by=0.025),B=200,classifier="MD")
plot(seq(0.5,0.975,by=0.025), result$means, type = "b")
```

### Real Data
Set working directory to "[pathofrepo]/vigenettes". The packages imported at the
beginning of each file (other than StableMCD, which is this package) are all 
readily available on CRAN.
1. Run vigenettes/star_mcd.R to produce Figure 3 for the StarsCYG data.
2. Run vigenettes/banknote.R to produce Figure 4 for the Bank Note data.
3. Run vigenettes/fruit.R to produce Figure 5 for the Fruit data (warning, it will
take about 5 minutes).
4. Run vigenettes/glass_robpca.R to produce Figure 6 for the Glass Data.
5. Run vigenettes/breast.R to produce Figure 7 for the Breast Cancer Data (warning, 
it will take about 20 minutes).


## Quick Start for "Bootstrap Estimation of the Proportion of Outliers In Robust Regression"

### Simulated Data
```R
library(StableMCD)
# generate data
X = matrix(rnorm(400*5),nrow=400)
beta = rep(1,5)+runif(5)
y = X%*%beta + rnorm(100)

#add outlier corruption to y
y[1:100,1] =  y[1:100,1] + 10

# Bootstrap!
result = bootstrap_lts(X,y,seq(0.5,0.975,by=0.025))
plot(seq(0.5,0.975,by=0.025), result$means, type = "b")
```

### Real Data
Set working directory to "[pathofrepo]/vigenettes".
1. Run vigenettes/star_lts.R to produce the instability path for the StarsCYG data.
2. Run vigenettes/glass_lts.R to produce the instability path for the Glass data.
