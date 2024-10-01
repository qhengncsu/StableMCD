library(DetMCD)
library(mclust)
library(StableMCD)
library(ddalpha)
library(ggplot2)
library(gridExtra)

data(banknote)
x = as.matrix(banknote[101:200,2:7])

ptm <- proc.time()
bootstrap_result = bootstrap_mcd(x,seq(0.5,0.99,0.01),B=100,classifier='depth')
time <- proc.time() - ptm

data1 = data.frame(h = seq(0.5,0.99,0.01)*dim(x)[1],insta=as.vector(bootstrap_result$means),sd = as.vector(bootstrap_result$sds))
plot1 = ggplot(data1, aes(x=h, y=insta)) + 
  geom_line(colour="navyblue") +
  geom_point(colour="navyblue")+
  geom_errorbar(aes(ymin=insta-sd,ymax=insta+sd),width=0.01)+
  labs(y = "Instability", x = "h", title="Instability on Banknote Data")+theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),text = element_text(size=12))


subset = mcd(x,0.75)
mds = mahalanobis(x, subset$muhat, subset$Sigmahat)
mds = sort(mds)
data2 = data.frame(index = 51:dim(x)[1],mds=mds[51:100])
plot2 = ggplot(data2, aes(x=index, y=mds)) + 
  geom_bar(stat = "identity",color='lightblue',fill='lightblue')+
  labs(title = "Mahalanobis Distances", x = "Observation Index", y = "MD") +theme_bw()+
  geom_vline(data=data2, aes(xintercept=bootstrap_result$best_alpha*dim(x)[1]),color='red')+
  theme(plot.title = element_text(hjust = 0.5),text = element_text(size=12))

grid.arrange(plot1, plot2, nrow=1, ncol=2)
