library(ggplot2)
library(gridExtra)
library(StableMCD)
library(rrcov)

data(fruit)
x = data.matrix(fruit[,2:257])

ptm <- proc.time()
result = bootstrap_mcd(x,seq(0.5,0.975,by=0.025),B=50,classifier="MD")
time <- proc.time() - ptm

data1 = data.frame(alphas = seq(0.5,0.975,by=0.025),insta=as.vector(result$means),sd = as.vector(result$sds))
plot1 = ggplot(data1, aes(x=alphas, y=insta)) + 
  geom_line(colour="navyblue") +
  geom_point(colour="navyblue")+
  geom_errorbar(aes(ymin=insta-sd,ymax=insta+sd),width=0.01)+
  labs(y = "Instability", x = "h/n", title="Instability on Fruit Data")+theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),text = element_text(size=12))

subset = mcd(x,0.825)
mds = log(mahalanobis(x, subset$muhat, subset$Sigmahat))
mds = sort(mds)
data2 = data.frame(index = 549:dim(x)[1],mds=mds[549:dim(x)[1]])
plot2 = ggplot(data2, aes(x=index, y=mds)) + 
  geom_bar(stat = "identity",color='lightblue')+
  labs(title = "Mahalanobis Distances", x = "Observation Index", y = "MD (log scale)") +theme_bw()+
  geom_vline(data=data2, aes(xintercept=result$best_alpha*dim(x)[1]),color='red')+
  theme(plot.title = element_text(hjust = 0.5),text = element_text(size=12))

grid.arrange(plot1, plot2, nrow=1, ncol=2)
