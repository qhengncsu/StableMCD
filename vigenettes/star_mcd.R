library(StableMCD)
library(DetMCD)
library(ggplot2)
library(gridExtra)
star = readRDS("star.rds")
x = as.matrix(star)
alphas = seq(25,46)/47+0.001

ptm <- proc.time()
result = bootstrap_mcd(x,alphas,B=100,classifier="depth")
time <- proc.time() - ptm

data = data.frame(h=25:46,insta_mean=result$means,insta_sd=result$sds)
plot1 = ggplot(data,aes(x=h))+
        geom_point(aes(y=insta_mean,color="bootstrap"))+
        geom_line(aes(y=insta_mean,color="bootstrap"))+
        geom_errorbar(aes(ymin=insta_mean-insta_sd,ymax=insta_mean+insta_sd,color="bootstrap"),width=0.01)+
        scale_color_manual("Method",breaks=c("bootstrap"),
                     values=c("navyblue"))+
        labs(y = "Instability", x = "h", title="Instability on Star Data")+theme_bw()+
        theme(plot.title = element_text(hjust = 0.5), legend.position = "none", text = element_text(size=12))

result = mcd(x,0.8946170)
is_outlier = rep("outlier",dim(x)[1])
is_outlier[result$index] = "inlier"
data2 = data.frame(X1 = x[,1],X2=x[,2])
data2["Class"] = is_outlier
plot2 = ggplot()+ geom_point(data=data2, aes(x=X1,y=X2,color=Class),size=3) + 
  labs(y = "Log Light", x = "Log Temperature", title="Inlier/Outlier Map for Star Data") + theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none", text = element_text(size=12)) +
  scale_color_manual("inlier/outlier",breaks=c("inlier","outlier"),
                     values=c("blue","red"))

grid.arrange(plot1, plot2, nrow=1, ncol=2)
