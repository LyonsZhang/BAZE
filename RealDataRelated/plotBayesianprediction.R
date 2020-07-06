
setwd('C:\Users\lzhang27\Desktop\BayesianCompositionSelection\RealDataRelated')
library("R.matlab")
mydata<-readMat("realestimation.mat")

###########do the predition plot of Bayesian method

MSE<-21.2288
#MSE=mean(mydata$MSE[15000:20000])
df<-data.frame(cbind(mydata$Ytestobs,mydata$Ytestest))
names(df)<-c("Yobs","Yest")
r=cor(df$Yobs,df$Yest)
library(ggplot2)
library(Cairo)
ggplot(df,aes(x=Yobs,y=Yest))+
#  xlim(12,45)+
#  ylim(12,45)+
  scale_y_continuous(breaks=seq(15,45,5))+
  scale_x_continuous(breaks=seq(15,45,5))+
  geom_point(fill="blue",size=5,pch=21,stroke=1,color = "black",alpha=0.2)+
  geom_segment(aes(x = 12, xend = 45, y = 12, yend = 45),
               colour = "red",lty = "dashed",size=2)+
  annotate("text", x=40, y = 20, label = paste("MSE = ", round(MSE,2),sep=""),colour = "black", size = 8,fontface=2)+
  annotate("text", x=40, y = 22, label = paste("r = ", round(r,2),sep=""),colour = "black", size = 8,fontface=2)+
  labs(x = "Observed BMI", 
       y = "Fitted BMI")+
  theme_bw()+
  # theme(plot.title = element_text(size = 30, face = "bold"),
  #       legend.title=element_text(size=20),
  #       legend.text=element_text(size=16))+
  theme(axis.line = element_line(colour = "black"),
        #panel.grid.major = element_line(size = 0.2),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.line.x = element_line(color="black", size = 1.5),
        #axis.line.x.top = element_line(color="black", size = 1.5),
        axis.line.y = element_line(color="black", size = 1.5),
        axis.text=element_text(size=16,face="bold"),
        axis.title=element_text(size=20,face="bold"))
ggsave(filename=paste("bayesianrealfittedplot",".eps", sep=""),width = 10, height = 7,device = cairo_ps)

###########do the predition plot of LASSO method

setwd("\\LASSO\\")
df<-data.frame(cbind(ytest,fitted))
names(df)<-c("Yobs","Yest")
MSE<-1/length(df$Yobs)*t(df$Yobs-df$Yest)%*%(df$Yobs-df$Yest)
library(ggplot2)
library(Cairo)
r=cor(df$Yobs,df$Yest)
ggplot(df,aes(x=Yobs,y=Yest))+
  #  xlim(12,45)+
  #  ylim(12,45)+
  scale_y_continuous(breaks=seq(15,45,5))+
  scale_x_continuous(breaks=seq(15,45,5))+
  geom_point(fill="blue",size=5,pch=21,stroke=1,color = "black",alpha=0.2)+
  geom_segment(aes(x = 12, xend = 45, y = 12, yend = 45),
               colour = "red",lty = "dashed",size=2)+
  annotate("text", x=40, y = 20, label = paste("MSE = ", round(MSE,2),sep=""),colour = "black", size = 8,fontface=2)+
  annotate("text", x=40, y = 22, label = paste("r = ", round(r,2),sep=""),colour = "black", size = 8,fontface=2)+
  labs(x = "Observed BMI", 
       y = "Fitted BMI")+
  theme_bw()+
  # theme(plot.title = element_text(size = 30, face = "bold"),
  #       legend.title=element_text(size=20),
  #       legend.text=element_text(size=16))+
  theme(axis.line = element_line(colour = "black"),
        #panel.grid.major = element_line(size = 0.2),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.line.x = element_line(color="black", size = 1.5),
        #axis.line.x.top = element_line(color="black", size = 1.5),
        axis.line.y = element_line(color="black", size = 1.5),
        axis.text=element_text(size=16,face="bold"),
        axis.title=element_text(size=20,face="bold"))
ggsave(filename=paste("weilinLASSOrealfittedplot",".eps", sep=""),width = 10, height = 7,device = cairo_ps)


