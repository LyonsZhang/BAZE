
setwd('C:\Users\lzhang27\Desktop\BayesianCompositionSelection\RealDataRelated')
library("ggcorrplot")
library("R.matlab")
library("Cairo")
Qorg<-readMat("Qorg.mat")
Qinv<-readMat("Qinv.mat")



quantile(Qorg$Q.full)
f<-ecdf(Qorg$Q.full)
f(0.5)
# te<-newftree$tip.label
# te[is.na(te)] <- ""
# colnames(Qorg$Q.full)<-te

###############plot the correlation and precision matrix

setEPS()
postscript("Euclideancorrelation.eps")
hist(Qorg$Q.full,freq=T,col="blue",xlab="Euclidean correlation",cex.lab=1.5,cex.axis=1.5,font.lab=2,main=NULL)
abline(v = 0.5, col="red", lwd=3, lty=2)
text(0.52, 4e+05,"88% quantile",cex=1.5)
dev.off()

tt<-Qinv$Q.full
tt[tt>1]<-1
tt[tt< -1]<- -1
te<-newftree$tip.label
te[is.na(te)] <- ""
#colnames(tt)<-te
te[793]<-""
te[205]<-""
te[250]<-"Coprococcus"
te[326]<-""
te[360]<-"Faecalibacterium"
quantile(tt)
f<-ecdf(tt)
f(0)
setEPS()
postscript("Euclideanprecision.eps")
hist(tt,freq=T,col="blue",xlab="Euclidean precision",cex.lab=1.5,cex.axis=1.5,font.lab=2,main=NULL)
abline(v = 0, col="red", lwd=3, lty=2)
text(0, 4e+05,"56.7% quantile",cex=1.5)
dev.off()




#ggcorrplot(Qorg$Q.full,outline.col = "white",insig = "blank",colors = c("red", "white", "#E46726"))
p1<-pheatmap::pheatmap(Qorg$Q.full, cluster_rows=F,cluster_cols=F,border_color = NA,
                   color = colorRampPalette(c("white", "mintcream",  "gold", "orange", "red"))(100),fontsize = 18,labels_col=te,angle_col = "270")
ggsave(filename="Correlationheatmap.eps",plot=p1,width = 10, height = 10,dpi=150,device = cairo_ps)
ggsave(filename="Correlationheatmap.pdf",plot=p1,dpi=150,width = 10, height = 10)
ggsave(filename="Correlationheatmap.jpg",plot=p1,width = 10, height = 10)

p2<-pheatmap::pheatmap(tt, cluster_rows=F,cluster_cols=F,border_color = NA,
                   color = colorRampPalette(c("darkslateblue", "blue",  "white", "orange", "red"))(100),fontsize = 18,labels_col=te,angle_col = "270")
ggsave(filename="Precisionheatmap.eps",plot=p2,width = 10, height = 10,dpi=150,device = cairo_ps)
ggsave(filename="Precisionheatmap.pdf",plot=p2,width = 10, height = 10)
ggsave(filename="Precisionheatmap.jpg",plot=p2,width = 10, height = 10)
