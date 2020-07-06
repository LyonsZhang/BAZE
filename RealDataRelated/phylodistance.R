rm(list=ls())

library("phyloseq")
#packageVersion("phyloseq")
library(ape)
#library(adephylo)
#library(ggplot2)
#library(RRphylo)
#library(ggplot2)
library(ggtree)

setwd("C:\Users\lzhang27\Desktop\BayesianCompositionSelection\RealDataRelated")
load("combo.data.obj.RData")

myotu<-data.obj$otu.tab

write.table(row.names(myotu),"seqs.to.keep.txt",quote=F,row.names=F,col.names=F)

#################################
#######raw tree distance#####
################################
###But the taxonomies of some OTUs are missing
mytree<-data.obj$tree

dist.all<-vcv(mytree)
dist.share<-diag(dist.all)

corr<-vcv(mytree,corr=T)
quantile(corr)
f<-ecdf(corr)
##identify the 50% quantile of the distribution of correlations
f(0.5)
hist(corr,freq=T,main="Histogram of correlations",col="blue",xlab="Tree based correlations")
abline(v = 0.5, col="red", lwd=3, lty=2)
text(0.5, 5e+06,"89% quantile",cex=1.5)

write.table(corr,"combocorr.txt", row.names = F, col.names = F)

###We reclassify the representative sequences of each OTU into taxonomies
###with reference to Silva release 128

mytax<-read.csv("tax.levels.csv")
newtax<-mytax[,2:7]
row.names(newtax)<-mytax$X

###calculate shared distance
mymeta<-data.obj$meta.dat
myotu<-data.obj$otu.tab
myotu[1:6,1:6]
rowotu<-rowSums(myotu)
hist(rowotu)
quantile(rowotu)

otu.tab = otu_table(as.matrix(myotu),taxa_are_rows = TRUE)
tax.tab = tax_table(as.matrix(newtax))
meta.tab = sample_data(mymeta)

####check names
taxa_names(otu.tab)
sample_names(otu.tab)
taxa_names(tax.tab)
sample_names(tax.tab)
sample_names(meta.tab)
#####merge into an object of class phyloseq
physeq = phyloseq(otu.tab, tax.tab, meta.tab)
physeq.tree = merge_phyloseq(physeq, mytree)

#####Sample selection and taxa pruning
suppressPackageStartupMessages(library(phyloseq))
quantile(taxa_sums(physeq.tree))
physeq.tree= prune_taxa(taxa_sums(physeq.tree) > 10, physeq.tree)
levels(factor(tax_table(physeq.tree)[,"Genus"]))

newotu<-t(otu_table(physeq.tree))
newtax<-tax_table(physeq.tree)
newmeta<-sample_data(physeq.tree)
newtree<-phy_tree(physeq.tree)

##############################################################
####calculate the distance in the new truncated tree##########
##############################################################

dist.all<-vcv(newtree)
dist.share<-diag(dist.all)

sortcorr<-vcv(newtree,corr=T)
quantile(sortcorr)
f<-ecdf(sortcorr)
##identify the 50% quantile of the distribution of correlations
f(0.5)
setEPS()
postscript("Euclideancorrelation.eps")
hist(sortcorr,freq=T,col="blue",xlab="Euclidean correlation",cex.lab=1.5,cex.axis=1.5,font.lab=2,main=NULL)
abline(v = 0.5, col="red", lwd=3, lty=2)
text(0.52, 4e+05,"88% quantile",cex=1.5)
dev.off()

sum(colnames(sortcorr)!=colnames(newotu))
sum(colnames(sortcorr)!=rownames(newtax))

write.table(sortcorr,"sortcorr.csv",row.names = F, col.names = F,quote=F,sep=",")

# is.rooted(mytree)
# tr <- root(mytree, 1, resolve.root = TRUE)
# dist.mat<-cophenetic(tr)

###calculate distance between tips
dist.tip<-cophenetic(newtree)
rho=1.05
dist.mat<-exp(-2*rho*dist.tip)
quantile(dist.mat)
f<-ecdf(dist.mat)
##identify the 50% quantile of the distribution of correlations
f(0.5)
setEPS()
postscript("Exponentialcorrelation.eps")
hist(dist.mat,freq=T,col="blue",xlab="Exponential correlation",cex.lab=1.5,cex.axis=1.5,font.lab=2,main=NULL)
abline(v = 0.5, col="red", lwd=3, lty=2)
text(0.52, 4e+05,"88% quantile",cex=1.5)
dev.off()

dist.mat[dist.mat>0.9]<-0.9
dist.mat<-dist.mat+diag(rep(0.1,dim(dist.mat)[1]))
ttt<-solve(dist.mat)
ttt[ttt>1]<-1
ttt[ttt< -1]<- -1
f<-ecdf(ttt)
f(0)
setEPS()
postscript("Exponentialprecision.eps")
hist(ttt,freq=T,col="blue",xlab="Exponential precision",cex.lab=1.5,cex.axis=1.5,font.lab=2,main=NULL)
abline(v = 0, col="red", lwd=3, lty=2)
text(0, 4e+05,"53% quantile",cex=1.5)
dev.off()

#########################################################
####obtain relative abundances and save them
########################################################
load("combo.demo.RData")
bmi<-as.data.frame(bmi.c)
bmi$pid<-row.names(bmi)


newotu[newotu==0]<-0.5
newotu<-newotu/rowSums(newotu)
plot(rowSums(newotu))
newotu<-as.data.frame(newotu)
fulldata<-cbind(bmi,newotu)


otuid<-row.names(newtax)
newtax<-cbind(otuid,newtax)
#newtax[,7]<-gsub("_.*","",newtax[,7])
#unique(newtax[,7])
write.table(newtax,"tax.csv",row.names=F,col.names=T,quote=F,sep=",")
write.table(fulldata,"combodata.csv",row.names = F, col.names = T,quote=F,sep=",")
write.table(fulldata$bmi,"bmi.csv",row.names = F, col.names = F,quote=F,sep=",")
write.table(fulldata[,c(-1,-2)],"otu.csv",row.names = F, col.names = F,quote=F,sep=",")


###########filter the data by threshholding correlations###########
tt<-(sortcorr>0.9)
ttt<-colSums(tt)>500
plot(colSums(tt),ylab="Number of correlation more than 0.9")


hicorname<-names(ttt)[ttt]
important<-newdata[,hicorname]
barplot(colSums(important>0))

ignorename<-names(newdata[,c(-1,-2,-3)])[colSums(newdata[,c(-1,-2,-3)]>0)<4]
sortcorr[ignorename,]=rep(0,dim(sortcorr)[2])
sortcorr[,ignorename]=rep(0,dim(sortcorr)[2])
sortcorr[outer(rownames(sortcorr), colnames(sortcorr), "==")]=1
#sortcorr[sortcorr<0.5]<-0
hist(sortcorr,freq=T,main="Histogram of correlations",col="blue",xlab="Tree based correlations")
abline(v = 0.6, col="red", lwd=3, lty=2)
text(0.6, 4e+05,"80% quantile",cex=1.5)
quantile(sortcorr)
ecdf(sortcorr,0.6)
write.table(sortcorr,"ignoresortcorr.csv",row.names = F,col.names = F,quote=F,sep=",")
write.table(sortcorr,"ignoresortcorr.txt",row.names = F, col.names = F,quote=F,sep=" ")

##################################################
#############selected OTUs and their taxonomy#####
##################################################
select<-c(163, 168, 214, 268, 290, 357, 422, 470, 493, 513, 522, 605, 744, 770,834,844, 899, 1048,
1149,1234,1598,1922,1982,2096,2307,2621, 3081, 3231,3239,3299)

tax<-read.table("tax.levels.csv",header=TRUE,sep=",")
new<-tax[select+1,]
write.table(new,"selected.csv",row.names = F,quote=F,sep=",")
