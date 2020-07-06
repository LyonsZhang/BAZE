closeAllConnections()
rm(list=ls())
# install.packages("devtools")
#devtools::install_github("USCBiostats/rphyloxml")

library(ape)
library(phylobase)
library("phyloseq")
#packageVersion("phyloseq")

#library(phyloseq)
set.seed(12)

#########load taxonomy data######
setwd("C:\Users\lzhang27\Desktop\BayesianCompositionSelection\RealDataRelated")
tax.clean<-read.table("tax.levels.csv",header=TRUE,sep=",",row.names = 1)
tax.clean<-tax.clean[,names(tax.clean)[-8]]

species<-as.character(tax.clean$Species)

temp<-tax.clean$Genus
temp<-sub("(^[^_]+_[^_]+)_(.*)$","\\1",temp)
# #temp<-gsub("_\\w+", "", mytree$tip.label)
temp<-gsub("-\\w+","",temp)
genus<-gsub("_.*","",temp)
genus[! (genus %in% c("Ruminococcaceae","Ruminiclostridium","Roseburia",
             "Lachnospiraceae","Blautia","Bacteroides",
             "Alistipes","Clostridiales","Christensenellaceae",
             "Erysipelotrichaceae","Faecalibacterium","Lachnoclostridium",
             "Coprococcus","Acidaminococcus","Parasutterella","Allisonella"))]<-NA
#"Bacteroidales","Prevotella", "Parabacteroides"

temp<-tax.clean$Family
temp<-sub("(^[^_]+_[^_]+)_(.*)$","\\1",temp)
# #temp<-gsub("_\\w+", "", mytree$tip.label)
temp<-gsub("-\\w+","",temp)
family<-gsub("_.*","",temp)

temp<-tax.clean$Order
order<-gsub("_.*","",temp)
temp<-tax.clean$Class
class<-gsub("_.*","",temp)
temp<-tax.clean$Phylum
phylum<-gsub("_.*","",temp)
kingdom<-tax.clean$Kingdom


# typeof(tax.clean)
# kingdom<-paste("k",tax.clean$Kingdom,sep="_")
# phylum<-paste("p",tax.clean$Phylum,sep="_")
# class<-paste("c",tax.clean$Class,sep="_")
# order<-paste("o",tax.clean$Order,sep="_")
# family<-paste("f",tax.clean$Family,sep="_")
# genus<-paste("g",tax.clean$Genus,sep="_")
# species<-paste("s",tax.clean$Species,sep="_")
# tax.clean$taxonomy<-paste(kingdom,phylum,class,order,family,genus,species, sep=";")
# tax.clean$taxonomy<-gsub("Woesearchaeota_","Woesearchaeota",tax.clean$taxonomy)
# tax.clean$taxonomy<-gsub("SR1_","SR1",tax.clean$taxonomy)
# tax.clean$taxonomy<-gsub("TM6_","TM6",tax.clean$taxonomy)
tax<-data.frame(cbind(phylum,class,order,family,genus,species))
rownames(tax)<-rownames(tax.clean)

load("combo.data.obj.RData")
mymeta<-data.obj$meta.dat
myotu<-data.obj$otu.tab
#mytree<-read_tree("total.10.tree")
mytree<-data.obj$tree
meta.tab = sample_data(mymeta)
otu.tab = otu_table(as.matrix(myotu),taxa_are_rows = TRUE)
tax.tab = tax_table(as.matrix(tax))
physeq = phyloseq(otu.tab, tax.tab, meta.tab, mytree)
physeq.tree= prune_taxa(taxa_sums(physeq) > 10, physeq)

temtax<-data.frame(tax_table(physeq.tree))
newtax<-data.frame(temtax[,c("phylum","genus")])
newtax$phylum[which(!(newtax$phylum=="Firmicutes" | newtax$phylum=="Bacteroidetes"))]<-NA
newtree<-phy_tree(physeq.tree)


library(geiger)
tphy=phylo.lookup(newtax, ncores=4)
finaltree<-nodelabel.phylo(newtree, newtax, strict=F, ncores=4)

grep("Firmicutes",finaltree$node.label,value=T)


#########load tree#######
mytree$tip.label
phy=mytree
taxonomy=tax

mytree$node.label<-as.character(1:tree$Nnode)

library(rphyloxml)
phylo4.dat<-phylo4(tree)
ancestors(phylo4.dat, c(1,2), type = "all")


tree$tip.label
tree$node.label
str(tree)
tree$edge
tree$node.label
plot(tree)
nodelabels
tiplabels
z <- write_phyloxml(finaltree)

xml2::write_xml(z, "mynicetree.xml")

  
taxn<-temtax
taxn$tip.label<-rownames(temtax)
taxn$tip.label<-as.character(taxn$tip.label)
tip.label<-finaltree$tip.label
tips<-as.data.frame(tip.label,make.names = F)
tips$tip.label<-as.character(tips$tip.label)
library(plyr)
tem<-join(tips,taxn,by="tip.label")

#######find common ancester of phylum
library(phytools)
temptree<-newtree
temptree$node.label<-as.character(1:length(temptree$node.label))
T1<-(tem$phylum=="Firmicutes")
T1[is.na(T1)]<-FALSE
t1<-findMRCA(tree=temptree,tips=tem$tip.label[T1],type="node")
T2<-tem$phylum=="Bacteroidetes"
T2[is.na(T2)]<-FALSE
t2<-findMRCA(tree=temptree,tips=tem$tip.label[T2],type="node")

library(ggtree)

finaltree$tip.label<-as.character(tem$genus)
temp<-finaltree$tip.label
groupInfo <- split(finaltree$tip.label, temp)
newftree <- groupOTU(finaltree, groupInfo)
indx<-1:length(temp)
indx1<-indx[!duplicated(temp)]
indx0<-indx[duplicated(temp)]
temp[indx0]<-NA
newftree$tip.label<-temp



ggtree(newftree, aes(color=group),layout='circular',branch.length="none",linetype="solid")+
  geom_tiplab(size=0, aes(angle=angle),align=T,offset=3,linesize=6,linetype="1F")+
  geom_nodelab(size=7,aes(angle=angle),linetype="blank",align=T,nudge_x=14)+
  geom_hilight(node=1824,fill='aquamarine4',alpha=.3)+
  geom_hilight(node=3063,fill="darkgoldenrod4",alpha=.3)
ggsave("grouptreeplot.eps",device=cairo_ps,width=15,height=15)

#   theme(legend.position = "right") 
#   
#   
#   
#   
#   geom_hilight(node=2396,fill="darkgreen",alpha=.4)+
#   
#   geom_hilight(node=1824,fill="darkgreen",alpha=.4)+
#   
#   geom_hilight(node=2090,fill="darkgreen",alpha=.4)
#   
#   geom_hilight(node=1826,fill="darkgreen",alpha=.4)
#   
#   
#   
#   geom_hilight(node=which(newftree$node.label=="\"Bacteroides\""),fill="darkgreen",alpha=.6)
#   geom_nodelab(size=5,aes(angle=angle),linetype="blank",align=T,nudge_x=13)
# 
#   geom_hilight(node=t1,fill="grey",alpha=.5)+
#   geom_hilight(node=t2,fill="grey",alpha=.5)
#   
# 
# theme(legend.position = "right") 
#   geom_hilight(node=which(newftree$node.label=="\"Bacteroides\""),fill="darkgreen",alpha=.6)
# 
#   
#   geom_cladelabel(node=5,label="Streptococcus",offset=10)+
#   geom_cladelabel(node=1332,label="Bacteroides",offset=5)
# +
#   geom_hilight(node=1332,fill="black",alpha=.6)
# 
#  
#   geom_nodelab(size=4, aes(angle=angle),color="black",linetype="blank")+
#   geom_hilight(node="Bacteroides",fill="darkgreen",alpha=.6)
#   
  
ggsave("grouptreeplot.eps",device=cairo_ps,width=15,height=15)

