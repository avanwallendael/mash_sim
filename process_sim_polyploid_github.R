###R script for processing simulated loci for polyploids
##05/20/20

#load libraries
library(phylotools)
library(phangorn)
library(adegenet)
library(ape)
library(vcfR)

#input files are output from mash_sim_github.ipynb

setwd("~/Desktop/ipcoal-sims/poly_mini")

#make vectors of filenames
phyiles<-list.files(pattern=".phy")
philesa<-list.files(pattern = "mini1")
philesb<-list.files(pattern = "mini2")

#read in genomes a & b
physa<-vector("list", 1000)
for (i in 1:1000){
  physa[[i]]<-read.phylip(philesa[[i]])
}

physb<-vector("list", 1000)
for (i in 1:1000){
  physb[[i]]<-read.phylip(philesb[[i]])
}

#name them
names(physa)<-paste("loca",as.character(1000:1999), sep="_")
names(physb)<-paste("locb",as.character(1000:1999), sep="_")

locnamesa<-paste("loca",as.character(1000:1999), sep="_")
locnamesb<-paste("locb",as.character(1000:1999), sep="_")

#add names as a column
for(i in 1:1000){physa[[i]]$locname<-locnamesa[i]}
for(i in 1:1000){physb[[i]]$locname<-locnamesb[i]}

#bind all
physlonga<-do.call(rbind, physa)
physlongb<-do.call(rbind, physb)

#now split into list with seq.name as splitter
physa_ind<-split(physlonga, physlonga$seq.name)
physb_ind<-split(physlongb, physlongb$seq.name)

#clean to just seqs
for(i in 1:50){
  physa_ind[[i]]<-physa_ind[[i]][,c("locname", "seq.text")]
  colnames(physa_ind[[i]])<-c("seq.name", "seq.text")
}
for(i in 1:50){
  physb_ind[[i]]<-physb_ind[[i]][,c("locname", "seq.text")]
  colnames(physb_ind[[i]])<-c("seq.name", "seq.text")
}

#combine subgenomes
physcomb_ind<-vector("list", 50)
for(i in 1:50){
  physcomb_ind[[i]]<-rbind(physa_ind[[i]], physb_ind[[i]])
}

fa_names<-paste(names(physa_ind), ".fa",sep="")

#write as fasta files
for(i in 1:50){
  dat2fasta(physcomb_ind[[i]], outfile = fa_names[i])
}

#output is processed by mash and alignment tools in bash scripts. 
# examples in capsella_mash.bash & polymiss_github.bash