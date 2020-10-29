###General R script for calculating genetic distance from aligned reads
##performed on the MSU HPCC cluster
##06/05/20

#load libraries
library(adegenet)
library(ape)
library(vcfR)
library(dplyr)

#read output vcf
outall<-read.vcfR("outall.vcf.gz")
vcf_genind <- vcfR2genind(outall)
genind<-vcf_genind

#add sample names as pops
popvector<-attr(genind$tab, "dimnames")[1][[1]]
popvector2<-as.data.frame(substr(popvector, 5,(nchar(popvector)-14)))
colnames(popvector2)<-"pops"
genind$pop<-as.factor(popvector2$pops)

#convert to genpop
genpop<-genind2genpop(genind) 

#calculate nei's distance 
dist_nei <- dist.genpop(genpop, method=1)
save(dist_nei, file="dist_poly_time.rda")