###R script for randomly removing reads from simulated data
##05/20/20

library(phylotools)

#load data from process_sum_polyploid_github.R
setwd("~/Desktop/ipcoal-sims/poly_mini")
load("~/Desktop/ipcoal-sims/poly_mini/sim_poly.RData")

#what proportion of reads do you want to remove
rm <- 0.50
reads_rm <- 2000*0.5

#loop through all
#set up position of random samples
cutsall<-vector("list", 50)

for(i in 1:50){
set.seed(385804+i)
cutsall[[i]]<-  base::sample.int(2000, reads_rm, replace = F)
}

#replace those positions with NA
for(i in 1:50){
  physcomb_miss[[i]]$seq.text[cutsall[[i]]]<-NA
}

#rename to missing
fa_names_miss50<-paste(names(physa_ind), "_miss50.fa",sep="")

#write output
for(i in 1:50){
  dat2fasta(physcomb_miss[[i]], outfile = fa_names_miss50[i])
}

#output files run through mash and alignment
