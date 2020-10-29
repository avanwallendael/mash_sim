###General R script for processing and visualizing mash outputs
##10/25/20

#load libraries
library(readr)
library(tidyverse)
library(vegan)
library(adegenet)
library("maps")


setwd("~/Downloads")

#read in mash output
capsellatrim <- read_delim("tbl1_capsellatrim_tailcrop.tab", 
                         "\t", escape_double = FALSE, col_names = FALSE, 
                         trim_ws = TRUE)

#function to clean dists
clean_dist<-function(mashdist){
  #rename cols
  colnames(mashdist)<-c("f1", "f2","dist","num", "frac")
  #rm extra cols
  mashdist$num<-NULL
  mashdist$frac<-NULL
  #spread to matrix
  in_wide<-spread(mashdist, key=f2, value = dist)
  #set rownames for mat
  row.names(imp_wide)<-in_wide$f1
  #rm col
  in_wide$f1<-NULL
  #make into numeric matrix
  in_mat<-data.matrix(in_wide)
  #make dist object
  in_dist<-as.dist(in_mat)
  in_dist
}

dists<-clean_dist(capsellatrim)

#read in metadata from SRA
phenolist <- read_delim("~/Downloads/SraRunTable (29).txt", ",", escape_double = FALSE, trim_ws = TRUE)

#perform principal coordinate analysis
cailliez<-cailliez(dists)
pcoa<-dudi.pco(cailliez, scannf = T, full=T)

#pull pcs
PCs<-data.frame(PC1=pcoa$tab[,1],
                PC2=pcoa$tab[,2],
                PC3=pcoa$tab[,3],
                PC4=pcoa$tab[,4],
                PC5=pcoa$tab[,5],
                PC6=pcoa$tab[,6],
                file=substr(attr(mds$points, 'dimnames')[[1]], 1, 10), stringsAsFactors = F)
colnames(PCs)[7]<-"Run"

#add in metadata
PCall<-inner_join(PCs, phenolist, by="Run")

#save(PCall, file = "PCall_mash_oct25.rda")


###Visualization####

#clean pop names
pops <- strsplit(PCall$Isolation_Source, "_")
pops2 <- NA
for (i in 1:length(pops)){
  pops2[i] <- pops[[i]][1]
}

PCall$pops<-pops2

#read in metadata from paper supplement
pop_loc_capsella <- read_delim("~/Downloads/pop_loc_capsella.txt", 
                               "\t", escape_double = FALSE, col_names = FALSE, 
                               trim_ws = TRUE)
unique(pop_loc_capsella$X2)
unique(PCall$pops)

#replace erroneous pop names
PCall$pops[PCall$pops %in% c("6")]<-"FR6"
PCall$pops[PCall$pops %in% c("53")]<-"SP53"
PCall$pops[PCall$pops %in% c("39")]<-"IT39"
PCall$pops[PCall$pops %in% c("22")]<-"IT22"
PCall$pops[PCall$pops %in% c("IRRU2","IRRU3")]<-"IRRU"
PCall$pops[PCall$pops %in% c("JO56","JO59")]<-"JO"
PCall$pops[PCall$pops %in% c("OBL")]<-"OBL-RU5"
PCall$pops[PCall$pops %in% c("SABO10", "SABO4",  "SABO5",  "SABO6",  "SABO7",  "SABO9")]<-"SABO"
PCall$pops[PCall$pops %in% c("SE42","SE43")]<-"SE4x"
PCall$pops[PCall$pops %in% c("SY61",   "SY64",   "SY67" ,  "SY68" ,  "SY69")]<-"SY6x" 
PCall$pops[PCall$pops %in% c("SY70")]<-"SY7x" 
PCall$pops[PCall$pops %in% c( "TR71"  , "TR73" ,  "TR75" ,  "TR79"  , "TR83"  )]<-"TR" 
PCall$pops[PCall$pops %in% c( "VORU0",  "VORU1" , "VORU2" , "VORU3")]<-"VORU" 
PCall$pops[PCall$pops %in% c("CHS")]<-"CSH" 
PCall$pops[PCall$pops %in% c("SY70")]<-"SY7x" 

pop_loc_capsella$pops<-pop_loc_capsella$X2

named_pops <- left_join(PCall, pop_loc_capsella, by = "pops") %>% filter(!is.na(X1))

#assign countries to regions for visualization
countrykey <- data.frame(country = unique(named_pops$geo_loc_name),
                         region = c("Europe", "Europe","Europe","Africa","Russia",rep("Europe", 2), "ME", "ME", rep("Europe", 3),"ME", "ME", 
                                    "USA",  "China", "Taiwan"))
named_pops$region <- countrykey$region[match(named_pops$geo_loc_name, countrykey$country)]

#PCoA visualization
ggplot(named_pops)+
  geom_jitter(aes(x=PC1, y=PC2, col=region), size=2, width = 0.005)+
  geom_hline(aes(yintercept=-0.001))+
  geom_vline(aes(xintercept=0))+
  theme_classic()

#assign clusters from PC1&2
named_pops$group <- "Europe1"
named_pops$group[named_pops$PC1<0&abs(named_pops$PC2)<0.0025]<-"Asia1"
named_pops$group[named_pops$PC1<0&abs(named_pops$PC2)>0.007]<-"Asia2"
named_pops$group[named_pops$PC1>0&named_pops$PC2>(-0.001)]<-"Europe2"
named_pops$group[named_pops$PC1>0&named_pops$PC2<(-0.001)]<-"ME"


#download worldmap
worlddata <- map_data("world")

#plot worldmap with assigned pops
ggplot(worlddata)+
  geom_polygon(aes(x=long, y=lat, group=group), fill="grey90", col="black")+
  geom_jitter(data=named_pops, aes(x=X5, y=X4, col=group), size=4)+
  scale_color_manual(values=c("chartreuse3","forestgreen", "red", "firebrick", "royalblue4"))+
  coord_cartesian(xlim=c(-15,165), ylim=c(20,70))+
  theme_classic()
ggsave("capsella_map_main.png", height = 3, width=5.5)

#plot USA portion
ggplot(worlddata)+
  geom_polygon(aes(x=long, y=lat, group=group), fill="grey90", col="black")+
  geom_jitter(data=named_pops, aes(x=X5, y=X4, col=group), size=8)+
  scale_color_manual(values=c("chartreuse3","forestgreen", "red", "firebrick", "royalblue4"), guide="none")+
  coord_cartesian(xlim=c(-125,-70), ylim=c(25,50))+
  theme_classic()
ggsave("capsella_map_US.png", height = 3, width=5.5)

