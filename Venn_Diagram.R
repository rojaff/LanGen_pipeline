
#####################
#### Venn diagram ###
#####################

library(VennDiagram)
library(tidyverse)

######################################
############# DATA ENTRY ############# 
######################################

library(seqinr)

## The Processed Data and Diagrams (NUMBER ADAPTATIVE LOCI) ##

rda_gea <-read.fasta(file="adapt_var_mapping/GEA-RDA-Mahalanobis/SEQ_gea_RDA_mahalanobis_ipomoea.fasta")
length(rda_gea) ## 408

lfmm_gea <-read.fasta(file="adapt_var_mapping/EAA_LFMM/SEQ_total_lfmm_ipomoea.fasta")
length(lfmm_gea) ## 112


BM1 <- rda_gea
BM2 <- lfmm_gea

## Creating a Venn Diagram with FOUR circles

w1 <- venn.diagram(list(RDA=BM1,LFMM=BM2), lty = c("blank", "blank"),
                  fill = c("red", "light blue"),
                  alpha = c(0.5, 0.5), cat.cex = 1.2, cex= 1.5,  cat.pos = 0, scaled = F,
                  filename=NULL)

## The default plot

grid.newpage()
grid.draw(w1)

rda_gea
lfmm_gea

### save in pdf