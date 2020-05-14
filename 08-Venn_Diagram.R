### Script prepared by Carolina S. Carvalho, Jeronymo Dalapicolla, Luciana C. Resende-Moreira, Jamille C. Veiga, and Rodolfo Jaffé ###

#------------------------------------------------------------------------------
#                               PRE-ANALYSIS 
#------------------------------------------------------------------------------
##1. GOALS FOR THIS STEP:
#A. PROCESSING DATA AND CREATE DIAGRAMS FOR THE NUMBER ADAPTATIVE LOCI

##2. REMOVE ANY OBJECT OR FUNCTION IN THE ENVIRONMENT:
rm(list=ls())

##3. INSTALL AND LOAD THE PACKAGES

# devtools::install_github("cran/VennDiagram", force=T)
# devtools::install_github("tidyverse/tidyverse", force=T)
# devtools::install_github("cran/seqinr", force=T)

library(VennDiagram)
library(tidyverse)
library(seqinr)

##4. INPUTS FOR THIS STEP:
#A. FASTA FILES WITH GEA RESULTS USING RDA AND LFMM ANALYSES.


#------------------------------------------------------------------------------
#                            Load Fasta Files 
#------------------------------------------------------------------------------
## Load RDA results
rda_gea <-read.fasta(file="adapt_var_mapping/GEA-RDA-Mahalanobis/SEQ_gea_RDA_mahalanobis_ipomoea.fasta")
length(rda_gea) ## 408

## Load LFMM results
lfmm_gea <-read.fasta(file="adapt_var_mapping/EAA_LFMM/SEQ_total_lfmm_ipomoea.fasta")
length(lfmm_gea) ## 112

BM1 <- rda_gea
BM2 <- lfmm_gea


#------------------------------------------------------------------------------
#                                Create the Diagram 
#------------------------------------------------------------------------------

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

## Save in pdf

## END
