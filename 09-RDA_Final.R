### Script prepared by Carolina S. Carvalho, Jeronymo Dalapicolla, Luciana C. Resende-Moreira, Jamille C. Veiga, and Rodolfo Jaffé ###

#------------------------------------------------------------------------------
#                               PRE-ANALYSIS 
#------------------------------------------------------------------------------
##1. GOALS FOR THIS STEP:
#A. PERFORM VISUAL ASSESSMENT OF GENOTYPE-ENVIRONMENT ASSOCIATIONS USING REDUNDANCY ANALYSES (RDA) ON THE CANDIDATE LOCI

##2. REMOVE ANY OBJECT OR FUNCTION IN THE ENVIRONMENT:
rm(list=ls())

##3. INSTALL AND LOAD THE PACKAGES

# devtools::install_github("bcm-uga/LEA", force=T)
# devtools::install_github("nspope/r2vcftools", force=T)
#devtools::install_github("cran/usdm", force=T)
# devtools::install_github("vegandevs/vegan", force=T)
# devtools::install_github("cran/psych", force=T)
# devtools::install_github("cran/seqinr", force=T)
# devtools::install_github("cran/robust", force=T)
# devtools::install_github("sdray/ade4", force=T)
# devtools::install_github("sdray/adespatial", force=T)
# devtools::install_github("r-spatial/spdep", force=T)
#install.packages("spacemakeR", #repos="http://R-Forge.R-#project.org")

library(LEA)
library(r2vcftools)
library(vegan)    
library(usdm)
library(spdep)
library(adespatial)
library(spacemakeR) 
library(ade4)
library(psych)
library(seqinr)
library(robust)  

##4. INPUTS FOR THIS STEP:
#A. THE FILE ".VCF" FOR ADAPATATION ANALYSES CREATED IN FILTERING STEP 1, NAMED AS "_filtered_within_ld_test2.vcf".
#B. GENOTYPE FILE WITHOUT MISSING DATA CREATED IN FILTERING STEP 1, NAMED AS "_filtered_within_ld_imputed.geno".
#C. GEA RESULTS USING RDA AND LFMM ANALYSES FROM STEPS 6 AND 7.
#D. ENVIRONMENTAL INFORMATION FILTERED AND SELECTED IN STEP 3.
#E. POPULATIONS ID PER INDIVIDUAL FROM STEP 2. 

##5. LOAD FUNCTIONS TO BE USED ON THIS STEP.
VCFsummary <- function(snps){
  Nid <- nrow(snps@meta)
  Nsnps <- length(snps@site_id)
  cat(paste(Nid, "individuals and", Nsnps, "SNPs."), sep="\n")
}


#------------------------------------------------------------------------------
#                        Load Genotypes Files and GEA Results 
#------------------------------------------------------------------------------
#Load snps data set to take the snps and individuals ID
snps <-  vcfLink("vcf/ipomoea_filtered_within_ld_test2.vcf", overwriteID=T)
snps@meta
VCFsummary(snps) ## 115 individuals and 17025 SNPs.

## Read geno file
gen.imp <- read.geno("vcf/ipomoea_filtered_within_ld_imputed.geno")
gen.imp[1:10,1:10]
colnames(gen.imp) <- snps@site_id
rownames(gen.imp) <- snps@sample_id
dim(gen.imp)

## Load RDA and LFMM results - genotype vs enviroment
rda_snps <- vcfLink("adapt_var_mapping/GEA-RDA-Mahalanobis/snps_candidate_gea_RDA_mahalanobis_ipomoea.vcf", overwriteID=F) # here overwrite need to be F, otherwise we will loss the snps name
lfmm_snps <- vcfLink("adapt_var_mapping/EAA_LFMM/snps_candidate_lfmm_ipomoea.vcf", overwriteID=F)

all_candidate <- unique(c(rda_snps@site_id,lfmm_snps@site_id))
all_candidate_gen <- gen.imp[,all_candidate]
dim(all_candidate_gen)


#------------------------------------------------------------------------------
#                           Load Environmental Dataset 
#------------------------------------------------------------------------------
## Load environmental data set
dados_env <- read.table("adapt_var_mapping/Dados_ambientais/climatic_data_ipomoea.txt", head=T)
head(dados_env)
str(dados_env)

## Set population ID from step 2 
pop_lea <- read.csv("metadados/info_serras_115samples.csv", head=T)
head(pop_lea)

## Confirm that genotypes and environmental data are in the same order
data.frame(rownames(gen.imp),as.character(dados_env$ID))

## Scale environmental dataset
pred_scale <- scale(dados_env[,-c(1:3)]) ## remove information of lat, long, location, etc
rownames(pred_scale) <-dados_env[,1]
head(pred_scale)


#------------------------------------------------------------------------------
#                               Perform RDA Analysis 
#------------------------------------------------------------------------------
## Run RDA
rda <- rda(all_candidate_gen ~ ., data=as.data.frame(pred_scale))  
summary(rda)        
RsquareAdj(rda)#Our constrained ordination explains 17% of the variation

summary(eigenvals(rda, model = "constrained")) #The eigenvalues for the constrained axes reflect the variance explained by each canonical axis
screeplot(rda) 

## Here, we can see that the first two constrained axes explain most of the variance. The screeplot provides an informal (and quick) way to determine how many constrained axes to include when we search for candidate SNPs.
signif.full <- anova(rda)  #  check our RDA model for significance using formal tests.
signif.axis <- anova(rda, by="axis") #  check our for significance of rda axis using formal tests.
signif.terms <-anova(rda, by="terms")


#------------------------------------------------------------------------------
#                                  Plot Results 
#------------------------------------------------------------------------------
## Plots
plot(rda, scaling=3)          # default is axes 1 and 2

eco<-as.factor(pop_lea$Canga) ## color by populations
bg <- c("#ff7f00","#1f78b4","#ffff33", "green") # 4 nice colors for our ecotypes
plot(rda, type="n", scaling=3)
points(rda, display="species", pch=20, cex=0.7, col="gray32", scaling=3) # the SNPs
points(rda, display="sites", pch=21, cex=1.3, col="gray32", scaling=3, bg=bg[eco]) # the individuos
text(rda, scaling=3, display="bp", col="#0868ac", cex=1) # the predictors
legend("bottomright", legend=levels(eco), bty="n", col="gray32", pch=21, cex=1, pt.bg=bg)

## END
