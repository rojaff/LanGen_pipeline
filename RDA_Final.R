## This code will be used to run a RDA using only the candidate SNPs

# Load packages -----------------------------------------------------------
library(LEA)
library(r2vcftools)
library(vegan)    # Used to run RDA
library(usdm)
library(spdep)
library(adespatial)
#library(spacemakeR) ## Nao instalei
library(ade4)
library(psych)
library(seqinr)
library(robust)  

# Load functions ----------------------------------------------------------

##Load functions that will be used throughout the code

### Summary stats
VCFsummary <- function(snps){
  Nid <- nrow(snps@meta)
  Nsnps <- length(snps@site_id)
  cat(paste(Nid, "individuals and", Nsnps, "SNPs."), sep="\n")
}

# RDA of GEA candidates SNPs --------------------------------------------------

## genotype data and impute data file

#Load snps data set to take the snps and individuals ID

snps <-  vcfLink("vcf/ipomoea_filtered_within_ld_test2.vcf", overwriteID=T)
snps@meta
VCFsummary(snps) ## 115 individuals and 17025 SNPs.

## read geno file

gen.imp <- read.geno("vcf/ipomoea_filtered_within_ld_imputed.geno")
gen.imp[1:10,1:10]
colnames(gen.imp) <- snps@site_id
rownames(gen.imp) <- snps@sample_id
dim(gen.imp)

# load rda and lfmm candidate snps - genotype vs enviroment
rda_snps <- vcfLink("adapt_var_mapping/GEA-RDA-Mahalanobis/snps_candidate_gea_RDA_mahalanobis_ipomoea.vcf", overwriteID=F) # here overwrite need to be F, otherwise we will loss the snps name
lfmm_snps <- vcfLink("adapt_var_mapping/EAA_LFMM/snps_candidate_lfmm_ipomoea.vcf", overwriteID=F)

all_candidate <- unique(c(rda_snps@site_id,lfmm_snps@site_id))
all_candidate_gen <- gen.imp[,all_candidate]
dim(all_candidate_gen)

### load environmental data set

dados_env <- read.table("adapt_var_mapping/Dados_ambientais/climatic_data_ipomoea.txt", head=T)
head(dados_env)
str(dados_env)

# population dataset 

pop_lea <- read.csv("metadados/info_serras_115samples.csv", head=T)
head(pop_lea)


# Confirm that genotypes and environmental data are in the same order

data.frame(rownames(gen.imp),as.character(dados_env$ID))

## scale environmental dataset

pred_scale <- scale(dados_env[,-c(1:3)]) ## remove information of lat, long, location, etc
rownames(pred_scale) <-dados_env[,1]
head(pred_scale)


### Run RDA

rda <- rda(all_candidate_gen ~ ., data=as.data.frame(pred_scale))  
summary(rda)        
RsquareAdj(rda)#Our constrained ordination explains 17% of the variation

summary(eigenvals(rda, model = "constrained")) #The eigenvalues for the constrained axes reflect the variance explained by each canonical axis
screeplot(rda) 

#Here, we can see that the first two constrained axes explain most of the variance. The screeplot provides an informal (and quick) way to determine how many constrained axes to include when we search for candidate SNPs.

signif.full <- anova(rda)  #  check our RDA model for significance using formal tests.

signif.axis <- anova(rda, by="axis") #  check our for significance of rda axis using formal tests.

signif.terms <-anova(rda, by="terms")

## plots

plot(rda, scaling=3)          # default is axes 1 and 2

eco<-as.factor(pop_lea$Canga) ## color by populations
bg <- c("#ff7f00","#1f78b4","#ffff33", "green") # 4 nice colors for our ecotypes
plot(rda, type="n", scaling=3)
points(rda, display="species", pch=20, cex=0.7, col="gray32", scaling=3) # the SNPs
points(rda, display="sites", pch=21, cex=1.3, col="gray32", scaling=3, bg=bg[eco]) # the individuos
text(rda, scaling=3, display="bp", col="#0868ac", cex=1) # the predictors
legend("bottomright", legend=levels(eco), bty="n", col="gray32", pch=21, cex=1, pt.bg=bg)

