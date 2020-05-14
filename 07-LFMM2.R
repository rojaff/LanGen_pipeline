### Script prepared by Carolina S. Carvalho, Jeronymo Dalapicolla, Luciana C. Resende-Moreira, Jamille C. Veiga, and Rodolfo Jaffé ###

#------------------------------------------------------------------------------
#                               PRE-ANALYSIS 
#------------------------------------------------------------------------------

##1. GOALS FOR THIS STEP:
#A. PERFORM LFMM2 ANALYSIS TO IDENTIFY LOCI UNDER SELECTION (OUTLIER DETECTION)

##2. REMOVE ANY OBJECT OR FUNCTION IN THE ENVIRONMENT:
rm(list=ls())

##3. INSTALL AND LOAD THE PACKAGES

# devtools::install_github("bcm-uga/LEA", force=T)
# devtools::install_github("nspope/r2vcftools", force=T)
# devtools::install_github("bcm-uga/lfmm", force=T)
# devtools::install_github("cran/usdm", force=T)
# devtools::install_github("vegandevs/vegan", force=T)
# devtools::install_github("cran/seqinr", force=T)
# if (!# requireNamespace("BiocManager", # quietly = TRUE))
#    install.packages("BiocManager")
# BiocManager::install("qvalue")

library(r2vcftools)
library(LEA)
library(lfmm)
library(usdm)
library(vegan)
library(seqinr)
library(qvalue)

##4. INPUTS FOR THIS STEP:
#A. THE FILE ".VCF" FOR ADAPATATION ANALYSES CREATED IN FILTERING STEP 1, NAMED AS "_filtered_within_ld_test2.vcf".
#B. GENOTYPE FILE WITHOUT MISSING DATA CREATED IN FILTERING STEP 1, NAMED AS "_filtered_within_ld_imputed.geno".
#C. REFERENCE SEQUENCES IN FASTA FORMAT. YOU CAN DOWNLOAD AS EXAMPLE "Icavalcantei.fasta" FROM FIGSHARE: https://doi.org/10.6084/m9.figshare.6100004.v1
#D. ENVIRONMENTAL INFORMATION FILTERED AND SELECTED IN STEP 3.

##5. LOAD FUNCTIONS TO BE USED ON THIS STEP.
VCFsummary <- function(snps){
  Nid <- nrow(snps@meta)
  Nsnps <- length(snps@site_id)
  cat(paste(Nid, "individuals and", Nsnps, "SNPs."), sep="\n")
}


#------------------------------------------------------------------------------
#                            Load Datasets 
#------------------------------------------------------------------------------
## Load snps data set to take the snps and individuals ID
snps <-  vcfLink("vcf/ipomoea_filtered_within_ld_test2.vcf", overwriteID=T)
snps@meta
VCFsummary(snps) ## 115 individuals and 17025 SNPs.

## Read geno file
ipomoea_geno <- read.geno("vcf/ipomoea_filtered_within_ld_imputed.geno")
ipomoea_geno[1:10,1:10]
colnames(ipomoea_geno) <- snps@site_id
rownames(ipomoea_geno) <- snps@sample_id
dim(ipomoea_geno)

## Number of missing data
sum(is.na(ipomoea_geno))

## Load environmental data set
dados_env <- read.table("adapt_var_mapping/Dados_ambientais/climatic_data_ipomoea.txt", head=T)
head(dados_env)
str(dados_env)

## Check if individuals in genotype and environmental dataset are in the same order
data.frame(rownames(ipomoea_geno),as.character(dados_env$ID))

## Check the correlation among enviromental variables
vif(dados_env[,-c(1:3)])

## scale environmental dataset
pred_scale <- scale(dados_env[,-c(1:3)]) ## remove information of lat, long, location, etc
rownames(pred_scale) <-dados_env[,1]
head(pred_scale)


#------------------------------------------------------------------------------
#                            LFMM Analysis 
#------------------------------------------------------------------------------
## Carry out the lfmm using the package lfmm
mod.lfmm <- lfmm_ridge(Y = ipomoea_geno, 
                       X = pred_scale, 
                       K = 1)  ## K is the number of populations

## Performs association testing using the fitted model:
pv <- lfmm_test(Y = ipomoea_geno, 
                X = pred_scale, 
                lfmm = mod.lfmm,
                calibrate="gif")


#------------------------------------------------------------------------------
#                             Outlier Detection 
#------------------------------------------------------------------------------
## Estimate adjusted p-values based on gif
pv$gif# the GIFs for the predictors

# Reminder:
# GIF of 1=well calibrated, >1=liberal (too many small p-values), <1=conservative (too few small p-values)
# Note: GIFs > 2 indicate poor test calibration; try increasing K...and be skeptical!

## Estimate adjusted p-values
pvalues <- pv$pvalue 
zs <- pv$score

p.PC1 <- pchisq(zs[,1]^2/pv$gif[1], df = 1, lower = FALSE)  ## too conservative
adj.PC1 <- pchisq(zs[,1]^2/1.2, df = 1, lower = FALSE) 
hist(p.PC1, main="Histogram of p-values PC1")
hist(adj.PC1, main="Histogram of adjusted p-values PC1")

p.PC2 <- pchisq(zs[,2]^2/pv$gif[2], df = 1, lower = FALSE) 
adj.PC2 <- pchisq(zs[,2]^2/0.9, df = 1, lower = FALSE) 
hist(p.PC2, main="Histogram of p-values PC2")
hist(adj.PC2, main="Histogram of adjusted p-values PC2")

p.PC3 <- pchisq(zs[,3]^2/pv$gif[3], df = 1, lower = FALSE) 
adj.PC3 <- pchisq(zs[,3]^2/2.3, df = 1, lower = FALSE) 
hist(p.PC3, main="Histogram of p-values PC3")
hist(adj.PC3, main="Histogram of adjusted p-values PC3")

qv_adj.PC1 <- which(qvalue::qvalue(adj.PC1, fdr=0.05)$signif)
length(qv_adj.PC1)

qv_adj.PC2 <- which(qvalue::qvalue(adj.PC2, fdr=0.05)$signif)
length(qv_adj.PC2)

qv_adj.PC3 <- which(qvalue::qvalue(adj.PC3, fdr=0.05)$signif)
length(qv_adj.PC3)


snps_fil_lfmm_candidate1 <- Subset(snps, sites=qv_adj.PC1)
snps_fil_lfmm_candidate1@site_id ##These are the candidate SNPs associated with PC1. 

snps_fil_lfmm_candidate2 <- Subset(snps, sites=qv_adj.PC2)
snps_fil_lfmm_candidate2@site_id ##These are the candidate SNPs associated with PC2.

snps_fil_lfmm_candidate3 <- Subset(snps, sites=qv_adj.PC3)
snps_fil_lfmm_candidate3@site_id ##These are the candidate SNPs associated with PC3.

## Remove duplicated snps names (some SNPs might be associated with more than on PC)
snps_candidate <- unique(c(snps_fil_lfmm_candidate1@site_id,snps_fil_lfmm_candidate2@site_id,
                           snps_fil_lfmm_candidate3@site_id)) 

length(snps_candidate)  ## 116
snps_fil_lfmm_candidate_all <- Subset(snps, sites=snps_candidate)
snps_fil_lfmm_candidate_all@site_id ##These are all the candidate SNPs. From here to FASTA files and BLAST

#### Save filtered vcf
Save(snps_fil_lfmm_candidate_all, "adapt_var_mapping/EAA_LFMM/snps_candidate_lfmm_ipomoea.vcf")


#------------------------------------------------------------------------------
#                     Save Candidates Loci and their Fasta
#------------------------------------------------------------------------------  
## Retrieve chromosome IDS
CANDIDATES_d1 <- Chrom(snps_fil_lfmm_candidate1)
CH1 <- unique(CANDIDATES_d1[, 1])

CANDIDATES_d2 <- Chrom(snps_fil_lfmm_candidate2)
CH2 <- unique(CANDIDATES_d2[, 1])

CANDIDATES_d3 <- Chrom(snps_fil_lfmm_candidate3)
CH3 <- unique(CANDIDATES_d3[, 1])

contigs <- unique(c(CH1,CH2,CH3))
length(contigs)

## Retrieve sequences from fasta file
fastafile <- read.fasta(file = "fasta/Icavalcantei.fasta", seqtype = "DNA", as.string = TRUE, set.attributes = FALSE)
head(fastafile)
length(fastafile)

SEQ1 <- fastafile[names(fastafile) %in% CH1]
SEQ2 <- fastafile[names(fastafile) %in% CH2]
SEQ3 <- fastafile[names(fastafile) %in% CH3]

total <- fastafile[names(fastafile) %in% contigs]

## Save new fasta files containing candidate sequences
write.fasta(sequences = SEQ1, names = names(SEQ1), nbchar = 150, file.out = "adapt_var_mapping/EAA_LFMM/SEQ1_lfmm_ipomoea.fasta")
write.fasta(sequences = SEQ2, names = names(SEQ2), nbchar = 150, file.out = "adapt_var_mapping/EAA_LFMM/SEQ2_lfmm_ipomoea.fasta")
write.fasta(sequences = SEQ3, names = names(SEQ3), nbchar = 150, file.out = "adapt_var_mapping/EAA_LFMM/SEQ3_lfmm_ipomoea.fasta")
write.fasta(sequences = total, names = names(total), nbchar = 150, file.out = "adapt_var_mapping/EAA_LFMM/SEQ_total_lfmm_ipomoea.fasta")

#END
