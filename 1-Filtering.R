###############################################################################
####################### VALE INSTITUTE OF TECHNOLOGY ##########################
############### LABORATORIO DE GENETICA DA PAISAGEM - GENPAI ##################
###############################################################################


###############################################################################
##### COMBINING GENOTYPE, PHENOTYPE, AND ENVIRONMENTAL DATA TO DELINEATE ######
######  SITE-AJUSTED PROVENANCE STRATEGIES FOR ECOLOGICAL RESTORATION #########
###############################################################################
### AUTHORED BY: CAROLINA S. CARVALHO, BRENNA R. FORESTER, SIMONE K. MITRE, ###
########## RONNIE ALVES, VERA L. IMPERATRIZ-FONSECA, SILVIO J. RAMOS, #########
##### LUCIANA C. RESENDE-MOREIRA, JOSÉ 0. SIQUEIRA, LEONARDO C. TREVELIN, #####
############# CECILIO F. CALDEIRA, MARKUS GASTAUER, RODOLFO JAFFÉ #############
###############################################################################


###############################################################################
############################### PRE-ANALYSIS ##################################
###############################################################################


##1. GOALS FOR THIS STEP:
#A. FILTER RAW SNPS DATASET INTO NEUTRAL SNPS DATASET. 
#B. CREATE A DATASET THAT WILL BE USED IN RDA AND LFMM ANALYSES.


##2. REMOVE ANY OBJECT OR FUNCTION IN THE ENVIRONMENT:
rm(list=ls())


##3. INSTALL AND LOAD THE PACKAGES
library(r2vcftools)
library(LEA)


##4. LOAD FUNCTIONS TO BE USED ON THIS STEP.
VCFsummary <- function(snps){
  Nid <- nrow(snps@meta)
  Nsnps <- length(snps@site_id)
  cat(paste(Nid, "individuals and", Nsnps, "SNPs."), sep="\n")
}

Diference_Final <- function(All_SNPs, P1, P2, P3){ 
  Diference_1 <- intersect (All_SNPs@site_id, P1@site_id)
  Diference_2 <- intersect (All_SNPs@site_id, P2@site_id)
  Diference_3 <- intersect (All_SNPs@site_id, P3@site_id)
  
  Diference_All_1 <- intersect (Diference_1, Diference_2) 
  Diference_Final <- intersect (Diference_All_1, Diference_3)
  
  return(Diference_Final)
}


##5. DOWNLOAD A VCF FILE AS EXAMPLE "Icavalcantei.vcf" FROM FIGSHARE: https://doi.org/10.6084/m9.figshare.6100004.v1
#A. CREATE A FOLDER NAMED "vcf" IN YOUR WORKING DIRECTORY AND SAVE THE .vcf FILE THERE.


#######################################################################################
##################################### ANALYSES ########################################
#######################################################################################


#--------------------------------------------------------------------------------------
#    Load the VCF file and verify the quality of data and filter by biallelic snps 
#--------------------------------------------------------------------------------------
snps_raw <-  vcfLink("vcf/Icavalcantei.vcf", overwriteID=T)
snps_raw

## Remove loci that are not SNPs
snps <- Filter(snps_raw, filterOptions(max.alleles=2, min.alleles=2), indels="remove")

###Basic stats
VCFsummary(snps_raw) 
VCFsummary(snps)
snps@sample_id
snps@site_id
snps@meta
Chrom(snps) ##Chromosome, possitions, and IDs

##Remove repeated individual
snps@meta
UNIND <- snps@sample_id[snps@sample_id != "I096_rep_sorted.bam"]
snps_unind <- Subset(snps, samples=UNIND)

snps_unind@meta
VCFsummary(snps_unind) ## 122 individuals and 40229 SNPs.

# See GenotypeMatrix and remove loci and individuals with missing genotypes
genotypes <- GenotypeMatrix(snps_unind) # only returns biallelic
genotypes[1:10, 1:10] ## -1 is missing; otherwise, gives the number of derived alleles in the genotype -- ie. a 0 and 2 are both homozygotes

# Missing per locus
Missing <- apply(GenotypeMatrix(snps_unind), 2, function(x) sum(x < 0)/length(x)*100) ## Actual percentage of missing data
summary(Missing) ## Max missing = 30
hist(Missing)

### Look at depth, quality, HWE, HE, allele frequencies, and Pi
site.depth <- Query(snps_unind, type="site-mean-depth")
summary(site.depth$MEAN_DEPTH)
hist(site.depth$MEAN_DEPTH, breaks=30) 
hist(site.depth$MEAN_DEPTH[site.depth$MEAN_DEPTH <100]) ### >20 
hist(site.depth$MEAN_DEPTH[site.depth$MEAN_DEPTH <250]) ### <200

quality <- Query(snps_unind, type="site-quality")
summary(quality$QUAL)
hist(quality$QUAL) ## >30

PI <- Query(snps_unind, type="site-pi")
mean(PI$PI) ## Mean nucleotide divergency per-site
hist(PI$PI)

HWE <- Query(snps_unind, type="hardy")
summary(HWE$P_HWE)
hist(HWE$P_HWE)
hist(HWE$P_HWE[HWE$P_HWE<0.0001])

HE <- Query(snps_unind, type="het")
hist(HE$O.HOM) ## O.HOM, E.HOM, N_SITES, F
hist(HE$E.HOM) 
hist(HE$N_SITES) 
hist(HE$F) ## Inbreeding coefficient

freq <- Query(snps_unind, type="freq2")
hist(freq$X.FREQ.)
hist(freq$X.FREQ.1)
head(freq)

#--------------------------------------------------------------------------------------
#    FILTER ADAPTIVE loci by depth, quality, maf, snps with missing data and
#                 removing high LD loci within same contig
#--------------------------------------------------------------------------------------

snps_fil <- Filter(snps_unind, filterOptions(minQ=30, max.missing = 0.8, maf=0.05, min.meanDP=20, max.meanDP=200)) 
VCFsummary(snps_fil) 

### LD within contigs
#ld_within <- Linkage(snps_fil, type="geno-r2", linkageOptions(min.r2=0.0))
#ld_within <- read.csv("adapt_var_mapping/Filtering/ld_within_0.csv")
#head(ld_within)
#hist(ld_within$R.2)
#write.csv(ld_within, file="adapt_var_mapping/Filtering/ld_within_0.csv")

ld_within <- Linkage(snps_fil, type="geno-r2", linkageOptions(min.r2=0.8))
ld_within <- read.csv("adapt_var_mapping/Filtering/ld_within_0_8_test2.csv")
head(ld_within)
hist(ld_within$R.2)
#write.csv(ld_within, file="adapt_var_mapping/Filtering/ld_within_0_8_test2.csv")


#### Select one set of the correlated snps (ID1 or ID2)

ld_snps <- ld_within$ID1
nold_snps <- snps_fil@site_id[!(snps_fil@site_id %in% ld_snps)] 
snps_fil_ld <- Subset(snps_fil, sites=nold_snps) # Keep snps that are not in LD.
VCFsummary(snps_fil_ld) ## 122 individuals and 17025 SNPs.


# Missing per individual

Missing_ind <- apply(GenotypeMatrix(snps_fil_ld),1, function(x) sum(x<0)/length(x)*100)
summary(Missing_ind) ## Max missing = 99
hist(Missing_ind)
Missingind.df <- as.data.frame(Missing_ind)
Missingind.df$ID <- row.names(Missingind.df)
Missingind.df[Missingind.df$Missing_ind>70,] # define threshold for individual missing data tolerance - shows individuals that should be removed from ID data

## Remove individuals with large amounts of missing (set to max 50% missing)

indtokeep <- Missingind.df[Missingind.df$Missing_ind <= 70,]
snps_lowindmiss <- Subset(snps_fil_ld, samples = indtokeep$ID)
snps_lowindmiss@meta
VCFsummary(snps_unind)
VCFsummary(snps_lowindmiss) ## 115 individuals and 17025 SNPs.


#### Save filtered vcf
Save(snps_lowindmiss, "vcf/ipomoea_filtered_within_ld_test2.vcf")

#--------------------------------------------------------------------------------------
#           FILTER NEUTRAL loci by depth, quality, maf, missing data and hwe
#--------------------------------------------------------------------------------------

snps_fil_hwe <- Filter(snps_unind, filterOptions(minQ=30, max.missing = 0.8, maf=0.05, min.meanDP=20, max.meanDP=200, hwe=0.0001)) 
VCFsummary(snps_fil_hwe)  ## 122 individuals and 15077 SNPs.

## If you have more than one population, use the code below. This code identifies SNPs deviating from HW equilibrium
## within each population, and then removes those SNPs that are in desequilibrium in all populations.

### LD within contigs
#ld_within <- Linkage(snps_fil_hwe, type="geno-r2", linkageOptions(min.r2=0.0))
#ld_within <- read.csv("adapt_var_mapping/Filtering/ld_within_0_hwe_test2.csv")
#head(ld_within)
#hist(ld_within$R.2)
#write.csv(ld_within, file="adapt_var_mapping/Filtering/ld_within_0_hwe_test2.csv")

#ld_within <- Linkage(snps_fil_hwe, type="geno-r2", linkageOptions(min.r2=0.8))
ld_within <- read.csv("adapt_var_mapping/Filtering/ld_within_0_8_hwe_test2.csv")
head(ld_within)
hist(ld_within$R.2)
#write.csv(ld_within, file="adapt_var_mapping/Filtering/ld_within_0_8_hwe_test2.csv")

#### Select one set of the correlated snps (ID1 or ID2)
ld_snps <- ld_within$ID1
nold_snps <- snps_fil_hwe@site_id[!(snps_fil_hwe@site_id %in% ld_snps)] 
snps_fil_ld <- Subset(snps_fil_hwe, sites=nold_snps) # Keep snps that are not in LD.
VCFsummary(snps_fil_ld) ## 13604 SNPs

## LD between contigs
#ld_between <- Linkage(snps_fil_ld, type="interchrom-geno-r2", linkageOptions(min.r2=0)) 
#ld_between <- read.csv("adapt_var_mapping/Filtering/ld_between_0_hwe_test2.csv")
#head(ld_between)
#hist(ld_between$R.2)
#write.csv(ld_between, file="adapt_var_mapping/Filtering/ld_between_0_hwe_test2.csv")

#ld_between <- Linkage(snps_fil_ld, type="interchrom-geno-r2", linkageOptions(min.r2=0.8)) 
ld_between <- read.csv("adapt_var_mapping/Filtering/ld_between_0_8_hwe_test2.csv")
head(ld_between)
hist(ld_between$R.2)
#write.csv(ld_between, file="adapt_var_mapping/Filtering/ld_between_0_8_hwe_test2.csv")

ld2_snps <- ld_between$ID1
nold2_snps <- snps_fil_ld@site_id[!(snps_fil_ld@site_id %in% ld2_snps)]
snps_fil_ldF <- Subset(snps_fil_ld, sites=nold2_snps) # Keep snps that are not in LD.
VCFsummary(snps_fil_ldF) ## 122 individuals and 13425 SNPs.

######## Check if filtering worked fine
site.depth2 <- Query(snps_fil_ldF, "site-mean-depth")
quality2 <- Query(snps_fil_ldF, "site-quality")
HWE2 <- Query(snps_fil_ldF, type="hardy")

hist(site.depth2$MEAN_DEPTH, breaks=30)
hist(site.depth2$MEAN_DEPTH[site.depth2$MEAN_DEPTH <200])
hist(quality2$QUAL)
hist(HWE2$P_HWE[HWE2$P_HWE<200])


Missing_ind <- apply(GenotypeMatrix(snps_fil_ldF),1, function(x) sum(x<0)/length(x)*100)
summary(Missing_ind) ## Max missing = 99
hist(Missing_ind)
Missingind.df <- as.data.frame(Missing_ind)
Missingind.df$ID <- row.names(Missingind.df)
Missingind.df[Missingind.df$Missing_ind>70,] # define threshold for individual missing data tolerance - shows individuals that should be removed from ID data

## Remove individuals with large amounts of missing (set to max 50% missing)

indtokeep <- Missingind.df[Missingind.df$Missing_ind <= 70,]
snps_lowindmiss <- Subset(snps_fil_ldF, samples = indtokeep$ID)
snps_lowindmiss@meta
VCFsummary(snps_unind)
VCFsummary(snps_lowindmiss)


#### Save filtered vcf
Save(snps_lowindmiss, "vcf/ipomoea_filtered_ld_hwe_test2.vcf")

#--------------------------------------------------------------------------------------
#                       Dataset for RDA and LFMM analyses
#--------------------------------------------------------------------------------------
#Some analysis as RDA and LFMM2 does not allow missing data, so we used LFMM to
#impute genetic missing data, based on the population that each individual
#belongs, using the package LEA

snps <-  vcfLink("vcf/ipomoea_filtered_within_ld_test2.vcf", overwriteID=T)
snps@meta
VCFsummary(snps)

gen <- GenotypeMatrix(snps)
gen[gen==-1] <- NA
sum(is.na(gen))
dim(gen)

###Create LFMM file

snps_fil_lfmm <- Lfmm(snps, "vcf/ipomoea_filtered_within_ld.lfmm")
project.snmf = snmf("vcf/ipomoea_filtered_within_ld.lfmm", K = 1, 
                    entropy = TRUE, repetitions = 10,
                    project = "new")

# select the run with the lowest cross-entropy value
best = which.min(LEA::cross.entropy(project.snmf, K = 1))

# Impute the missing genotypes
impute(project.snmf, "vcf/ipomoea_filtered_within_ld.lfmm", method = 'mode', K = 1, run = best)

# Convert lfmm to geno and save
lfmm2geno("vcf/ipomoea_filtered_within_ld.lfmm_imputed.lfmm", output.file = "vcf/ipomoea_filtered_within_ld_imputed.geno", force = TRUE)

#END
