# This code serves to carry out PCAdapt to outlier detection. This method is not a GEA
# This script was first developed by Brenna Forester and adapted to our data.

# pcadapt is a differentiation-based outlier method (no environmental data required)
# It can be run on individuals or pool-seq data.
# It does not require the identification of discrete populations, and performs well
# even in cases of hierarchical population structure and individual admixture.

# Like almost all differentiation-based methods, it is a univariate test which means
# we need to correct for multiple testing...more on this below.

rm(list=ls())

## Load packages

library(pcadapt)    # outlier detection & PCA
library(qvalue)     # FDR calculation
library(psych)      # nice pair plots
library(seqinr)

##Load functions that will be used throughout the code

### Summary stats
VCFsummary <- function(snps){
  Nid <- nrow(snps@meta)
  Nsnps <- length(snps@site_id)
  cat(paste(Nid, "individuals and", Nsnps, "SNPs."), sep="\n")
}

##First of all, load imputed SNP and the environmetal datasets

### Load snps data set to take the snps and individuals ID

snps <-  vcfLink("vcf/ipomoea_filtered_within_ld_test2.vcf", overwriteID=T)
snps@meta
VCFsummary(snps)

## read geno file

gen <- read.geno("vcf/ipomoea_filtered_within_ld_imputed.geno")
gen[1:10,1:10]
colnames(gen) <- snps@site_id
rownames(gen) <- snps@sample_id
dim(gen)

# Convert our genetic data to a pcadapt matrix
# use type="lfmm" when individuals are in rows and snps are in columns
gen.pcadapt <- read.pcadapt(gen, type="lfmm")

# check on what we're doing here:
?pcadapt

# PCA to determine K
# ------------------
ipomoea.pca <- pcadapt(gen.pcadapt, K=20, ploidy=2, min.maf=0.05)

# PCA screeplot
plot(ipomoea.pca, option="screeplot")
(ipomoea.pca$singular.values)^2       # proportion of variance explained by each PC
# The authors recommend "Cattell's Rule" for interpreting the scree plot (PC to the left of the flat line)

# Run pcadapt
# -----------
ipomoea.pcadapt <- pcadapt(gen.pcadapt, K=3, min.maf=0.05, 
                        method="mahalanobis", ploidy=2)


# Use Mahalanobis distance to compute the (rescaled) test statistics (z-scores in this case).
# The robust Mahalanobis distance is a metric that identifies outliers in multidimensional space.
# "Robust" means the estimation is not sensitive to outliers in the covariance matrix of the z-scores.

names(ipomoea.pcadapt)

# The Genomic Inflation Factor is stored in "gif" (more on that below)
# The (rescaled = stat/gif) test statistics are stored in "chi2.stat"
# The recalibrated p-values are stored in "pvalues"

plot(ipomoea.pcadapt, option="manhattan")
plot(ipomoea.pcadapt, option="qqplot")

# The GIF indicates how well the test is "calibrated".
# It corrects for inflation of the test score at each locus, which can occur when population
# structure or other confounding factors are not appropriately accounted for in the model.

ipomoea.pcadapt$gif

# GIF of 1=well calibrated, >1=liberal (too many small p-values), <1=conservative (too few small p-values)
# Note: GIFs > 2 indicate poor test calibration

hist(ipomoea.pcadapt$pvalues, xlab="p-values", main=NULL, breaks=50)
# Note: these p-values are already rescaled with the GIF...


# Q-values are calculated from p-values. They measure the expected proportion of false positives
# among the list of positive tests, a.k.a the "False Discovery Rate".

qv_adj<- which(qvalue::qvalue(ipomoea.pcadapt$pvalues, fdr=0.05)$signif)
length(qv_adj)

## save candidates loci

snps_pcadapt_candidate <- Subset(snps, sites=qv_adj)
snps_pcadapt_candidate@site_id ##These are all the candidate SNPs. From here to FASTA files and BLAST

#### Save filtered vcf
Save(snps_pcadapt_candidate, "adapt_var_mapping/PCAdapt/snps_candidate_pcadapt_ipomoea.vcf")

#### Retrieve chromosome IDS
CANDIDATES_d1 <- Chrom(snps_pcadapt_candidate)
CH1 <- unique(CANDIDATES_d1[, 1])

## Retrieve sequences from fasta file
fastafile <- read.fasta(file = "fasta/Icavalcantei.fasta", seqtype = "DNA", as.string = TRUE, set.attributes = FALSE)
head(fastafile)
length(fastafile)

SEQ1 <- fastafile[names(fastafile) %in% CH1]

### Save new fasta files containing candidate sequences

seqinr::write.fasta(sequences = SEQ1, names = names(SEQ1), nbchar = 150, file.out = "adapt_var_mapping/PCAdapt/SEQ1_pcadapt_ipomoea.fasta")
