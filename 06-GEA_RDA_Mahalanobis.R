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


#------------------------------------------------------------------------------
#                               PRE-ANALYSIS 
#------------------------------------------------------------------------------

##1. GOALS FOR THIS STEP:
#A. PERFORM RDA ANALYSIS TO IDENTIFY LOCI UNDER SELECTION (OUTLIER DETECTION)

##2. REMOVE ANY OBJECT OR FUNCTION IN THE ENVIRONMENT:
rm(list=ls())

##3. INSTALL AND LOAD THE PACKAGES
library(vegan)
library(usdm)
library(ade4)
library(psych)
library(seqinr)
library(robust)  
library(qvalue)
library(r2vcftools)
library(LEA)

##4. INPUTS FOR THIS STEP:
#A. THE FILE ".VCF" CLEANED AFTER FILTERING STEP 1.
#B. DOWNLOAD A VCF FILE AS EXAMPLE "Icavalcantei.vcf" FROM FIGSHARE: https://doi.org/10.6084/m9.figshare.6100004.v1
#C. CREATE A FOLDER NAMED "vcf" IN YOUR WORKING DIRECTORY AND SAVE THE .vcf FILE THERE.
#D. ".CSV" FILE WITH GEOGRAPHICAL INFORMATION FOR THE GENETIC SAMPLES. YOU CAN DOWNLOAD IT THIS EXAMPLE IN https://www.frontiersin.org/articles/10.3389/fpls.2018.00532/full#supplementary-material

##5. LOAD FUNCTIONS TO BE USED ON THIS STEP.
VCFsummary <- function(snps){
  Nid <- nrow(snps@meta)
  Nsnps <- length(snps@site_id)
  cat(paste(Nid, "individuals and", Nsnps, "SNPs."), sep="\n")
}


#------------------------------------------------------------------------------
#                            Load SNPs Dataset 
#------------------------------------------------------------------------------
## This script was developed by Brenna Forester and adapted to our data

## Load snps data set to take the snps and individuals ID
snps <-  vcfLink("vcf/ipomoea_filtered_within_ld_test2.vcf", overwriteID=T)
snps@meta
VCFsummary(snps) ## 115 individuals and 17025 SNPs.

## Read geno file
gen.imp <- read.geno("vcf/ipomoea_filtered_within_ld_imputed.geno")
gen.imp[1:10,1:10]
colnames(gen.imp) <- snps@site_id
rownames(gen.imp) <- snps@sample_id
dim(gen.imp)

## Number of missing data
sum(is.na(gen.imp))


#------------------------------------------------------------------------------
#                          Load Environmental Dataset 
#------------------------------------------------------------------------------
## Load environmental data set
dados_env <- read.table("adapt_var_mapping/Dados_ambientais/climatic_data_ipomoea.txt", head=T)
head(dados_env)
str(dados_env)

## Confirm that genotypes and environmental data are in the same order
data.frame(rownames(gen.imp),as.character(dados_env$ID))

## Check the correlation among enviromental variables (VIF < 3)
vif(dados_env[,-c(1:3)])

## Scale environmental dataset
pred_scale <- scale(dados_env[,-c(1:3)]) ## remove information of lat, long, location, etc
rownames(pred_scale) <-dados_env[,1]
head(pred_scale)


#------------------------------------------------------------------------------
#                             Perform RDA Analysis 
#------------------------------------------------------------------------------
#If you have more than one population, use the population assigned for each individual 
#in the Condition argument of RDA to control for population structure
# m.rda <- rda(gen.imp ~ pred_scale + Condition(as.factor(pop_lea))) 

## RDA with mahalanobis distance
m.rda <- rda(gen.imp ~ pred_scale) 
m.rda

RsquareAdj(m.rda)
screeplot(m.rda)
plot(m.rda, scaling=3)

summary(eigenvals(m.rda, model= "constrained"))

## Check axis significance
#signif.axis <- anova(m.rda, by= "axis", nperm= 99) 
## we chose the first 2 axes, which explained more than 80% of the variation

## Save loading of the first and second pcs
load.rda <- summary(m.rda)$species[,1:2]
K <- 2  # the number of RDA axes you're looking at


#------------------------------------------------------------------------------
#                             Outlier Detection 
#------------------------------------------------------------------------------
## Estimate adjusted p-values based on gif
zscale <- apply(load.rda, 2, scale)        # I'm not sure scaling is the best approach here...
mscale <- covRob(zscale, distance=TRUE, na.action=na.omit, estim="pairwiseGK")$dist
(gif <- median(mscale)/qchisq(0.5, df=K)) #1.45

## Look at histogram
rda.pval <- pchisq(mscale/gif, df=K, lower.tail=FALSE)   # Remember: you can always change the GIF!!
hist(rda.pval)

## Try different gif values if necessary
#rda.pval2 <- pchisq(mscale/1.1, df=K, lower.tail=FALSE)  # Less conservative GIF (smaller)
#hist(rda.pval2)

## Select snps with FDR < 0.05
rda.qval <- qvalue(rda.pval)$qvalues
m.snps <- (rda.FDR <- colnames(gen.imp)[which(rda.qval < 0.05)])

# Let's look at where these candidates are in the ordination space:
snp.color <- as.data.frame(colnames(gen.imp), stringsAsFactors=F)
snp.color[,2] <- apply(snp.color, 1, function(x) if(x %in% rda.FDR) 'gray32' else '#00000000')
snp.color[,3] <- apply(snp.color, 1, function(x) if(x %in% rda.FDR) '#e31a1c' else '#00000000')

## Axes 1 & 2
plot(m.rda, type="n", scaling=3, xlim=c(-1,1), ylim=c(-1,1))
points(m.rda, display="species", pch=21, cex=1, col="gray32", bg='#f1eef6', scaling=3)
points(m.rda, display="species", pch=21, cex=1, col=snp.color[,2], bg=snp.color[,3], scaling=3)
text(m.rda, scaling=3, display="bp", col="#0868ac", cex=1)

## Axes 2 & 3
plot(m.rda, type="n", scaling=3, xlim=c(-1,1), ylim=c(-1,1), choices=c(2,3))
points(m.rda, display="species", pch=21, cex=1, col="gray32", bg='#f1eef6', scaling=3, choices=c(2,3))
points(m.rda, display="species", pch=21, cex=1, col=snp.color[,2], bg=snp.color[,3], scaling=3, choices=c(2,3))
text(m.rda, scaling=3, display="bp", col="#0868ac", cex=1, choices=c(2,3))

## Save the snps names to compare with other approaches
write.table(m.snps, "adapt_var_mapping/GEA-RDA-Mahalanobis/ipomoea_snps_mahalanobis_gea.txt")

                       
#------------------------------------------------------------------------------
#        Correlation among Candidate SNPs and Environmental Predictors 
#------------------------------------------------------------------------------                                   
## Add in the correlations of each candidate SNP with the environmental predictors:
foo <- matrix(nrow=length(m.snps), ncol=3)  # 3 columns for 3 predictors
colnames(foo) <- c("bio06","bio18", "bio16")

for (i in 1:length(m.snps)) {
  nam <- m.snps[i]
  snp.gen <- gen.imp[,nam]
  foo[i,] <- apply(pred_scale,2,function(x) cor(x,snp.gen, method="spearman"))
}

cand <- cbind.data.frame(m.snps,foo)  
head(cand)

## Looking for duplicate detections
length(cand$m.snps[duplicated(cand$m.snps)])
cand <- cand[!duplicated(cand$m.snps),] # remove duplicate detections

## See which of the predictors each candidate SNP is most strongly correlated with
for (i in 1:length(cand$m.snps)) {
  bar <- cand[i,]
  cand[i,5] <- names(which.max(abs(bar[2:4]))) # gives the variable
  cand[i,6] <- max(abs(bar[2:4]))              # gives the absolute correlation value
  cand[i,7] <- ifelse (abs(max((bar[2:4])))>abs(min((bar[2:4]))),max((bar[2:4])), min((bar[2:4]))) # gives the raw correlation value
  
}

colnames(cand)[5] <- "predictor"
colnames(cand)[6] <- "correlation"
colnames(cand)[7] <- "correlation_real"

table(cand$predictor) 
head(cand)

## Plot the SNPs
sel <- cand$m.snps
env <- cand$predictor

env[env=="bio06"] <- '#1f78b4'
env[env=="bio18"] <- '#e31a1c'
env[env=="bio16"] <- '#b2df8a'


## Color by predictor:
col.pred <- rownames(m.rda$CCA$v) # pull the SNP names

for (i in 1:length(sel)) {           # color code candidate SNPs
  foo <- match(sel[i],col.pred)
  col.pred[foo] <- env[i]
}

col.pred[grep("snp",col.pred)] <- '#f1eef6' # non-candidate SNPs
empty <- col.pred
empty[grep("#f1eef6",empty)] <- rgb(0,1,0, alpha=0) # transparent
empty.outline <- ifelse(empty=="#00FF0000","#00FF0000","gray32")
bg <- c('#1f78b4','#e31a1c','#b2df8a')

# Axes 1 & 2
plot(m.rda, type="n", scaling=3, ylim=c(-1,1), xlim=c(-1,1))
points(m.rda, display="species", pch=21, cex=1, col="gray32", bg=col.pred, scaling=3)
points(m.rda, display="species", pch=21, cex=1, col=empty.outline, bg=empty, scaling=3)
text(m.rda, scaling=3, display="bp", col="#0868ac", cex=1)
legend("bottomright", legend=c("bio06","bio18","bio16"), bty="n", col="gray32", pch=21, cex=1, pt.bg=bg)


#------------------------------------------------------------------------------
#                     Save Candidates Loci and their Fasta
#------------------------------------------------------------------------------                  
## Load snps data set 
snps <-  vcfLink("vcf/ipomoea_filtered_within_ld_test2.vcf", overwriteID=T)
snps@meta
VCFsummary(snps) ## 115 individuals and 17025 SNPs.

sel ##These are all the candidate SNPs. From here to FASTA files and BLAST

snps_fil_rda_candidate <- Subset(snps, sites=as.character(sel))
snps_fil_rda_candidate@site_id ##These are all the candidate SNPs. From here to FASTA files and BLAST

## Save filtered vcf
Save(snps_fil_rda_candidate, "adapt_var_mapping/GEA-RDA-Mahalanobis/snps_candidate_gea_RDA_mahalanobis_ipomoea.vcf")

## Retrieve chromosome IDS
snp_chrom <- Chrom(snps_fil_rda_candidate)
CH1 <- unique(snp_chrom[, 1]) #225

## Retrieve sequences from fasta file 
fastafile <- read.fasta(file = "fasta/Icavalcantei.fasta", seqtype = "DNA", as.string = TRUE, set.attributes = FALSE)
head(fastafile)
length(fastafile)

SEQ1 <- fastafile[names(fastafile) %in% CH1]
length(SEQ1)

## Save new fasta files containing candidate sequences

seqinr::write.fasta(sequences = SEQ1, names = names(SEQ1), nbchar = 150, file.out = "adapt_var_mapping/GEA-RDA-Mahalanobis/SEQ_gea_RDA_mahalanobis_ipomoea.fasta")

##END
