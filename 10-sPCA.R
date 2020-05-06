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
#A. MAPPING ADAPTATIVE GENETIC VARIATION USING SPATIAL PRINCIPAL COMPONENTS ANALYSIS (sPCA)

##2. REMOVE ANY OBJECT OR FUNCTION IN THE ENVIRONMENT:
rm(list=ls())

##3. INSTALL AND LOAD THE PACKAGES
library(adegenet)
library(splancs)
library(r2vcftools)
library(usdm)

##4. INPUTS FOR THIS STEP:
#A. THE FILE ".VCF" FOR ADAPATATION ANALYSES CREATED IN FILTERING STEP 1, NAMED AS "_filtered_within_ld_test2.vcf".
#B. GENOTYPE FILE WITHOUT MISSING DATA CREATED IN FILTERING STEP 1, NAMED AS "_filtered_within_ld_imputed.geno".
#C. GEA RESULTS USING RDA AND LFMM ANALYSES FROM STEPS 6 AND 7.
#D. ENVIRONMENTAL INFORMATION FILTERED AND SELECTED IN STEP 3.

##5. LOAD FUNCTIONS TO BE USED ON THIS STEP.
VCFsummary <- function(snps){
  Nid <- nrow(snps@meta)
  Nsnps <- length(snps@site_id)
  cat(paste(Nid, "individuals and", Nsnps, "SNPs."), sep="\n")
}


#------------------------------------------------------------------------------
#                               Load Inputs 
#------------------------------------------------------------------------------
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

## RDA candidate snps
rda_snps <- vcfLink("adapt_var_mapping/GEA-RDA-Mahalanobis/snps_candidate_gea_RDA_mahalanobis_ipomoea.vcf", overwriteID=F) # here overwrite need to be F, otherwise we will loss the snps name
rda_snps@site_id
rda_snps_gen <- gen.imp[,rda_snps@site_id]
dim(rda_snps_gen)

## LFMM candidate snps
lfmm_snps <- vcfLink("adapt_var_mapping/EAA_LFMM/snps_candidate_lfmm_ipomoea.vcf", overwriteID=F)
lfmm_snps@site_id
lfmm_snps_gen <- gen.imp[,lfmm_snps@site_id]
dim(lfmm_snps_gen)

## LFMM and RDA candidates
all_candidates <- unique(c(rda_snps@site_id,lfmm_snps@site_id))
all_candidate_gen <- gen.imp[,all_candidates] 
dim(all_candidate_gen)

## Load coordinates
coord <- read.table("adapt_var_mapping/Dados_ambientais/climatic_data_ipomoea.txt", head=T)
espaco <- coord[,c(2,3)]
head(espaco)
espaco <- as.matrix(espaco)


#------------------------------------------------------------------------------
#                                Perform sPCA Analysis 
#------------------------------------------------------------------------------
#Here we will run sPCA following the tutorial provide by Jombart. We chose neighbourhood by distance to run sPCA

## All candidates
##1. Run sPCA using the argument "scannf=F" and the Neighbourhood by distance method
mySpca <- spca(all_candidate_gen, xy=espaco, type=5, d1=0, d2="dmin", ask=F, scannf = F) ## type =5 used neighbourhood by distance. I chose 4 positive axis and 0 negatve axis
plot(mySpca)
summary(mySpca)
str(mySpca)

##2. Test the local and global structure
myGtest <- global.rtest(espaco,mySpca$lw,nperm=999)
myGtest
plot(myGtest) ## significant

myLtest <- local.rtest(espaco,mySpca$lw,nperm=999)
myLtest
plot(myLtest) ## not significant

##3. Run sPCA using the argument "scannf=T" and choose the number of axis - local and global
mySpca2 <- spca(all_candidate_gen,xy=espaco, ask=F, type=5, d1=0, d2="dmin", scannf=T) 
# Select the first number of axes (>=1): 3
# Select the second number of axes (>=0): 0
## I used 3 global axis and 0 local axis

## Verify results
plot(mySpca2)
summary(mySpca2)
head(mySpca2$li)


#------------------------------------------------------------------------------
#                                Plot sPCA Results 
#------------------------------------------------------------------------------
## Plot of Eigenvalues
barplot(mySpca2$eig, main="A variant of the plot\n of sPCA eigenvalues",
        col=spectral(length(mySpca2$eig)))
legend("topright", fill=spectral(2),
       leg=c("Global structures", "Local structures"))
abline(h=0,col="black")

## Colorplot of mySpca_candidate
colorplot(mySpca2$xy,as.data.frame(mySpca2$ls[,1:3]), cex=3, main="Colorplot of mySpca_candidate, three first global score")

## Extract to Interplolation
head(mySpca2$ls)
head(cbind(mySpca2$ls,espaco))
write.csv(cbind(mySpca2$ls,espaco), "adapt_var_mapping/sPCA_adaptation_maps/spca_axis_all_candidate_loci_ipomoea.csv")


## Make Adaptation Maps using QGIS 3.4      
#1. Open QGIS and install the Processing plugin
#2. Load the file containing the three sPCA axes scores as "Delimited Text" and save it as a shapefile
#3. Interpolate each axis separately using the "IDW interpolation" tool from the Interpolation option in the Processing Toolbox
#4. Create an RGB composite raster combining the three created rasters, using Raster>Micelaneous>Merge and clicking the option "put each file in a separatly band"
#5. Customize maps as desired

## END
