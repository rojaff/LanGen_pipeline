#This code serves to carry out sPCA analysis to build adaptation maps. 
#To run this code, you need first select adaptation loci. 
#Here we have 3 dataset: candidates loci selected by RDA, candidates loci selected by LFMM2, and combining RDA and LFMM2 candidate loci dataset.

rm(list=ls())

#Load packages

library("adegenet")
library("splancs")
library(r2vcftools)
library(usdm)

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

## RDA candidate snps
rda_snps <- vcfLink("adapt_var_mapping/GEA-RDA-Mahalanobis/snps_candidate_gea_RDA_mahalanobis_ipomoea.vcf", overwriteID=F) # here overwrite need to be F, otherwise we will loss the snps name
rda_snps@site_id
rda_snps_gen <- gen.imp[,rda_snps@site_id]
dim(rda_snps_gen)

#lfmm candidate snps
lfmm_snps <- vcfLink("adapt_var_mapping/EAA_LFMM/snps_candidate_lfmm_ipomoea.vcf", overwriteID=F)
lfmm_snps@site_id
lfmm_snps_gen <- gen.imp[,lfmm_snps@site_id]
dim(lfmm_snps_gen)

# lfmm and rda candidates
all_candidates <- unique(c(rda_snps@site_id,lfmm_snps@site_id))
all_candidate_gen <- gen.imp[,all_candidates] 
dim(all_candidate_gen)

### load coordinates

coord <- read.table("adapt_var_mapping/Dados_ambientais/climatic_data_ipomoea.txt", head=T)
espaco <- coord[,c(2,3)]
head(espaco)
espaco <- as.matrix(espaco)

#Here we will run sPCA following the tutorial provide by Jombart. I chose neighbourhood by distance to run sPCA

## Run spca

# all candidate loci

#1- run spca using the argument "scannf=F" and the Neighbourhood by distance method

mySpca <- spca(all_candidate_gen, xy=espaco, type=5, d1=0, d2="dmin", ask=F, scannf = F) ## type =5 used neighbourhood by distance. I chose 4 positive axis and 0 negatve axis
plot(mySpca)
summary(mySpca)
str(mySpca)

#2- test the local and global structure

myGtest <- global.rtest(espaco,mySpca$lw,nperm=999)
myGtest
plot(myGtest) ## significant

myLtest <- local.rtest(espaco,mySpca$lw,nperm=999)
myLtest
plot(myLtest) ## not significant

#3- run spca using the argument "scannf=T" and choose the number of axis - local and global

mySpca2 <- spca(all_candidate_gen,xy=espaco, ask=F, type=5, d1=0, d2="dmin", scannf=T) 
# Select the first number of axes (>=1): 3
# Select the second number of axes (>=0): 0
## I used 3 global axis and 0 local axis

plot(mySpca2)
summary(mySpca2)
head(mySpca2$li)

## plot of eigenvalues

barplot(mySpca2$eig, main="A variant of the plot\n of sPCA eigenvalues",
        col=spectral(length(mySpca2$eig)))
legend("topright", fill=spectral(2),
       leg=c("Global structures", "Local structures"))
abline(h=0,col="black")

# Colorplot of mySpca_candidate

colorplot(mySpca2$xy,as.data.frame(mySpca2$ls[,1:3]), cex=3, main="Colorplot of mySpca_candidate, three first global score")

# Extract to Interplolation

head(mySpca2$ls)
head(cbind(mySpca2$ls,espaco))
write.csv(cbind(mySpca2$ls,espaco), "adapt_var_mapping/sPCA_adaptation_maps/spca_axis_all_candidate_loci_ipomoea.csv")


## Make Adaptation Maps using QGIS      

#1. Open QGIS
#2. Go to "Add a delimited text file" and open the file with spca axes
#3. Save the opened file as a shapefile
#4. Interpolate each axis separately: Raster>Analysis>Grid(Interpolation)
#5. At the Grid(Interpolation) window, a) select the shapefile with spca axes, 
    # b)select "Z field" and choose one of the axis, c) select the "extension"
#6. Repeat the previous steps to each axis
#7. Combine the three interpolated axis raster in a rgb raster using Raster>Micelaneous>Mosaic. 
    # a) Open the three interpolated axis raster, b) select "put each file in a separatly band".
