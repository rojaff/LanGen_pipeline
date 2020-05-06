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
#A. PREDICT ADAPTATIVE GENOTYPES FROM RDA MODELS AND ENVIRONMENTAL DATA

##2. REMOVE ANY OBJECT OR FUNCTION IN THE ENVIRONMENT:
rm(list=ls())

##3. INSTALL AND LOAD THE PACKAGES
library(r2vcftools)
library(vegan)    
library(usdm)
library(factoextra)
library(NbClust)
library(ggplot2)
library(dplyr)
library(LEA)

##4. INPUTS FOR THIS STEP:
#A. THE FILE ".VCF" FOR ADAPATATION ANALYSES CREATED IN FILTERING STEP 1, NAMED AS "_filtered_within_ld_test2.vcf".
#B. GENOTYPE FILE WITHOUT MISSING DATA CREATED IN FILTERING STEP 1, NAMED AS "_filtered_within_ld_imputed.geno".
#C. GEA RESULTS USING RDA AND LFMM ANALYSES FROM STEPS 6 AND 7.
#D. ENVIRONMENTAL INFORMATION FILTERED AND SELECTED IN STEP 3.
#E. BIOCLIMATIC VARIABLES FROM WORLDCLIM.
#F. INFORMATION TO PERFORM THE PREDICTIONS IN THIS CASE BASED ON MINE INFORMATION.

##5. LOAD FUNCTIONS TO BE USED ON THIS STEP.
VCFsummary <- function(snps){
  Nid <- nrow(snps@meta)
  Nsnps <- length(snps@site_id)
  cat(paste(Nid, "individuals and", Nsnps, "SNPs."), sep="\n")
}


#------------------------------------------------------------------------------
#                               Load GEA Results 
#------------------------------------------------------------------------------
## This code serves to predict which plants are most adapted to climatic conditions using the function predict.cca with RDA.

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

## Load rda and lfmm candidate snps - genotype vs enviroment
rda_snps <- vcfLink("adapt_var_mapping/GEA-RDA-Mahalanobis/snps_candidate_gea_RDA_mahalanobis_ipomoea.vcf", overwriteID=F) # here overwrite need to be F, otherwise we will loss the snps name
lfmm_snps <- vcfLink("adapt_var_mapping/EAA_LFMM/snps_candidate_lfmm_ipomoea.vcf", overwriteID=F)

## GEA candidate snps
gea <- unique(c(rda_snps@site_id,lfmm_snps@site_id))
gea_gen <- gen.imp[,gea]
dim(gea_gen)


#------------------------------------------------------------------------------
#                          Load Environmental Dataset 
#------------------------------------------------------------------------------
## Load environment data set
dados_env <- read.table("adapt_var_mapping/Dados_ambientais/climatic_data_ipomoea.txt", head=T)
head(dados_env)
str(dados_env)
dados_env$Descricao2 <- rep("SN", nrow(dados_env))

## Load mining data
dados_mine <- read.table("adapt_var_mapping/RDA_predict/data_mining3", head=T)
head(dados_mine)
str(dados_mine)

## Select the coordinates
I <-  dados_mine
head(I)
str(I)
coordinates(I) <- I[,c(3,2)]
plot(I)

## Build raster stack with all climatic layers
dir <- "adapt_var_mapping/Dados_ambientais/wc0.5/"
current.list <- list.files(path=dir, pattern =".bil", full.names=TRUE)
Bioclimate <- stack(current.list)
names(Bioclimate)
projection(Bioclimate) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +units=m"

## Extract climatic data
mining_climate <- raster::extract(Bioclimate, I)
head(mining_climate)

## Create dataframe with climatic and soil data of the mining points 
dados_mine$bio06 <- mining_climate[,c("bio6_34")]
dados_mine$bio18 <- mining_climate[,c("bio18_34")]
dados_mine$bio16 <- mining_climate[,c("bio16_34")]

## Scaling dataset
ID <- c(dados_env$ID,dados_mine$Ponto)
Descricao <- c(as.character(dados_env$Descricao2), as.character(dados_mine$Descricao2))
Long <- c(dados_env$longDD, dados_mine$Longitude)
Lat <- c(dados_env$latDD, dados_mine$Latitude)
bio06 <- scale(c(dados_env$bio06, dados_mine$bio06))
bio18 <- scale(c(dados_env$bio18, dados_mine$bio18))
bio16 <- scale(c(dados_env$bio16, dados_mine$bio16))

dados_scaled <- data.frame(ID, Descricao, Long, Lat, bio06, bio18, bio16)
head(dados_scaled)

dados_env <- subset(dados_scaled, Descricao=="SN")
str(dados_env)
dados_mine <- subset(dados_scaled, Descricao!="SN")
str(dados_mine)

#------------------------------------------------------------------------------
#                               Perform RDA Analysis 
#------------------------------------------------------------------------------
## Now we will run a RDA using the gea candidates loci,
#create a new data set and use the function predict to predict values from rda models.
# It is not possible to use predict.cca with condition argument. So we will run rda without take account with population structure

## Rub RDA analysis
gea_rda <- rda(gea_gen ~ bio06+bio18+bio16, data=as.data.frame(dados_env[,-c(1:4)]))
summary(gea_rda)        
RsquareAdj(gea_rda)

## Our constrained ordination explains about 16% of the variation
summary(eigenvals(gea_rda, model = "constrained")) #The eigenvalues for the constrained axes reflect the variance explained by each canonical axis
screeplot(gea_rda) 
signif.full <- anova.cca(gea_rda, by="axis", parallel=getOption("mc.cores"), permutations = 999) #  check our RDA model for significance using formal tests.

## Predict data to new data set using the function predict
predict_rda <- predict(gea_rda, newdata=as.data.frame(dados_mine[,-c(1:4)]), type = "lc",scal=2)
predict_rda

## data.frame with predicted scores for minepit and observed score for Ipomoea individuals
scores_rda <- rbind(predict_rda[,1:2],scores(gea_rda, choices = 1:2)$site)
areas <- c(as.character(dados_mine$Descricao),as.character(dados_env$Descricao))
row.names(scores_rda)<-paste(areas,1:nrow(scores_rda))

## Kmeans

## Evaluate initial number of clusters: indices for choosing the best number of clusters
NbClust(data = scores_rda, distance = "euclidean",
        min.nc = 2, max.nc = 10, method ="kmeans")

#According to the majority rule, the best number of clusters is  3

## Function for k-means clustering ####
km.res <- kmeans(scores_rda, centers =  3, nstart = 25)

# centers = number of clusters
# nstarts = number of random centroids to start procedure

print(km.res)


#------------------------------------------------------------------------------
#                                  Plot Results 
#------------------------------------------------------------------------------
## Visualizing clustering
p <- fviz_cluster(km.res, scores_rda, ellipse.type = "norm", geom="point", stand = F)
data <- p$data
data$area <- areas

ggplot(data, aes(x, y)) + geom_point(size=4,aes(color=cluster,shape=area)) +
  stat_ellipse(aes(x=x, y=y,color=cluster),type="norm")+
  theme_minimal()+ theme(axis.title.x = element_text(size=20, color = "black", face = "bold"),
                         axis.title.y = element_text(size=20, color = "black", face = "bold"),
                         axis.text.y = element_text(size=20, color = "black", face = "bold"),
                         axis.text.x = element_text(size=20, color = "black", face = "bold"),
                         legend.position="none")+
  theme(legend.position = "right")

## Correlation axes
predict_rda2 <- predict(gea_rda, newdata=as.data.frame(pred_scale), type = "lc",scal=2)
predict_rda2

scores(gea_rda, choices = 1)$site
  
plot(predict_rda2[,1],scores(gea_rda, choices = 1)$site)
plot(predict_rda2[,2],scores(gea_rda, choices = 2)$site)
plot(predict_rda2[,3],scores(gea_rda, choices = 3)$site)

##END
