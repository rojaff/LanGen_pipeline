#This code was made to perform fine-scale spatial genetic structure with local polynomial fitting (LOESS)


library(ggplot2)
library(r2vcftools)
library(geosphere)
library(raster)
library(rgdal)
library(gridExtra)
library("dplyr")

#######################################################
######## STEP 1: CALCULATE GENETIC DISTANCES ##########
#######################################################


######################################################################
#                                                                    #
#                      LOSS MANTEL            
#                                                                    #
######################################################################


# mat1 is a geographic distance matrix (a full matrix NOT a dist object)
# mat2 is a kinship matrix (a full matrix NOT a dist object)

mDistoLoess <- function(mat1,mat2,nperm=999){
  ltd <- mat1[lower.tri(mat1)]
  ltk <- mat2[lower.tri(mat2)]
  disto <- matrix(nrow = nperm, ncol = length(ltk))
  for(n in 1:nperm){
    perm <- sample(1:nrow(mat2))
    mat2_p <- mat2[perm,perm]
    lt <- mat2_p[lower.tri(mat2_p)]
    disto[n,] <- predict(loess(lt~ltd))
    if(n %% 100 == 0) print(n)
  }
  obs_fit <- loess(ltk~ltd)
  obs <- predict(obs_fit)
  CI <- apply(disto, 2, quantile, probs = c(0.025,0.975))
  attr(disto, "obs") <- obs
  attr(disto, "CI") <- CI
  return(list(disto=disto))
}



## Relatedness dataframe

snps_neutral <- vcfLink("vcf/ipomoea_filtered_ld_hw_neutral.vcf", overwriteID=T)

REL <- Relatedness(snps_neutral)
head(REL)

###Transform Relatedness dataframe into matrix

nb_SN <- length(unique(REL$INDV1))
nb_SN

rb_SN <- matrix(0, nb_SN, nb_SN)
rb_SN[lower.tri(rb_SN, diag=T)] <- REL$RELATEDNESS_AJK
rb_SN <- rb_SN + t(rb_SN) # to make symmetric
diag(rb_SN) <- diag(rb_SN)/2 # to correct diagonal for previous line
rb_SN

class(rb_SN)
dim(rb_SN)
nrow(rb_SN)
ncol(rb_SN)

###################################
# CALCULATE GEOGRAFPHIC DISTANCES #
###################################

coord_SN <- read.table("adapt_var_mapping/Dados_ambientais/climatic_data_ipomoea.txt", head=T)
head(coord_SN)
coordinates(coord_SN) <- coord_SN[,c(2,3)]
projection(coord_SN) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84")

sdist_SN  <- distm(coord_SN, fun=distGeo)
sdist_SN  <- sdist_SN /1000
nrow(sdist_SN)
ncol(sdist_SN)
max(sdist_SN)
dim(sdist_SN)
sdist_SN

ot_SN <- mDistoLoess(sdist_SN, rb_SN, 999)

dd_SN <- data.frame(cov = sdist_SN[lower.tri(sdist_SN)], obs = attr(ot_SN$disto, "obs"), lci = attr(ot_SN$disto, "CI")[1,], uci = attr(ot_SN$disto, "CI")[2,])

LP3_SN <- ggplot(dd_SN, aes(x=cov,y=obs)) + geom_line(size=1.2) + geom_ribbon(aes(ymax = uci, ymin = lci), alpha = 0.3) + theme_minimal() + 
  geom_hline(yintercept = mean(ot_SN$disto), lty = 3) + ylab("Mean Relatedness") + xlab("Distance (Km)") + geom_rug(sides = "b", alpha = 0.02) +
  annotate("text", x = 2.2, y = 0.10, label = "", size=10) +
  theme(axis.title.y = element_text(size=20, color = "black", face = "bold"),
        axis.title.x = element_text(size=20, color = "black", face = "bold")) + 
  theme(axis.text = element_text(face = "bold", color = "black", size = 12))

#tiff(filename = "Loss_Dioclea_SN.jpeg",
#     width = 20, height = 20, units = "cm",
 #    compression = "lzw", bg = "white", res = 300)
#LP3_SN
#dev.off()