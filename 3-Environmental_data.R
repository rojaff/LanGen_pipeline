#This code was made to evaluate the environmental variables
#that we will be used to select loci under selection for Ipomoea cavalcantei

# Loading packages

require(missMDA)
require(ggfortify)
library(raster)
library(rgdal)
library(sp)
library(rgeos)


#################################################################################
### Load coords of samples with genotype data

## Download coordinate data from the repository: coords_genotype_icav_UTM.txt
## Remove samples with missing data

coords <- read.table("data_prep/coords_genotype_icav_UTM", sep="\t", h=T)
head(coords)
str(coords)

coords_gen <- coords[coords$ID!="I13" & coords$ID!="I41" & coords$ID!="I48" & 
                       coords$ID!="I52" & coords$ID!="I101" & coords$ID!="I110" & coords$ID!="I116",]

str(coords_gen) ## 115 obs

## Transform into Spatial Points Dataframe
coordinates(coords_gen) <- coords_gen[ ,2:3]
nrow(coords_gen) #122 individuals
class(coords_gen) #"SpatialPointsDataFrame"
projection(coords_gen) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +units=m")

plot(coords_gen, axes=TRUE)


#Now, we selected the climatic variables from WorlClim dataset. 
#For that, we carried out a PCA to select the variables that are more correlated 
#to the first three axes. Then we built a climatic dataframe that will be used 
#in the subsequent analysis.


### WorldClim data
Mlong <- mean(coords_gen@coords[,1])
Mlat <- mean(coords_gen@coords[,2])
Bioclimate = getData('worldclim', var='bio', res=0.5, path = "adapt_var_mapping/Dados_ambientais", lon=Mlong, lat=Mlat) # Para baixar as camadas com alta resolução precisa definir as coordenadas: getData('worldclim', var='bio', res=0.5, lon=5, lat=45)

#Alt = Altitude
#BIO1 = Annual Mean Temperature
#BIO2 = Mean Diurnal Range (Mean of monthly (max temp - min temp))
#BIO3 = Isothermality (BIO2/BIO7) (* 100)
#BIO4 = Temperature Seasonality (standard deviation *100)
#BIO5 = Max Temperature of Warmest Month
#BIO6 = Min Temperature of Coldest Month
#BIO7 = Temperature Annual Range (BIO5-BIO6)
#BIO8 = Mean Temperature of Wettest Quarter
#BIO9 = Mean Temperature of Driest Quarter
#BIO10 = Mean Temperature of Warmest Quarter
#BIO11 = Mean Temperature of Coldest Quarter
#BIO12 = Annual Precipitation
#BIO13 = Precipitation of Wettest Month
#BIO14 = Precipitation of Driest Month
#BIO15 = Precipitation Seasonality (Coefficient of Variation)
#BIO16 = Precipitation of Wettest Quarter
#BIO17 = Precipitation of Driest Quarter
#BIO18 = Precipitation of Warmest Quarter
#BIO19 = Precipitation of Coldest Quarter


#### Build raster stack with all layers
getwd()
dir <- "adapt_var_mapping/Dados_ambientais/wc0.5/" # indicate the path to bioclimatic data
current.list <- list.files(path=dir, pattern =".bil", full.names=TRUE)
Bioclimate <- stack(current.list)
names(Bioclimate)
projection(Bioclimate) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +units=m"


##### Extract data
Idata <- extract(Bioclimate, coords_gen)
Idata <- as.data.frame(Idata)
class(Idata)
head(Idata)
summary(Idata)

######### Run PCA and extract the first components
Ipc <- princomp(scale(Idata))
plot(Ipc)
autoplot(Ipc, loadings=T, loadings.label =T)
summary(Ipc)

Idata$F1 <- Ipc$scores[,1]
Idata$F2 <- Ipc$scores[,2]
Idata$F3 <- Ipc$scores[,3]

##### Identify the most relevant environemntal variables
Icor_clima <- as.data.frame(cor(Idata))
Icor_clima

Nf1 <- which.max(abs(Icor_clima[Icor_clima$F1!=1, "F1"])) 
Nf2 <- which.max(abs(Icor_clima[Icor_clima$F2!=1, "F2"])) 
Nf3 <- which.max(abs(Icor_clima[Icor_clima$F3!=1, "F3"])) 

row.names(Icor_clima[Nf1, ]) 
row.names(Icor_clima[Nf2, ]) 
row.names(Icor_clima[Nf3, ]) 

#Saving the climatic variables for Ipomoea

climatic_data <- coords_gen@data[,c(1:3)]
climatic_data$bio06 <- Idata$bio6_34      ## here, use the name of your selected variables 
climatic_data$bio18 <- Idata$bio18_34
climatic_data$bio16 <- Idata$bio16_34

str(climatic_data)  

write.table(climatic_data, "adapt_var_mapping/Dados_ambientais/climatic_data_ipomoea.txt")

