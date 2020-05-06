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
#A. ESTIMATE POPULATION GENETIC DIVERSITY
#B. ESTIMATE POPULATION STRUCTURE

##2. REMOVE ANY OBJECT OR FUNCTION IN THE ENVIRONMENT:
rm(list=ls())

##3. INSTALL AND LOAD THE PACKAGES
library(r2vcftools)
library(LEA)
library(tess3r) 
library(seqinr)
library(adegenet)
library(vcfR)
library(SNPRelate)
library(dartR)

##4. INPUTS FOR THIS STEP:
#A. THE FILE ".VCF" AFTER FILTERING STEP 1.
#B. DOWNLOAD A VCF FILE AS EXAMPLE "Icavalcantei.vcf" FROM FIGSHARE: https://doi.org/10.6084/m9.figshare.6100004.v1
#C. KEEP VCF FILES IN A FOLDER NAMED "vcf" IN YOUR WORKING DIRECTORY.

##5. LOAD FUNCTIONS TO BE USED ON THIS STEP.
VCFsummary <- function(snps){
  Nid <- nrow(snps@meta)
  Nsnps <- length(snps@site_id)
  cat(paste(Nid, "individuals and", Nsnps, "SNPs."), sep="\n")
}

GenDiv <- function(vcf){
  HE <- Query(vcf, type="het")
  
  HE$HO <- (HE$N_SITES-HE$O.HOM.)/(HE$N_SITES) ## Observed heterozygosity (HO)
  error_ho <- qt(0.975,df=length(HE$HO)-1)*sd(HE$HO)/sqrt(length(HE$HO))
  ho <- mean(HE$HO)
  left_ho <- ho-error_ho
  right_ho <- ho+error_ho
  HO.df <- data.frame(He_O=ho, low_He_O=left_ho, up_He_O=right_ho)
  
  HE$HE <- (HE$N_SITES-HE$E.HOM.)/(HE$N_SITES) ## Expected heterozygosity (HE)
  error_he <- qt(0.975,df=length(HE$HE)-1)*sd(HE$HE)/sqrt(length(HE$HE))
  he <- mean(HE$HE)
  left_he <- he-error_he
  right_he <- he+error_he
  HE.df <- data.frame(He_E=he, low_He_E=left_he, up_He_E=right_he)
  
  error_f <- qt(0.975,df=length(HE$F)-1)*sd(HE$F)/sqrt(length(HE$F))
  f <- mean(HE$F)
  left_f <- f-error_f
  right_f <- f+error_f
  F.df <- data.frame(F=f, low_F=left_f, up_F=right_f)
  
  PI <- Query(vcf, type="site-pi")
  error_pi <- qt(0.975,df=length(PI$PI)-1)*sd(PI$PI)/sqrt(length(PI$PI))
  pi <- mean(PI$PI)
  left_pi <- pi-error_pi
  right_pi <- pi+error_pi
  PI.df <- data.frame(PI=pi, low_PI=left_pi, up_PI=right_pi)
  
  print(paste0("OUTPUT NUCLEOTIDE DIVERGENCE STATISTICS")) ## DIVERGENCE STATISTICS
  print(paste0("Observed heterozygosity (He_O)"))
  print(paste0("Expected heterozygosity (He_E)"))
  print(paste0("Coefficient of inbreeding (F)"))
  print(paste0("Measures nucleotide divergency on a per-site basis (PI)"))
  print(paste0("The 95% confidence interval - Lower limit (low_)"))
  print(paste0("The 95% confidence interval - Upper limit (up_)"))
  RES <- cbind(HO.df,  HE.df,  F.df,  PI.df)
  return(RES)
}

PlotK <- function(snmfproject){
  
  Mk1 <- mean(cross.entropy(snmfproject, K=1))
  SDk1 <- sd(cross.entropy(snmfproject, K=1))
  Mk2 <- mean(cross.entropy(snmfproject, K=2))
  SDk2 <- sd(cross.entropy(snmfproject, K=2))
  Mk3 <- mean(cross.entropy(snmfproject, K=3))
  SDk3 <- sd(cross.entropy(snmfproject, K=3))
  Mk4 <- mean(cross.entropy(snmfproject, K=4))
  SDk4 <- sd(cross.entropy(snmfproject, K=4))
  Mk5 <- mean(cross.entropy(snmfproject, K=5))
  SDk5 <- sd(cross.entropy(snmfproject, K=5))
  Mk6 <- mean(cross.entropy(snmfproject, K=6))
  SDk6 <- sd(cross.entropy(snmfproject, K=6))
  Mk7 <- mean(cross.entropy(snmfproject, K=7))
  SDk7 <- sd(cross.entropy(snmfproject, K=7))
  Mk8 <- mean(cross.entropy(snmfproject, K=8))
  SDk8 <- sd(cross.entropy(snmfproject, K=8))
  Mk9 <- mean(cross.entropy(snmfproject, K=9))
  SDk9 <- sd(cross.entropy(snmfproject, K=9))
  Mk10 <- mean(cross.entropy(snmfproject, K=10))
  SDk10 <- sd(cross.entropy(snmfproject, K=10))
  
  CE <- data.frame(K=c(1:10), Mean = c(Mk1,Mk2,Mk3,Mk4,Mk5,Mk6,Mk7,Mk8,Mk9,Mk10),
                   SD = c(SDk1,SDk2,SDk3,SDk4,SDk5,SDk6,SDk7,SDk8,SDk9,SDk10))
  
  library(ggplot2)
  
  ggplot(CE, aes(x=K, y=Mean)) + 
    geom_errorbar(aes(ymin=Mean-SD, ymax=Mean+SD), width=0.2)+
    geom_line() + 
    geom_point(size=4, shape=21, fill="red", color="darkred") + xlab("Number of ancestral populations") + ylab("Cross-entropy")+
    theme_bw() +
    theme(axis.text.x=element_text(angle=0, vjust=0.5), axis.text = element_text(size=15, face ="bold" , color = "black"), axis.title.x = element_text(size=15, face="bold", color="black"),axis.title.y = element_text(size=15, face="bold", color="black")) +
    scale_x_continuous(breaks = seq(0,10, by=2))
}

Best.run <- function(nrep, optimalK, p1, p2, p3, p4){
  ce1 = LEA::cross.entropy(p1, K = optimalK) # get the cross-entropy of each run for optimal K
  ce2 = LEA::cross.entropy(p2, K = optimalK) # get the cross-entropy of each run for optimal K
  ce3 = LEA::cross.entropy(p3, K = optimalK) # get the cross-entropy of each run for optimal K
  ce4 = LEA::cross.entropy(p4, K = optimalK) # get the cross-entropy of each run for optimal K
  AllProjects <- rbind(ce1, ce2, ce3, ce4)
  rownames(AllProjects) <- NULL
  AllProjects <- as.data.frame(AllProjects)
  AllProjects$Project <- c(rep(1, nrep), rep(2,nrep), rep(3, nrep), rep(4, nrep))
  Best_project <- AllProjects[AllProjects[, 1]==min(AllProjects[, 1]), 2]
  Best_runs <- AllProjects[AllProjects$Project==Best_project, ]
  Best_runs$Nrun <- 1:nrow(Best_runs)
  Best_run <- Best_runs[Best_runs[, 1]==min(Best_runs[, 1]), 3]
  print(paste0("Best run is: ", "project = ", Best_project, ", run = ", Best_run)) 
}

barplotK <- function(Qfile, Pop, Run_B){
  Q.matrix_1 <-LEA::Q(Qfile, K = Pop, run = Run_B)
  Q.matrix <- as.qmatrix(Q.matrix_1)
  barplot(Q.matrix, xlab = "Sampled individuals",
          ylab = "Ancestry coefficients",
          main = "", cex.axis = 1.5, cex.lab = 1)
}

fst = function(project, run = 1, K, ploidy = 2){
  library(LEA)
  l = dim(G(project, K = K, run = run))[1]
  q = apply(Q(project, K = K, run = run), MARGIN = 2, mean)
  if (ploidy == 2) {
    G1.t = G(project, K = K, run = run)[seq(2,l,by = 3),]
    G2.t = G(project, K = K, run = run)[seq(3,l,by = 3),]
    freq = G1.t/2 + G2.t}
  else {
    freq = G(project, K = K, run = run)[seq(2,l,by = 2),]}
  H.s = apply(freq*(1-freq), MARGIN = 1, FUN = function(x) sum(q*x))
  P.t = apply(freq, MARGIN = 1, FUN = function(x) sum(q*x))
  H.t = P.t*(1-P.t)
  return(1-H.s/H.t)
}

GIF <- function(project, run, K, fst.values){
  n = dim(Q(project, K, run))[1]
  fst.values[fst.values<0] = 0.000001
  z.scores = sqrt(fst.values*(n-K)/(1-fst.values))
  lambda = median(z.scores^2)/qchisq(1/2, df = K-1)
  print(lambda)
  return(lambda)
}

candidates <- function(alpha=0.05, adj.p.values){ ## alpha is significance level
  L = length(adj.p.values) ## Number of loci
  w = which(sort(adj.p.values) < alpha * (1:L) / L)
  candidates = order(adj.p.values)[w]
  print(length(candidates))
  return(candidates)
}

ManPlot <- function(adj.p.values, candidates, title){
  plot(-log10(adj.p.values), main=title, xlab = "Locus", cex = .7, col = "grey")
  points(candidates, -log10(adj.p.values)[candidates], pch = 19, cex = .7, col = "red")
}


#--------------------------------------------------------------------------------------
#                     Estimate the Genetic Diversity per Population 
#--------------------------------------------------------------------------------------

# In our example file there is only one population

## Load
snps_fil_ldF <-  vcfLink("vcf/ipomoea_filtered_ld_hwe_test2.vcf", overwriteID=T)
VCFsummary(snps_fil_ldF)  ## 115 individuals and 13604 SNPs.

## Overall genetic diversity
Overall <- GenDiv(snps_fil_ldF)
write.table(Overall, file="adapt_var_mapping/Pop_structure/Diversity_Overall_ipomoea.txt", quote=F, sep="/t")


#--------------------------------------------------------------------------------------
#                   Remove SNPs under selection using Fst Outlier approach 
#--------------------------------------------------------------------------------------           
              
##First we will carry out a population assignment analysis using sNMF to estimate number of K.

### Convert to geno object:
# You need to specify the output file. It will automatically subset the vcf file and assign it as a new object.
snps_fil_ldF <- Geno(snps_fil_ldF, output.file = "vcf/ipomoea_filtered_ld_hw.geno")
VCFsummary(snps_fil_ldF) ## 115 individuals and 13604 SNPs.

### Run SNMF (LEA) using different alpha
## Create different folders for each alpha (1, 10, 100, 500) and copy the GENO file in each folder

project_snmf1 = snmf("adapt_var_mapping/Pop_structure/Alpha1/ipomoea_filtered_ld_hw.geno", K = 1:10, rep = 10, alpha = 1, entropy = T, ploidy = 2, project = "new", CPU=4) #CPU=4 uses 4 CPUs 
project_snmf2 = snmf("adapt_var_mapping/Pop_structure/Alpha10/ipomoea_filtered_ld_hw.geno", K = 1:10, rep = 10, alpha = 10, entropy = T, ploidy = 2, project = "new", CPU=4) #CPU=4 uses 4 CPUs 
project_snmf3 = snmf("adapt_var_mapping/Pop_structure/Alpha100/ipomoea_filtered_ld_hw.geno", K = 1:10, rep = 10, alpha = 100, entropy = T, ploidy = 2, project = "new", CPU=4) #CPU=4 uses 4 CPUs 
project_snmf4 = snmf("adapt_var_mapping/Pop_structure/Alpha500/ipomoea_filtered_ld_hw.geno", K = 1:10, rep = 10, alpha = 500, entropy = T, ploidy = 2, project = "new", CPU=4) #CPU=4 uses 4 CPUs 

## To load the SNMF projects in a new R session (after quitting R), use:  project = load.snmfProject("Alpha1/Icavalcantei_filtered.snmfProject")
##This allows you to save time because you do not need to run sNMF again!
project1 = load.snmfProject("adapt_var_mapping/Pop_structure/Alpha1/ipomoea_filtered_ld_hw.snmfProject")
project2 = load.snmfProject("adapt_var_mapping/Pop_structure/Alpha10/ipomoea_filtered_ld_hw.snmfProject")
project3 = load.snmfProject("adapt_var_mapping/Pop_structure/Alpha100/ipomoea_filtered_ld_hw.snmfProject")
project4 = load.snmfProject("adapt_var_mapping/Pop_structure/Alpha500/ipomoea_filtered_ld_hw.snmfProject")

# summary of the project
summary(project1)
summary(project2)
summary(project3)
summary(project4)

##Cross-Entropy plot
plot(project1, lwd = 5, col = "red", pch=1)
plot(project2, lwd = 5, col = "red", pch=1)
plot(project3, lwd = 5, col = "red", pch=1)
plot(project4, lwd = 5, col = "red", pch=1)

##Select optimal K value
optimal_K = 1

## Bar plot of best K
ce1 = LEA::cross.entropy(project1, K = optimal_K) # get the cross-entropy of each run for optimal K
best1 = which.min(ce1) # select the run with the lowest cross-entropy

my.colors <- rainbow(optimal_K)
LEA::barchart(project1, K = optimal_K, run = best1, border = NA, space = 0, col = my.colors,
         xlab = "Individuals", ylab = "Ancestry proportions", main = "Ancestry matrix") -> bp
axis(1, at = 1:length(bp$order),labels = bp$order, las=1,cex.axis = .3)

## Bar plot of best K
ce2 = LEA::cross.entropy(project2, K = optimal_K) # get the cross-entropy of each run for optimal K
best2 = which.min(ce2) # select the run with the lowest cross-entropy

my.colors <- rainbow(optimal_K)
LEA::barchart(project2, K = optimal_K, run = best2, border = NA, space = 0, col = my.colors,
         xlab = "Individuals", ylab = "Ancestry proportions", main = "Ancestry matrix") -> bp
axis(1, at = 1:length(bp$order),labels = bp$order, las=1,cex.axis = .3)

## Bar plot of best K
ce3 = LEA::cross.entropy(project3, K = optimal_K) # get the cross-entropy of each run for optimal K
best3 = which.min(ce3) # select the run with the lowest cross-entropy

my.colors <- rainbow(optimal_K)
LEA::barchart(project3, K = optimal_K, run = best3, border = NA, space = 0, col = my.colors,
         xlab = "Individuals", ylab = "Ancestry proportions", main = "Ancestry matrix") -> bp
axis(1, at = 1:length(bp$order),labels = bp$order, las=1,cex.axis = .3)

## Bar plot of best K
ce4 = LEA::cross.entropy(project4, K = optimal_K) # get the cross-entropy of each run for optimal K
best4 = which.min(ce3) # select the run with the lowest cross-entropy

my.colors <- rainbow(optimal_K)
LEA::barchart(project4, K = optimal_K, run = best4, border = NA, space = 0, col = my.colors,
         xlab = "Individuals", ylab = "Ancestry proportions", main = "Ancestry matrix") -> bp
axis(1, at = 1:length(bp$order),labels = bp$order, las=1,cex.axis = .3)

## Select best run (lowest cross-entropy)
Best.run(nrep=10, optimalK=1, p1=project1, p2=project2, p3=project3, p4=project4)

## Final Barplot
barplotK(Qfile=project3, Pop = 1, Run_B = 6)

## Add population IDS to VCF file
Qmat <- Q(project3, run=6, K=optimal_K)

popIds = apply(Qmat, 1, which.max)
snps_fil_ldF@meta$PopID_snmf <- popIds
snps_fil_ldF@meta

## Save new vcf file and pop ID file
Save(snps_fil_ldF, "vcf/ipomoea_filtered_ld_hw_pops_snmf.vcf") 
              
## Compute the FST statistics using best run
project=project3
run=6
K=2

FST = fst(project, run, K) # you need at least 2 populations for a population-based test, so K>1.

## Compute the GIF
lambda <- GIF(project, run, K, fst.values=FST)
lambda

## Compute adjusted p-values from the combined z-scores and plot histogram of p-values
n = dim(Q(project, run, K))[1]
z.scores = sqrt(FST*(n-K)/(1-FST))
adj.p.values = pchisq(z.scores^2/lambda, df = K-1, lower = FALSE)
hist(adj.p.values, col = "red")

## Test different lambda values and plot histogram of p-values
adj.p.values = pchisq(z.scores^2/1.5, df = K-1, lower = FALSE) ## it is the best, but is still strange
hist(adj.p.values, col = "green")

## Candidate loci
## FDR control: Benjamini-Hochberg at level q
C_fst <- candidates(alpha=0.05, adj.p.values)

## Manhatan plot
ManPlot(adj.p.values, C_fst,"Fst")

## Exclude candidate FST outlier

snps_fil_ldF_candidate <- Subset(snps_fil_ldF, sites=C_fst)
snps_fil_ldF_candidate@site_id ## These are all the candidate SNPs
Chrom(snps_fil_ldF_candidate)

candidates_fst <- snps_fil_ldF_candidate@site_id
All_snp <- snps_fil_ldF@site_id
N_snp <- All_snp[!(All_snp %in% candidates_fst)] ###Exclude all candidate loci

snps_fil_ldF_neutral <- Subset(snps_fil_ldF, sites=N_snp)
VCFsummary(snps_fil_ldF_neutral)

length(N_snp)
length(snps_fil_ldF@site_id)-length(C_fst) 
length(snps_fil_ldF_neutral@site_id)

## Save neutral snp dataset
Save(snps_fil_ldF_neutral, "vcf/ipomoea_filtered_ld_hw_neutral.vcf")

              
#--------------------------------------------------------------------------------------
#                  Population Assignment Analysis using sNMF 
#-------------------------------------------------------------------------------------- 
              
#Now we will carry out population assignment again, but using only neutral loci dataset.

## Load
snps_fil_ldF_neutral <- vcfLink("vcf/ipomoea_filtered_ld_hw_neutral.vcf") 

## Convert to geno object:
## You need to specify the output file. It will automatically subset the vcf file and assign it as a new object.
project.snmf = NULL

snps_fil_ldF_neutral <- Geno(snps_fil_ldF_neutral, output.file = "vcf/ipomoea_filtered_ld_hw_neutral.geno")
VCFsummary(snps_fil_ldF)
VCFsummary(snps_fil_ldF_neutral) # 115 individuals and 13167 SNPs.

## Run SNMF (LEA) using different alpha
## Create different folders for each alpha (1, 10, 100, 500) and copy the GENO file in each folder
project_snmf1n = snmf("adapt_var_mapping/Pop_structure/Alpha1n/ipomoea_filtered_ld_hw_neutral.geno", K = 1:10, rep = 10, alpha = 1, entropy = T, ploidy = 2, project = "new", CPU=4) #CPU=4 uses 4 CPUs 
project_snmf2n = snmf("adapt_var_mapping/Pop_structure/Alpha10n/ipomoea_filtered_ld_hw_neutral.geno", K = 1:10, rep = 10, alpha = 10, entropy = T, ploidy = 2, project = "new", CPU=4) #CPU=4 uses 4 CPUs 
project_snmf3n = snmf("adapt_var_mapping/Pop_structure/Alpha100n/ipomoea_filtered_ld_hw_neutral.geno", K = 1:10, rep = 10, alpha = 100, entropy = T, ploidy = 2, project = "new", CPU=4) #CPU=4 uses 4 CPUs 
project_snmf4n = snmf("adapt_var_mapping/Pop_structure/Alpha500n/ipomoea_filtered_ld_hw_neutral.geno", K = 1:10, rep = 10, alpha = 500, entropy = T, ploidy = 2, project = "new", CPU=4) #CPU=4 uses 4 CPUs 

## This allows you to save time because you do not need to run SNMF again!
project1n = load.snmfProject("adapt_var_mapping/Pop_structure/Alpha1n/ipomoea_filtered_ld_hw_neutral.snmfProject")
project2n = load.snmfProject("adapt_var_mapping/Pop_structure/Alpha10n/ipomoea_filtered_ld_hw_neutral.snmfProject")
project3n = load.snmfProject("adapt_var_mapping/Pop_structure/Alpha100n/ipomoea_filtered_ld_hw_neutral.snmfProject")
project4n = load.snmfProject("adapt_var_mapping/Pop_structure/Alpha500n/ipomoea_filtered_ld_hw_neutral.snmfProject")

## Summary of the project
summary(project1n)
summary(project2n)
summary(project3n)
summary(project4n)

## Cross-Entropy plot
plot(project1n, lwd = 5, col = "red", pch=1)
plot(project2n, lwd = 5, col = "red", pch=1)
plot(project3n, lwd = 5, col = "red", pch=1)
plot(project4n, lwd = 5, col = "red", pch=1)

## Cross-Entropy plot with standard deviation error bars
PlotK(project1n)
PlotK(project2n)
PlotK(project3n)
PlotK(project4n)

## Select optimal K value
optimal_K = 1

## Bar plot of best K
ce1n = cross.entropy(project1n, K = optimal_K) # get the cross-entropy of each run for optimal K
best1n = which.min(ce1n) # select the run with the lowest cross-entropy

my.colors <- rainbow(optimal_K)
LEA::barchart(project1n, K = optimal_K, run = best1n, border = NA, space = 0, col = my.colors,
         xlab = "Individuals", ylab = "Ancestry proportions", main = "Ancestry matrix") -> bp
axis(1, at = 1:length(bp$order),labels = bp$order, las=1,cex.axis = .3)


ce2n = cross.entropy(project2n, K = optimal_K) # get the cross-entropy of each run for optimal K
best2n = which.min(ce2n) # select the run with the lowest cross-entropy
LEA::barchart(project2n, K = optimal_K, run = best2n, border = NA, space = 0, col = my.colors,
         xlab = "Individuals", ylab = "Ancestry proportions", main = "Ancestry matrix") -> bp
axis(1, at = 1:length(bp$order),labels = bp$order, las=1,cex.axis = .3)

ce3n = cross.entropy(project3n, K = optimal_K) # get the cross-entropy of each run for optimal K
best3n = which.min(ce3n) # select the run with the lowest cross-entropy
LEA::barchart(project3n, K = optimal_K, run = best3n, border = NA, space = 0, col = my.colors,
         xlab = "Individuals", ylab = "Ancestry proportions", main = "Ancestry matrix") -> bp
axis(1, at = 1:length(bp$order),labels = bp$order, las=1,cex.axis = .3)

ce4n = cross.entropy(project4n, K = optimal_K) # get the cross-entropy of each run for optimal K
best4n = which.min(ce4n) # select the run with the lowest cross-entropy
LEA::barchart(project4n, K = optimal_K, run = best4n, border = NA, space = 0, col = my.colors,
         xlab = "Individuals", ylab = "Ancestry proportions", main = "Ancestry matrix") -> bp
axis(1, at = 1:length(bp$order),labels = bp$order, las=1,cex.axis = .3)

## Select best run (lowest cross-entropy)
Best.run(nrep=10, optimalK=1, p1=project1n, p2=project2n, p3=project3n, p4=project4n)

## Barplot
barplotK(Qfile=project4n, Pop= 1, Run_B = 1)

## Add population IDS to VCF file
Qmat <- Q(project4n, run=1, K=1)

## Add Q-matrix and pop IDs to vcf file
cbind(snps_fil_ldF_neutral@meta, as.data.frame(Qmat))
popIds = apply(Qmat, 1, which.max)
snps_fil_ldF_neutral@meta$PopLEA <- popIds

## Save new vcf file and pop ID file
Save(snps_fil_ldF_neutral, "vcf/ipomoea_filtered_ld_hw_neutral_LEA.vcf") 
write.table(cbind(snps_fil_ldF_neutral@meta, as.data.frame(Qmat)), file="adapt_var_mapping/Pop_structure/ancestry_ipomoea.txt")

## Save csv file with pop ID of individuals assessed using LEA
write.csv(as.data.frame(snps_fil_ldF_neutral@meta)[,-(5)], file="adapt_var_mapping/Pop_structure/ipomoea__neutral_popsLEA.csv")


#--------------------------------------------------------------------------------------
#                  Population Assignment Analysis using DAPC 
#--------------------------------------------------------------------------------------              
              
#Here we will use other population clustering methods to assign individual to different population. This method is DAPC

## Load
vcf <- read.vcfR("vcf/ipomoea_filtered_ld_hw_neutral.vcf", verbose = FALSE)

## TRANSFORMATION "VCF" TO "GENLIGHT"
my_genind_pop <- vcfR2genlight(vcf) # It's used for alleles numbers counted by individuals
my_genind_pop

## Identification of the clusters (We specify that we want to evaluate up to k = 10 groups (max.n.clust=40):)
grp <- find.clusters(my_genind_pop, max.n.clust=10) # First run without informing these values (n.pca = 100 PCs  (All) and n.clust = 1 Cluster) 

## The ’best’ BIC is often indicated by an below in the curve of BIC values as a function of k
              
## Plot the best K
plot(grp$Kstat, type="o", xlab="Number of clusters (K)", ylab="BIC",
     col="blue", main="")

## Group (DPCA) of the first 10 individuals
head(grp$grp, 10)

## Size of Groups (DPCA) 
grp$size

## DAPC
dapc_p <- dapc(my_genind_pop, grp$grp)
dapc_p <- dapc(my_genind_pop, grp$grp, n.da=1, n.pca=100) #100 PCs and 1 discriminant functions saved
summary(dapc_p)

## Basic scatterplots (1)
scatter(dapc_p, posi.da="bottomleft", bg="white", pch=19:20)

## Plotted in a STRUCTURE-like 
compoplot(dapc_p)

              
#--------------------------------------------------------------------------------------
#                        Estimate FST among Genetic Clusters 
#--------------------------------------------------------------------------------------             

## Define clusters as POP_ID from sMNF
my_genind_pop@pop = as.factor(snps_fil_ldF_neutral@meta$PopID_smnf)

## Define clusters as POP_ID from DAPC
my_genind_pop@pop = as.factor(grp$grp)
  
## FST              
gl.fst.pop(my_genind_pop, nboots = 100, percent = 95, nclusters = 1)


#--------------------------------------------------------------------------------------
#                     Estimate the Genetic Diversity per Population 
#--------------------------------------------------------------------------------------

## Genetic diversity with neutral loci

## Load
snps_fil_ldF_neutral <- vcfLink("vcf/ipomoea_filtered_ld_hw_neutral.vcf", overwriteID=T)
VCFsummary(snps_fil_ldF_neutral)

## Overall genetic diversity
Overall <-GenDiv(snps_fil_ldF_neutral)
write.table(Overall, file="adapt_var_mapping/Pop_structure/Diversity_Overall_neutral_ipomoea.txt", quote=F, sep="/t")

## END
