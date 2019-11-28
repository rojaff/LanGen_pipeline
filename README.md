# landgenomic_analysis

These scripts were developed to perform landscape genomic analysis including: filtering SNPs dataset (1), population genetic structure (2), generate climatic data (3), fine-scale spatial genetic structure (4), select potential candidate loci using pcadapt (5), RDA with Mahalanobis distance (6) and LFMM2 (7), comparing results with Venn Diagram (8), run RDA using only the candidate SNPs (9), run sPCA to make adaptive genetic maps (10) and predict genetic variation using RDA (11).

We suggest to run the scripts in the following order, because new VCF files will be used in the subsequent steps:

1 - Filtering

2 - Population Structure

3 - Climatic Data

4 - LOESS

5- PCAdapt

6 - GEA_RDA_Mahalanobis

7 - LFMM2

8 - Venn Diagram

9 - RDA_final

10 - sPCA

11 - RDA_predict
