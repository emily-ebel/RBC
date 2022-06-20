# RBC

This repository contains data and code related to the manuscript "Common host variation drives malaria parasite fitness in healthy human red cells," now published at eLife (https://elifesciences.org/articles/69808).

All phenotypic data were normalized to correct for weekly batch effects. The raw and normalized data, as well as the normalization scripts, are provided. Example input data and scripts are also provided to run Lasso in R. 

Raw Illumina sequence data are available in the Sequence Read Archive (SRA) under BioProject accession number PRJNA683732.

Raw Sanger data and genotype calls for PIEZO1 E756del can be downloaded as a .zip file. This locus was amplified using primers described in Nguetse et al., 2020 (https://www.pnas.org/doi/abs/10.1073/pnas.1919843117).

Any questions? Contact Emily Ebel at ebel@stanford.edu.



## (1) Parasite fitness data

3D7_data_1.14.20.csv 

TH_data_1.14.20.csv 

parasite_normalization.R --> 

invasion_df.csv 

growth_df.csv 




## (2) RBC phenotype data  

###### (a) donor age, donor sex, Advia (CBC) ektacytometry summary statistics 

data_unnorm.csv 

median_normalization.R  --> 

data_norm.csv 


###### (b) osmotic fragility

raw_OF_data.csv

osmotic_fragility_normalization.R --> 

osmo50_unnorm.csv 

osmo50_norm.csv


###### (c) ektacytometry 

raw_ekta_data.csv

ektanormfactors.csv (from median_normalization.R)

sampleweek.csv (from median_normalization.R)

plot_normalized_ekta.R




## (3) LASSO files 

growth_geno.csv

invasion_geno.csv

Lasso target.R


## (4) Sanger data for PIEZO1 E756del

PIEZO1-E756del-Sanger.zip
