# Melbourne Veterinary Collection (MVC)
Scripts to run analysis on *E. coli* genomes from the Melbourne Veterinary Collection associated with the publication ["Genomic and temporal trends in canine ExPEC reflect those of human ExPEC"](https://doi.org/10.1128/spectrum.01291-22)   by Elankumaran *et al*, 2022

## Overview
The following repository allows conscientious readers of the manuscript to reproduce all the data processing, statistics and figures presented in the paper.

It comprises two directories: __`scripts`__ and __`data`__ (the contents of which should be self explanatory) and generates a __`data/outputs`__ folder with __`figures`__ and __`data`__ subdirectories.


## Installation
### Software requirements
These scripts are currently functional on mac OS Big Sur 11.5.2 using RStudio 1.4.1106 and R version 4.0.5. We cannot guarantee they will work on other distributions of R or RStudio. Your OS should not be an issue provided you use these versions of R and RStudio though.

### Packages required
You will need to install the following packages and versions to work with the scripts:
- tidytree_0.3.5     
- data.table_1.14.0  
- tidyverse_1.3.1    
- magrittr_2.0.1     
- plasmidmapR_0.1.0  
- abricateR_0.1.1   
- RColorBrewer_1.1-2 
- ggtree_3.1.0       
- pheatmap_1.0.12    
- reshape2_1.4.4     
- ggpubr_0.4.0       
- ggplot2_3.3.5     
- tibble_3.1.5       
- purrr_0.3.4        
- readr_2.0.1        
- stringr_1.4.0      
- forcats_0.5.1      
- tidyr_1.1.4       
- dplyr_1.0.7       


### Issues with ggtree
There have recently been some issues with ggtree in the way it interacts with dplyr. The solution is to install the latest version of ggtree directly from github instead of via BiocManager. You can do this in the console on RStudio with the __`remotes`__ package like so:
```
install.packages("remotes")
remotes::install_github("YuLab-SMU/ggtree")
```

## Usage
Clone this repository
```
git clone https://github.com/CJREID/MVC.git
cd MVC
pwd
```
Open the MVC_analyis.R script in a text editor or RStudio and set the variable __`wrkdir`__ on line 23 to the output of __`pwd`__ above and save the script.

Run the MVC_analyis.R script and your __`outputs`__ folder will populate with figures and tables.

## Outputs
### Figures
1. Figure 1. Characteristics of the genome collection
2. Figure 2. Core gene maximum-likelihood phylogenetic tree with metadata (Legends manually edited for publication)
3. Figure 3. F plasmid carriage by phylogroup and ST
4. Figure 4. Distribution of antimicrobial resistance genes and virulence genes by phylogroup and ST

### Supplementary
#### Tables
1. Table S1. Metadata, accession numbers and gene screening results for 377 canine E. coli isolates used in this study


#### Figures
1. Fig S1. Presence/absence of ARGs mapped to core gene phylogeny
2. Fig S2. Presence/absence of VAGs mapped to core gene phylogeny
3. Fig S3. Presence/absence of plasmid-associated genes mapped to core gene phylogeny
4. Fig S4. Alignment of all sequences to pUTI89 plasmid mapped to core gene phylogeny
5. Fig S5. Alignment of all sequences to pCERC4 ColV plasmid mapped to core gene phylogeny
