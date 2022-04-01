####------MVC ANALYSIS------####
####-----PACKAGES-----####
# Load required packages
library(dplyr)
library(tidyr)
library(forcats)
library(stringr)
library(readr)
library(purrr)
library(tibble)
library(ggplot2)
library(ggpubr)
library(reshape2)
library(pheatmap)
library(ggtree)
library(RColorBrewer)
library(abricateR)
library(plasmidmapR)

####-----SETUP-----####
# SET WORKING DIRECTORY - ADD THE PATH TO THE `ST58_project` REPOSITORY ON YOUR COMPUTER
# SCRIPT WILL NOT WORK IF YOU DO NOT CHANGE THIS
wrkdir <- "/Volumes/126451/WORK/Colleagues/Paarthiphan/MVC Collection/data/MVC_analysis/"
setwd(wrkdir)

# FUNCTION FOR SUBSETTING
'%notin%' <- Negate('%in%')

# Make output directory
if (dir.exists("outputs")){
} else {
  dir.create("outputs")
  dir.create("outputs/data")
  dir.create("outputs/figures")
}

####-----IMPORT DATA-----####
#Get input filenames
infiles <-c(list.files("data/pipelord_results/", pattern = "\\.txt"))

# Read in and name data.frames with whatever comes before .txt
for (file in infiles){
  indir <- c("data/pipelord_results/")
  f <- read_delim(paste(indir, file, sep = "/"), delim = "\t", col_names = TRUE, trim_ws = TRUE)
  assign(paste(substr(file, 1, nchar(file)-4), sep = ""), f)
  rm(f)
}

####------TREE-----####
## UNROOTED CORE
unrooted_core_tree <- read.tree("data/tree/core_gene_alignment.aln.treefile")

## ROOTED CORE
mproot_core <- phytools::midpoint.root(unrooted_core_tree)

####------METADATA-----####
names <- read.delim("data/meta/MVC_names.csv", sep = ",")
names <- rename(names, c(Plate = Initial_sample_name, Name = Current_sample_name))
names$Plate <- gsub("_P[0-9]", "", names$Plate)

raw_meta <- read.delim("data/meta/raw_data_020518.txt", na.strings = ".")

meta <- raw_meta %>% filter(Genus == "Escherichia", Species.type == "Canine")

meta <- meta %>% mutate(Specimen = case_when(Specimen == "skin and subcutis" | Specimen == "ear" | Specimen == "mammary" ~ "soft tissue",
                                             Specimen == "reproductive" | Specimen == "urinary" ~ "urogenital",
                                             Specimen == "Systemic" | Specimen == "systemic or large cavities" | Specimen == "cardiovascular" |
                                               Specimen == "musculo-skeletal" | Specimen == "lymphatic" | Specimen == "nervous" |
                                               Specimen == "respiratory"  ~ "general",
                                             TRUE ~ Specimen),
                        Specimen = str_to_title(Specimen))

meta <- left_join(names, meta)

# Filter metadata to tree strains only
meta <- meta %>% filter(Name %in% mproot_core$tip.label)

####------PHYLOGROUP-----####
raw_clermont <- read.delim("data/meta/MVC_ezclermont.txt", header = FALSE, col.names = c("Name", "Phylogroup"))

raw_clermont <- raw_clermont %>% 
  mutate(simple_phylo = gsub("U$", "Cryptic", Phylogroup)) %>% 
  mutate(simple_phylo = gsub("U/", "", simple_phylo)) %>%
  mutate(simple_phylo = str_to_sentence(simple_phylo)) %>% select(-Phylogroup, Phylogroup = simple_phylo)

meta <- left_join(meta, raw_clermont)

####------SEROTYPE-----####
raw_sero <- read_csv("data/meta/MVC_serotypes_summary.csv")

# Keep raw sero data for Table S1 later
full_OH <- raw_sero %>% select(Name, OH_type)

# Split into types with 10 or more, and less than 10 representatives
major_sero <- raw_sero %>% select(Name, OH_type) %>% group_by(OH_type) %>% filter(n() >= 10)
minor_sero <- raw_sero %>% select(Name, OH_type) %>% group_by(OH_type) %>% filter(n() < 10)

# Designate types with less than 10 as "Other"
minor_sero$OH_type <- "Other"

# Stick them back together
sero <- rbind(major_sero, minor_sero)
sero <- left_join(sero, raw_sero %>% select(-OH_type))

meta <- left_join(meta, sero)
meta <- left_join(meta, raw_sero)

####------MLST-----####
full_mlst <- mlst %>% select(Name = name, ST, -scheme) %>% mutate(ST = gsub("\\s.*", "", ST))

mlst <- full_mlst %>% mutate(ST = paste0("ST", ST), ST = gsub("ST-", "Other", ST))

# Split into types with 10 or more, and less than 10 representatives
major_ST <- mlst %>% group_by(ST) %>% filter(n() >= 5)
minor_ST <- mlst %>% group_by(ST) %>% filter(n() < 5)

# Designate types with less than 10 as "Other"
minor_ST$ST <- "Other"

# Stick them back together
mlst <- rbind(major_ST, minor_ST)

meta <- left_join(meta, mlst)

####------dfrA5-IS26 SIGNATURE----------####
# Process the dfrA5-IS26 signature data
abricateR("data/meta/MVC.848.summary.txt", output = "IS26_signature",identity = 99, length = 90, writecsv = FALSE)
signature_pos <- `IS26_signature_simple_summary_N99L90` %>% rename(IS26_signature = `848_gb|CP014489.1|:130176-132039`, Name = name) %>% select(-ColV)
signature_neg <- meta %>% select(Name) %>% filter(Name %notin% signature_pos$Name)
signature_neg$IS26_signature <- rep(0, nrow(signature_neg))
IS26_signature <- rbind(signature_pos, signature_neg) %>% mutate(IS26_signature = case_when(IS26_signature =="1"  ~ "Yes",
                                                                                            IS26_signature =="0"  ~ "No"))

####-----PROCESS ABRICATE DATA-----####
#Load paths to files needed for abricateR
abricate_path <- "data/pipelord_results/genotype.txt"
pointfinder_path <- "data/pipelord_results/pointfinder.txt"
pMLST_data <- "data/pipelord_results/pMLST.txt"

#Provide output names
output_name <- "MVC"

#Run abricateR
abricateR::abricateR(
  abricate_in = abricate_path,
  output = output_name,
  identity = 90,
  length = 95,
  output_directory = output_dir,
  writecsv = FALSE,
  pointfinder_data = pointfinder_path,
  pMLST_data = pMLST_data
)

####-----PROCESS POINTFINDER DATA-----####
# Read in phenotypic prediction data and make strains resistant to both NA and Cipro "Yes" for FQR
fqr_pheno <-
  read_tsv("data/pipelord_results/MVC.pointfinder.prediction.txt") %>%
  select(Name = `Sample ID`, `Nalidixic Acid` = `NALIDIXIC ACID`, `Ciprofloxacin`=`CIPROFLOXACIN`) %>% 
  mutate(FQR = case_when(`Nalidixic Acid` == 1 & Ciprofloxacin ==1 ~"Yes", TRUE ~ "No")) %>%
  select(Name, FQR)

####-----PROCESS PLASMID ALIGNMENTS-----####
# Get tree path and abricate alignment data
path_to_tree <- "data/tree/core_gene_alignment.aln.treefile"
path_to_abricate <- "data/pipelord_results/MVC.pCERC4.abricate.tab"
plasrefname <- "pCERC4"

pCERC4_ref_length <- read_delim("data/pipelord_results/MVC.pCERC4.abricate.tab", n_max = 1) %>% 
  mutate(COVERAGE = as.integer(gsub("^.*\\/", "", COVERAGE))) %>% pull(COVERAGE)

plasmid_mapR(path_to_abricate = path_to_abricate,
            plasmid_reference_name = plasrefname,
            output_directory = NULL,
            min_hit_id = 90,
            min_hit_length = 0.5,
            writecsv = FALSE,
            path_to_tree = path_to_tree)

# Generate dataframe with a summary statistic for how many bins were covered at >90%ID
pCERC4_ID <- as.data.frame(rowSums(pCERC4_binned_hits))

# Add names column
pCERC4_ID$Name <- rownames(pCERC4_binned_hits)

# Rename columns
colnames(pCERC4_ID) <- c("pCERC4_ID","Name")

# Convert value to a percentage
pCERC4_ID$pCERC4_ID <- round((pCERC4_ID$pCERC4_ID/pCERC4_ref_length) * 100)

path_to_abricate <- "data/pipelord_results/MVC.pUTI89.abricate.tab"
plasrefname <- "pUTI89"

pUTI89_ref_length <- read_delim("data/pipelord_results/MVC.pUTI89.abricate.tab", n_max = 1) %>% 
  mutate(COVERAGE = as.integer(gsub("^.*\\/", "", COVERAGE))) %>% pull(COVERAGE)

plasmid_mapR(path_to_abricate = path_to_abricate,
            plasmid_reference_name = plasrefname,
            output_directory = NULL,
            min_hit_id = 90,
            min_hit_length = 0.5,
            writecsv = FALSE,
            path_to_tree = path_to_tree)

# Generate dataframe with a summary statistic for how many bins were covered at >90%ID
pUTI89_ID <- as.data.frame(rowSums(pUTI89_binned_hits))

# Add names column
pUTI89_ID$Name <- rownames(pUTI89_binned_hits)

# Rename columns
colnames(pUTI89_ID) <- c("pUTI89_ID","Name")

# Convert value to a percentage
pUTI89_ID$pUTI89_ID <- round((pUTI89_ID$pUTI89_ID/pUTI89_ref_length) * 100)

####-----UPDATE METADATA-----####
#Select ColV data to add to meta
colv <- MVC_simple_summary_N90L95 %>% select(Name = name, ColV, IncF_RST)

# Join ColV and F-type to meta
meta <- left_join(meta, colv)

# Join pCERC4 alignment percentage
meta <- left_join(meta, pCERC4_ID)

# Categorise ColV
meta <- meta %>% mutate(ColV = case_when(ColV == "0" ~ "No", ColV == "1" ~ "Yes"))

# Categorise seqeunces that are F29:A-:B10 carriers or >90% coverage of pUTI89 as pUTI89+
meta <- left_join(meta, pUTI89_ID) %>%
  mutate(pUTI89 = case_when(pUTI89_ID >= 90 | IncF_RST == "F29:A-:B10" ~ "Yes", TRUE ~ "No"))

# Reduce IncF RST to only common types and re-add as a separate column
simple_F <- meta %>% select(Name, IncF_RST) 

# Get full F types for Table S1
full_F <- MVC_simple_summary_N90L95 %>% select(Name = name,`F Plasmid` =  IncF_RST)

# Split into types with 10 or more, and less than 10 representatives
major_F <- simple_F %>% group_by(IncF_RST) %>% filter(n() >= 10)
minor_F <- simple_F %>% group_by(IncF_RST) %>% filter(n() < 10)

# Designate types with less than 10 as "Other"
minor_F$IncF_RST <- "Other"

# Stick them back together
simple_F <- rbind(major_F, minor_F)
simple_F <- simple_F %>% rename(`F Plasmid` = IncF_RST)

simple_F$`F Plasmid` <- simple_F$`F Plasmid` %>% recode(`F-:A-:B-` = "No F Plasmid")

meta <- left_join(meta, simple_F)

# Read in fimH data and join to meta
fimH <- read_csv("data/meta/MVC_fimH.csv") %>% mutate(fimH = gsub("\\*", "Unknown", fimH))
meta <- left_join(meta, fimH)

# Join predicted fluoroquinolone resistance phenotype
meta <- left_join(meta, fqr_pheno)

# Join IS26 signature preence/absence data
meta <- left_join(meta, IS26_signature) 

# Extract working names for collection
working_names <- as.vector(meta$Name)

# Create 'geno_meta' df with all meta and gene screening data
geno_meta <- left_join(meta %>% select(Name, Specimen, Phylogroup, ST, OH_type, ColV, pUTI89, `F Plasmid`), 
                       MVC_simple_summary_N90L95 %>% select(-ColV, -IncF_RST), by = c("Name" = "name"))

# Define gene columns as those that are integers
gene_cols <- names(geno_meta %>% select(where(is.integer)))

# Recode multiple hits as a single hit
geno_meta <- geno_meta %>% mutate(across(where(is.integer), ~ gsub(2, 1, .x)))
geno_meta <- geno_meta %>% mutate(across(where(is.integer), ~ gsub(3, 1, .x)))
geno_meta <- geno_meta %>% mutate(across(where(is.integer), ~ gsub(4, 1, .x)))

# Convert back to integer (gsub makes everything character class)
geno_meta <- geno_meta %>% mutate(across(all_of(gene_cols), as.integer))

# Extract metadata column names for later use with the geno_meta dataframe
meta_cols <- names(meta %>% select(Name, Specimen, Phylogroup, ST, OH_type, ColV, pUTI89, `F Plasmid`))

####-----SCOARY-----####
# # Select and export categorical data
# # Select relevant columns
# scoary_ST <- meta %>% select(Name, ST) %>% mutate(value = rep(1, nrow(meta)))
# 
# # Cast into wide format
# scoary_ST <- scoary_ST %>% reshape2::dcast(Name ~ ST) %>% select(Name, ST372, ST73, ST127)
# 
# # Recode NAs to zeroes
# scoary_ST[is.na(scoary_ST)] <- 0
# 
# # Write file for Scoary
# write_csv(scoary_ST, "outputs/data/MVC.ST.scoary.csv")
# 
# # Generate file list of Scoary results for BAP groups
# filelist <- list.files(path="data/scoary", pattern="*.csv", full.names = TRUE)
# 
# # Read in thefiles that actually have data - used file size as a proxy for this
# for (f in 1:length(filelist)){
#   if(file.size(filelist[f]) > 229) {
#     assign(paste0(str_extract(filelist[f], "ST\\d{1,3}"), "_scoary"), read_csv(filelist[f]))
#   }
# }
# 
# # Combine into list
# scoary_list <- mget(ls(pattern = "ST\\d{1,3}_scoary"))
# 
# # Processing loop for each file
# for (f in 1:length(scoary_list)){
#   data <- scoary_list[[f]]
#   name <- str_extract(names(scoary_list[f]), "BAP_\\d")
#   # Filter hypotheticals and combine gene and non-unique gene name to create a unique name that splits paralogs
#   data <- data %>% 
#     mutate(BAP = rep(name, nrow(data)), 
#            Gene_Unique = paste(`Non-unique Gene name`, Gene, sep = "_")) %>% 
#     filter(!grepl("hypothetical", Annotation), Benjamini_H_p < 1E-30) %>% 
#     select(Gene_Unique, everything())
#   
#   # Get rid of NAs
#   data$Gene_Unique <- gsub("NA_", "", data$Gene_Unique)
#   
#   # Split into over and under-represented ColV clade genes based on the Odds Ratio
#   over_rep <- data %>% filter(Odds_ratio > 1 | Odds_ratio == "inf") %>%
#     # slice(0:20) %>%
#     dplyr::select(-Sensitivity,-Specificity,-Bonferroni_p) %>%
#     dplyr::rename(
#       Pos_present = Number_pos_present_in,
#       Neg_present = Number_neg_present_in,
#       Pos_absent = Number_pos_not_present_in,
#       Neg_absent = Number_neg_not_present_in
#     ) 
#   
#   under_rep <- data %>% filter(Odds_ratio < 1) %>%
#     # slice(0:20) %>%
#     dplyr::select(-Sensitivity,-Specificity,-Bonferroni_p) %>%
#     dplyr::rename(
#       Pos_present = Number_pos_present_in,
#       Neg_present = Number_neg_present_in,
#       Pos_absent = Number_pos_not_present_in,
#       Neg_absent = Number_neg_not_present_in
#     ) 
#   # Assign names to the outputs
#   assign(paste(str_extract(names(scoary_list[f]), "ST\\d{1,3}"), "over", sep = "_"), over_rep)
#   assign(paste(str_extract(names(scoary_list[f]), "ST\\d{1,3}"), "under", sep = "_"), under_rep)
#   
#   # Remove unnecessary objects
#   rm(data)
#   rm(over_rep)
#   rm(under_rep)
#   rm(name)
# }

####------COLOURS-----####
RColorBrewer::display.brewer.all()

# Source
source_vars <- unique(meta$Specimen)
source_clrs <- colorRampPalette(brewer.pal(12, "Paired"))(length(unique(meta$Specimen)))
names(source_clrs) <- source_vars
source_clrs["Other"] <- "#dbdbdb"

# ST
ST_vars <- unique(meta$ST)
ST_clrs <- colorRampPalette(brewer.pal(12, "Paired"))(length(unique(meta$ST)))
names(ST_clrs) <- sort(ST_vars)
ST_clrs["Other"] <- "#dbdbdb"

# Phylogroups
phylo_vars <- unique(meta$Phylogroup)
phylo_clrs <- brewer.pal(9, "Paired")
names(phylo_clrs) <- sort(phylo_vars)

# F-types
F_vars <- unique(meta$`F Plasmid`)
F_clrs <- brewer.pal(6, "Dark2")
names(F_clrs) <- F_vars 
F_clrs["Other"] <- "#dbdbdb"
F_clrs["No F Plasmid"] <- "#FFFFFF"

# Serotypes
OH_vars <- unique(meta$OH_type)
OH_clrs <- brewer.pal(10, "Paired")
names(OH_clrs) <- OH_vars
OH_clrs["Other"] <- "#dbdbdb"

# ColV
colv_clrs <- c("ColV-Pos" = "#17e6ae", "ColV-Neg" = 'white')

# pUTI89
puti89_clrs <- c("pUTI89-Pos" = "#E7298A", "pUTI89-Neg" = 'white')

# Combine for tree
tree_vars <- c(source_clrs, ST_clrs, OH_clrs, F_clrs, colv_clrs, puti89_clrs)

# Select phylogroup data for tips
phylo_tips <- meta %>% select(Name, Phylogroup)

# Select data for tree visualisation
tree_meta <- meta %>% select(Specimen, ST, Serotype = OH_type, `F Plasmid`, ColV, pUTI89)

# Recode ColV and pUTI89 values for clarity on the figure
tree_meta$ColV <- tree_meta$ColV %>% recode(Yes = "ColV-Pos", No = "ColV-Neg")
tree_meta$pUTI89 <- tree_meta$pUTI89 %>% recode(Yes = "pUTI89-Pos", No = "pUTI89-Neg")

# Add rownames to data for gheatmap
nems <- meta$Name
rownames(tree_meta) <- nems
rownames(phylo_tips) <- nems

####------HEATMAP PROCESSING-----####
# Split ABRicate gene hits into their functional groups
# Get all the hits from CARD database and intI1 and intI2. Fix all the messy names. Filter out genes present in >90% of isolates. These are housekeeping
# genes that sometimes mutate to confer AMR phenotypes but we are only concerned wiht acquired resistance genes
args <- geno_meta %>%
  select(all_of(meta_cols), starts_with("card"), contains("intI")) %>%
  rename_with(~ gsub("card_", "", .x, fixed = TRUE)) %>%
  rename_with(~ gsub("Escherichia_coli_", "", .x, fixed = TRUE)) %>%
  rename_with(~ gsub("_beta-lactamase", "", .x, fixed = TRUE)) %>%
  rename_with(~ gsub("(", "_", .x, fixed = TRUE)) %>%
  rename_with(~ gsub(")", "", .x, fixed = TRUE)) %>%
  rename_with(~ gsub("PC1__", "PC1_", .x, fixed = TRUE)) %>%
  rename_with(~ gsub("EC_custom_intI1.*", "intI1", .x)) %>%
  rename_with(~ gsub("EC_custom_intI2.*", "intI2", .x)) %>%
  rename_with(~ gsub("Shigella_flexneri_chloramphenicol_acetyltransferase", "catA1", .x, fixed = TRUE)) %>% 
  select(where(is.character), where( ~ is.integer(.x) && sum(.x) <= .9*nrow(geno_meta))) %>%
  select(sort(names(.))) %>%
  relocate(all_of(meta_cols), contains("intI"))

# Calculate total carriage of each gene in the collection
arg_totals <- t(args %>% summarise(across(where(is.integer), sum)))

arg_totals <- as_tibble(arg_totals, rownames = "Gene", .name_repair = "minimal") %>% 
  rename(Total = 2) %>% 
  mutate(Percentage = round(Total/nrow(geno_meta)*100, 2))

# Get hits from VFDB
vags <- geno_meta %>% 
  select(all_of(meta_cols), starts_with("vfdb")) %>%
  rename_with(~ gsub("vfdb_", "", .x))

# Select additional virulence genes from our custom database
custom_vags <- geno_meta %>% 
  select(Name, contains("EC_custom")) %>%
  rename_with(~ gsub("EC_custom_", "", .x, fixed = TRUE)) %>%
  mutate(usp = as.integer(rowSums(select(.,starts_with("usp"))))) %>%
  select(-starts_with("usp_")) %>%
  mutate(eitA = as.integer(rowSums(select(.,starts_with("eitA"))))) %>%
  select(-starts_with("eitA_")) %>% 
  select(Name,
         starts_with(c("cba", "cbi", "cjr", 
                       "cva", "cvi","eit",
                       "fecA", "hek", "hyl",
                       "iha","iss", "merA",
                       "ompT", "silA",
                       "terA", "traT", "usp"))) %>%
  rename_with(~ gsub("_[A-Z]{1,2}.*", "", .x)) %>%
  rename_with(~ gsub("_pAPEC-O1-ColBM", "", .x, fixed =TRUE)) %>%
  rename_with(~ gsub("_pUTI89", "", .x, fixed =TRUE)) %>%
  rename_with(~ gsub("_type3", "", .x, fixed =TRUE)) %>%
  rename_with(~ gsub("_|VFG1539", "", .x, fixed =TRUE)) %>%
  rename_with(~ gsub("_chromosomal", "_1", .x, fixed =TRUE)) %>%
  rename_with(~ gsub("_episomal", "_2", .x, fixed =TRUE))

# Join VFDB and custom hits
vags <- left_join(vags, custom_vags) %>% select(sort(names(.))) %>%
  relocate(all_of(meta_cols))

# Filter out genes we don't want due to rarity
# cps, fae, nle, yag
vags_preheat <- vags %>% select(-starts_with(c("cps", "fae", "nle", "yag")))

# Split out operons we want to filter to marker genes
# clbA, entB, escC, espL*, espX*, espY*, fepA, fimH, gspM, shuA, tssA, ybtA
# Select the markers
vag_operon_markers <- vags %>% 
  select(Name, clbA, entB, escC, espL1, espL4, espR1, espX1, espX4, espX5, starts_with("espY"), fepA, fimH, gspM, shuA, tssA, ybtA)

# Remove the operons to be replaced with markers
vags_preheat <- vags_preheat %>% 
  select(-starts_with(c("clb", "ent", "esc", "esp", "fep", "fim", "gsp", "shu", "tss", "ybt")))

# Rejoin the markers and filter out genes present in less than 10% of genomes
vags_preheat <- left_join(vags_preheat, vag_operon_markers) %>% 
  select(where(is.character), where( ~ is.integer(.x) && sum(.x) >= .10*nrow(vags))) %>% 
  select(sort(names(.))) %>%
  relocate(all_of(meta_cols))

# Show differences in VAGs when filtering out columns by total presence
  # vag_cols <- names(vags)
  # vag_cols5 <- names(vags %>% select(where(is.character), where( ~ is.integer(.x) && sum(.x) >= .05*nrow(vags))))
  # vag_cols10 <- names(vags %>% select(where(is.character), where( ~ is.integer(.x) && sum(.x) >= .1*nrow(vags))))
  # 
  # setdiff(vag_cols, vag_cols5)
  # setdiff(vag_cols5, vag_cols10)

# Calculate total carriage of each gene in the collection
vag_totals <- t(vags %>% summarise(across(where(is.integer), sum)))
vag_totals <- as_tibble(vag_totals, rownames = "Gene", .name_repair = "minimal") %>% 
  rename(Total = 2) %>% 
  mutate(Percentage = round(Total/nrow(geno_meta)*100, 2))

# Extract mobile genetic elements from the ISFinder database. NB: this data is not covered in the paper
mges <- geno_meta %>% 
  select(all_of(meta_cols), starts_with("ISfinder_")) %>%
  rename_with(~ gsub("ISfinder_Feb_2020_", "", .x)) %>%
  rename_with(~ gsub(":.*", "", .x))

# Extract plasmid related genes from the plasmidfinder database and tidy up messy names
plas <- geno_meta %>% 
  select(all_of(meta_cols), starts_with("plasmidfinder")) %>%
  rename_with(~ gsub("plasmidfinder_", "", .x)) %>%
  mutate(IncBOKZ = as.integer(rowSums(select(.,starts_with("IncB"))))) %>%
  select(-starts_with("IncB/")) %>%
  mutate(IncX1 = as.integer(rowSums(select(.,starts_with("IncX1"))))) %>%
  select(-starts_with("IncX1_")) %>%
  mutate(IncFII = as.integer(rowSums(select(.,starts_with("IncFII"))))) %>%
  select(-starts_with("IncFII("), -starts_with("IncFII_")) %>%
  mutate(IncFIA = as.integer(rowSums(select(.,starts_with("IncFIA"))))) %>%
  select(-starts_with("IncFIA("), -starts_with("IncFIA_")) %>%
  mutate(IncFIB = as.integer(rowSums(select(.,starts_with("IncFIB"))))) %>%
  select(-starts_with("IncFIB(")) %>%
  rename_with(~ gsub("_[0-9].*$", "", .x)) %>%
  rename_with(~ gsub("(", "_", .x, fixed = TRUE)) %>%
  rename_with(~ gsub(")", "", .x, fixed = TRUE)) %>%
  rename_with(~ gsub("/", "", .x, fixed = TRUE)) %>%
  rename_with(~ gsub("_FII", "", .x, fixed = TRUE)) %>%
  select(sort(names(.))) %>%
  relocate(all_of(meta_cols))

# Calculate total carriage of each gene in the collection
plas_totals <- t(plas %>% summarise(across(where(is.integer), sum)))
plas_totals <- as_tibble(plas_totals, rownames = "Gene") %>% 
  rename(Total = 2) %>% 
  mutate(Percentage = round(Total/nrow(geno_meta)*100, 2))

# Summary of IncF RSTs in the collection
F_RST_totals <- geno_meta %>% 
  filter(`F Plasmid` != "No F Plasmid") %>% 
  group_by(`F Plasmid`, ColV, pUTI89) %>%
  summarise(Total = n()) %>% 
  mutate(Percentage = round(Total/nrow(geno_meta)*100, 2))

# Put dataframes into a list so they can be converted to matrices for heatmap visualisation
gene_list <- list(args = args, vags = vags_preheat, mges = mges, plas = plas)

# Loop that converts the gene screening dataframes into binary heatmaps named `geneprefix_heat` e.g args heatmap is called args_heat
for (f in 1:length(gene_list)){
  heat <- gene_list[[f]]
  heat <- heat %>% select(where(is_integer)) %>% as.matrix()
  rownames(heat) <- gene_list[[f]]$Name
  heat[is.na(heat)] <- 0
  heat[heat >= 1] <- "Present"
  heat[heat == 0] <- "Absent"
  heat <- as.data.frame(cbind(heat, tree_meta %>% select(ST, ColV, pUTI89)))
  heat <- as.data.frame(cbind(heat, meta %>% select(Phylogroup)))
  heat <- heat %>% select(Phylogroup, ST, ColV, pUTI89, everything())
  assign(paste(names(gene_list[f]), "heat", sep="_"), heat)
}

# Add phenotypic data to args_heat
pheno <- meta %>% select(Amp.25.d:Ak.30.d) %>% 
  rename_with(~ gsub(".d", "", .x, fixed = TRUE))

pheno2 <- meta %>% select(Name, Amp.25.d:Ak.30.d) %>% 
  rename_with(~ gsub(".d", "", .x, fixed = TRUE))

pheno_sums <- pheno2 %>% mutate(`3GP` = rowSums(across(c(Amp.25, Amc.30), `%in%`, "R"))) %>% 
  mutate(`3GP` = case_when(`3GP` >= 1 ~ "R", `3GP` == 0 ~ "S")) %>% 
  select(Name, `3GP`, CL.100, Enr.5, Te.30, SF.300, W.5) %>%
  filter(CL.100 != "ND") %>%
  mutate(count = rowSums(across(c(`3GP`, CL.100, Enr.5, Te.30, SF.300, W.5),`%in%`, "R"))) %>% 
  mutate(MDR = case_when(count >= 3 ~ "Yes", count < 3 ~ "No"))

# Join MDR data to metadata
meta <- left_join(meta, pheno_sums %>% select(Name, MDR), by = "Name")

left_join(args2, args %>% select(Name, intI1), by = "Name") %>% 
  filter(ST == "ST1193"|ST == "ST38"|ST == "ST10"|ST == "ST117"|ST == "ST58") %>% 
  group_by(ST, intI1, ColV) %>% summarise(n())

pheno_names <- colnames(pheno)

pheno <- pheno %>% mutate(across(everything(), ~recode(., R = "Resistant", S = "Susceptible"))) %>%
  mutate(across(everything(), ~replace_na(., "ND")))

args_heat <- cbind(args_heat, pheno) %>%
  relocate(all_of(c("Phylogroup", "ST", "ColV", "pUTI89")), all_of(pheno_names), everything())

# Additional colours for visualisation
heat_clrs <- c("Present" = "#8dd3c7", "Absent" = "#ededed", "Yes" = "#df03fc","No" = "white")

####-----GENE TOTALS-----####
# Total ARGs
args2 <- args %>% select(-intI1, -intI2) %>% 
  mutate(`Total ARGs` = rowSums(across(where(is.integer)))) %>% 
  select(Name, Phylogroup, ST, ColV, pUTI89, `Total ARGs`)

# Average ARGs by ST
args2 %>% group_by(ST) %>% summarise(`Average ARGs` = mean(`Total ARGs`)) %>% arrange(desc(`Average ARGs`))

# Average ARGs by Phylogroup
args2 %>% group_by(Phylogroup) %>% summarise(`Average ARGs` = mean(`Total ARGs`)) %>% arrange(desc(`Average ARGs`))
args2 %>% group_by(Phylogroup) %>% summarise(`Max ARGs` = max(`Total ARGs`)) %>% arrange(desc(`Max ARGs`))
args2 %>% group_by(Phylogroup) %>% summarise(`Min ARGs` = min(`Total ARGs`)) %>% arrange(desc(`Min ARGs`))

# Average ARGs by ColV
args2 %>% group_by(ColV) %>% summarise(`Average ARGs` = mean(`Total ARGs`)) %>% arrange(desc(`Average ARGs`))

# Average ARGs by pUTI89
args2 %>% group_by(pUTI89) %>% summarise(`Average ARGs` = mean(`Total ARGs`)) %>% arrange(desc(`Average ARGs`))

# Total VAGs
vags2 <- vags %>% mutate(`Total VAGs` = rowSums(across(where(is.integer)))) %>% 
  select(Name, Phylogroup, ST, ColV, pUTI89, `Total VAGs`)

# Average VAGs by ST
vags2 %>% group_by(ST) %>% summarise(`Average VAGs` = mean(`Total VAGs`)) %>% arrange(desc(`Average VAGs`))

# Average VAGs by Phylogroup
vags2 %>% group_by(Phylogroup) %>% summarise(`Average VAGs` = mean(`Total VAGs`)) %>% arrange(desc(`Average VAGs`))

# Average VAGs by ColV
vags2 %>% group_by(ColV) %>% summarise(`Average VAGs` = mean(`Total VAGs`)) %>% arrange(desc(`Average VAGs`))

# Average VAGs by pUTI89
vags2 %>% group_by(pUTI89) %>% summarise(`Average VAGs` = mean(`Total VAGs`)) %>% arrange(desc(`Average VAGs`))

####-----INTI1 TRUNCATIONS-----####
intI1_trunc <- read_tsv("data/pipelord_results/MVC.intI1.blast.txt") %>% mutate(COV = round(NBASE/GENELEN*100, 2))

intI1_trunc <- left_join(intI1_trunc, full_mlst) %>% mutate(ST = gsub("^", "ST", ST))

intI1_trunc <- left_join(intI1_trunc, full_F)

intI1_trunc <- intI1_trunc %>% filter(COV <= 90, Name %in% meta$Name) %>% select(Name, ST, `F Plasmid`, intI1_length = NBASE) 

intI1_trunc_summary <- left_join(intI1_trunc %>% 
                                   group_by(ST, `F Plasmid`, intI1_length) %>% 
                                   summarise(Count_of_truncation = n()), 
                                 full_mlst %>% 
                                   mutate(ST = gsub("^", "ST", ST)) %>% 
                                   group_by(ST) %>% 
                                   summarise(Count_of_ST = n()), by = "ST") %>%
                        mutate(Percent_of_ST = round(Count_of_truncation/Count_of_ST*100,2))

####-----FIGURE 1 - METADATA SUMMARY GRAPHS-----####
fig1a <- ggplot(meta, aes(x = as.character(year))) +
  geom_bar(aes(fill = Specimen)) +
  scale_fill_manual(values = source_clrs) +
  scale_x_discrete(name = "Year")+
  scale_y_continuous(name = "Count", expand = c(0, 0), limits = c(0, 60), n.breaks = 6) +
  theme_classic()+
  theme(axis.text.x = element_text(color = "grey20", 
                                   size = 8,
                                   vjust = ,
                                   face = "plain"),
                                   legend.position = "right",
                                   legend.key.size = unit(3.5, "mm"),
                                   legend.title = element_text(size = 10),
                                   legend.text = element_text(size = 8),
                                   legend.direction = "vertical") 

fig1b <- ggplot(meta, aes(x = Phylogroup)) +
  geom_bar(aes(fill = ST)) +
  scale_fill_manual(values = ST_clrs) +
  scale_x_discrete(name = "Phylogroup")+
  scale_y_continuous(name = "Count", expand = c(0, 0), limits = c(0, 250), n.breaks = 6) +
  theme_classic()+
  theme(axis.text.x = element_text(color = "grey20", 
                                   size = 8, 
                                   vjust = , 
                                   face = "plain"), 
        legend.position = "right",
        legend.key.size = unit(3.5, "mm"),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8),
        legend.direction = "vertical") 
  

fig1c <- ggplot(meta, aes(x = Specimen)) +
  geom_bar(aes(fill = ST)) +
  scale_fill_manual(values = ST_clrs) +
  scale_x_discrete(name = "Specimen")+
  scale_y_continuous(name = "Count", expand = c(0, 0), limits = c(0, 220), n.breaks = 8) +
  theme_classic()+
  theme(axis.text.x = element_text(color = "grey20", 
                                   size = 7, 
                                   vjust = 1,
                                   hjust = 0.5,
                                   face = "plain",
                                   angle = ), 
        legend.position = "right",
        legend.key.size = unit(3.5, "mm"),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8),
        legend.direction = "vertical")

mt <- meta %>% mutate(Count=rep(1,nrow(meta)))
labs <- meta %>% group_by(year) %>% tally %>% mutate(ST=NA)

fig1d <- ggplot(mt, aes(x = as.character(year), y = Count)) +
  geom_bar(aes(fill = ST), stat = "identity", position = "fill") +
  geom_text(data=labs,aes(label=paste0("n=",n), y = rep(0.97,11))) +
  scale_fill_manual(values = ST_clrs) +
  scale_x_discrete(name = "Year")+
  scale_y_continuous(name = "Proportion", expand = c(0, 0), limits = c(0, 1.0), n.breaks = 5) +
  theme_classic()+
  theme(axis.text.x = element_text(color = "grey20", 
                                   size = 8,
                                   vjust = ,
                                   face = "plain"),
        legend.position = "right",
        legend.key.size = unit(3.5, "mm"),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8),
        legend.direction = "vertical")

fig1ab <- ggarrange(plotlist= list(fig1a, fig1b),
                    ncol =2, labels =c("a)", "b)"),
                    legend = "right")

fig1cd <- ggarrange(plotlist= list(fig1c, fig1d),
                    ncol =2, labels =c("c)", "d)"),
                    legend = "right", common.legend = TRUE, align = "v")

fig1 <- ggarrange(plotlist = list(fig1ab, fig1cd),
                  nrow =2)

ggsave("fig1.metadata.pdf", 
       fig1, 
       path = "outputs/figures", 
       device = "pdf", 
       width = 297,
       height = 210,
       unit ="mm")

####------FIGURE 2 - TREE + METADATA-----####
## Clade colours
a_node<- MRCA(mproot_core, c("MVC62", "MVC126"))
b1_node<- MRCA(mproot_core, c("MVC184", "MVC778"))
b2_node<- MRCA(mproot_core, c("MVC423", "MVC134"))
c_node<- MRCA(mproot_core, c("MVC561", "MVC652"))
d_node<- MRCA(mproot_core, c("MVC208", "MVC763"))
e_node<- MRCA(mproot_core, c("MVC73", "MVC307"))
f_node<- MRCA(mproot_core, c("MVC633", "MVC532"))
g_node<- MRCA(mproot_core, c("MVC509", "MVC562"))
cryptic_node1 <- MRCA(mproot_core, c("MVC382", "MVC526"))
cryptic_node2 <- MRCA(mproot_core, c("MVC370", "MVC421"))

## Draw tree
core_t_rotate_up <- rotate_tree(ggtree(mproot_core, layout = "fan", open.angle = 8.5) %<+%
                                  phylo_tips +
                                  geom_tippoint(size = .6, aes(color = Phylogroup)) +
                                  scale_color_manual(name = "Phylogroup", values = phylo_clrs) +
                                  theme(legend.key.size = unit(3, "mm")) +
                                  geom_cladelabel(node = a_node, label = "A", offset =.11, offset.text = .01) +
                                  geom_cladelabel(node = b1_node, label = "B1", offset =.111, offset.text = .015)+
                                  geom_cladelabel(node = b2_node, label = "B2", offset =.105, offset.text = .005)+
                                  geom_cladelabel(node = c_node, label = "C", offset =.1135, offset.text = .01)+
                                  geom_cladelabel(node = d_node, label = "D", offset =.11, offset.text = .01)+
                                  geom_cladelabel(node = e_node, label = "E", offset =.116, hjust = 0.5, offset.text = .01)+
                                  geom_cladelabel(node = f_node, label = "F", offset =.114, offset.text = .01)+
                                  geom_cladelabel(node = g_node, label = "G", offset =.114, offset.text = .012) +
                                  geom_cladelabel(node = cryptic_node1, label = "Cryptic 1", offset =.043, 
                                                  offset.text = .008, angle = 50, fontsize = 3.5) +
                                  geom_cladelabel(node = cryptic_node2, label = "Cryptic 2", offset =.043, 
                                                  offset.text = .008, angle = 50, fontsize = 3.5),
                                  91)

# Draw tree with metadata
fig2 <- gheatmap(p = core_t_rotate_up,
         data = tree_meta,
         colnames = TRUE,
         colnames_angle = ,
         colnames_offset_y = ,
         hjust = 0,
         font.size = 1.7,
         width = 0.25,
         color = rgb(0, 0, 0, alpha = .01)
) + 
  theme(legend.position = "right",
        legend.box = "vertical",
        legend.key.size = unit(3, "mm"),
        #legend.title = element_blank(),
        legend.text = element_text(size = 8))+
  scale_fill_manual(
    name = "Typing",
    aesthetics = c("fill"),
    values = tree_vars,
    na.value = 'white')

ggsave("fig2.tree.pdf", 
       fig2, 
       path = "outputs/figures", 
       device = "pdf", 
       width = 297,
       height = 210,
       unit ="mm")

# Extra tree with senB marker gene
plus_tree <- left_join(meta %>% select(Name, Specimen, ST, Serotype = OH_type, `F Plasmid`, ColV, pUTI89), vags %>% select(Name, senB, irp2, fyuA))

plus_tree <- plus_tree %>% mutate(HPI = case_when(fyuA + irp2 == "0" ~ "HPI-Neg", fyuA + irp2 == "2" ~ "HPI-Pos"),
         senB = case_when(senB == "0" ~ "senB-Neg", senB == "1" ~ "senB-Pos"),
         ColV = case_when(ColV == "Yes" ~ "ColV-Pos", ColV == "No" ~ "ColV-Neg"),
         pUTI89 = case_when(pUTI89 == "Yes" ~ "pUTI89-Pos", pUTI89 == "No" ~ "pUTI89-Neg")) %>% select(-irp2, -fyuA, -Name)

senB_clr <- c("senB-Neg"="#ffffff", "senB-Pos"="#E7298A")
HPI_clr <- c("HPI-Neg"="#ffffff", "HPI-Pos"="#38bcd6")

rownames(plus_tree) <- nems

tree_plus <- gheatmap(p = core_t_rotate_up,
                      data = plus_tree,
                      colnames = TRUE,
                      hjust = 0,
                      font.size = 1.7,
                      width = 0.25,
                      color = rgb(0, 0, 0, alpha = .2)
) + 
  theme(legend.position = "right",
        legend.box = "vertical",
        legend.key.size = unit(3, "mm"),
        #legend.title = element_blank(),
        legend.text = element_text(size = 8))+
  scale_fill_manual(
    name = "Typing",
    aesthetics = c("fill"),
    values = c(tree_vars, senB_clr, HPI_clr),
    na.value = 'white')

####----- FIGURE 3 - ColV and pUTI89-----####
fig3a <- ggplot(meta, aes(Phylogroup)) + geom_bar(aes(fill = ColV)) +
  scale_fill_manual(values = c("Yes" = "#00BFC4FF", "No" = "#F8766DFF")) +
  scale_y_continuous(name = "Count", expand = c(0, 0), limits = c(0,250), n.breaks = 6) +
  scale_x_discrete(name = "Phylogroup")+
  theme_classic()+
  theme(axis.text.x = element_text(color = "grey20", size = 10, vjust = , hjust = .5, face = "plain", angle = ),
        axis.text.y = element_text(color = "grey20", size = 10, vjust = ,face = "plain"),  
        axis.title.x = element_text(color = "grey20",vjust = , size = 12, face = "plain"),
        axis.title.y = element_text(color = "grey20",vjust = , size = 12, face = "plain"),
        legend.position = "right",
        legend.box = "vertical",
        legend.key.size = unit(5, "mm"),
        legend.title = element_text(size=12),
        legend.text = element_text(size = 10))

fig3b <- ggplot(meta, aes(x = fct_infreq(factor(ST)))) +
  geom_bar(aes(fill = ColV)) +
  scale_fill_manual(values = c("Yes" = "#00BFC4FF", "No" = "#F8766DFF")) +
  scale_x_discrete(name = "ST")+
  scale_y_continuous(name = "Count", expand = c(0, 0), limits = c(0, 125), n.breaks = 5) +
  theme_classic()+
  theme(axis.text.x = element_text(color = "grey20", size = 8, angle = 30, vjust = 1, hjust = .9, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 10, vjust = ,face = "plain"),  
        axis.title.x = element_text(color = "grey20", vjust = , size = 14, face = "plain"),
        axis.title.y = element_text(color = "grey20", vjust = , size = 14, face = "plain"), 
        legend.position = "right", 
        legend.key.size = unit(5, "mm"),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.direction = "vertical")

fig3c <- ggplot(meta, aes(Phylogroup)) + geom_bar(aes(fill = pUTI89)) +
  scale_fill_manual(values = c("Yes" = "#97fc9f", "No" = "#ffadd5")) +
  scale_y_continuous(name = "Count", expand = c(0, 0), limits = c(0,250), n.breaks = 6) +
  scale_x_discrete(name = "Phylogroup")+
  theme_classic()+
  theme(axis.text.x = element_text(color = "grey20", size = 10, vjust = , hjust = .5, face = "plain", angle = ),
        axis.text.y = element_text(color = "grey20", size = 10, vjust = ,face = "plain"),  
        axis.title.x = element_text(color = "grey20",vjust = , size = 12, face = "plain"),
        axis.title.y = element_text(color = "grey20",vjust = , size = 12, face = "plain"),
        legend.position = "right",
        legend.box = "vertical",
        legend.key.size = unit(5, "mm"),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10))

fig3d <- ggplot(meta, aes(x = fct_infreq(factor(ST)))) +
  geom_bar(aes(fill = pUTI89)) +
  scale_fill_manual(values = c("Yes" = "#97fc9f", "No" = "#ffadd5")) +
  scale_x_discrete(name = "ST")+
  scale_y_continuous(name = "Count", expand = c(0, 0), limits = c(0, 125), n.breaks = 5) +
  theme_classic()+
  theme(axis.text.x = element_text(color = "grey20", size = 8, angle = 30, vjust = 1, hjust = .9, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 10, vjust = ,face = "plain"),  
        axis.title.x = element_text(color = "grey20", vjust = , size = 14, face = "plain"),
        axis.title.y = element_text(color = "grey20", vjust = , size = 14, face = "plain"),
        legend.position = "right", 
        legend.key.size = unit(5, "mm"),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.direction = "vertical") 

f3ab <- ggarrange(plotlist = list(fig3a, fig3b),
          ncol = 2, nrow =1,
          labels =c("a)","b)"), legend = "right", align = "h", common.legend = TRUE)

f3cd <- ggarrange(plotlist = list(fig3c, fig3d),
                  ncol = 2, nrow =1,
                  labels =c("c)","d)"), legend = "right", align = "h", common.legend = TRUE)

fig3 <- ggarrange(f3ab, f3cd, nrow =2)

ggsave("fig3.colVpUTI.pdf", 
       fig3, 
       path = "outputs/figures", 
       device = "pdf", 
       width = 297,
       height = 210,
       unit ="mm")

####----- FIGURE 4 - ARG/VAG PLOTS -----####
fig4a <- ggboxplot(args2, x = "Phylogroup", y = "Total ARGs",
                   color = "Phylogroup",
                   order = c(sort(unique(args2$Phylogroup))),
                   ylab = "ARGs", xlab = "Phylogroup") +
  scale_y_continuous(name = "Total ARGs", expand = c(0, 0), limits = c(0,20), n.breaks = 5) +
  scale_color_manual(values = phylo_clrs) +
  theme_classic() +
  theme(axis.text.x = element_text(color = "grey20", size = 10, vjust = , hjust = .5, face = "plain", angle = ),
        axis.text.y = element_text(color = "grey20", size = 10, vjust = ,face = "plain"),  
        axis.title.x = element_text(color = "grey20",vjust = , size = 12, face = "plain"),
        axis.title.y = element_text(color = "grey20",vjust = , size = 12, face = "plain"),
        legend.position = "none",
        legend.box = "vertical",
        legend.key.size = unit(5, "mm"),
        legend.title = element_text(size= 12),
        legend.text = element_text(size = 10))

fig4b <- ggboxplot(args2 %>% filter(ST != "Other"), x = "ST", y = "Total ARGs",
                   color = "ST", 
                   order = c(levels(fct_infreq(factor(args2$ST)))),
                   ylab = "", xlab = "ST") +
  scale_y_continuous(expand = c(0, 0), limits = c(0,20), n.breaks = 5) +
  scale_color_manual(values = ST_clrs[2:length(ST_clrs)]) +
  theme_classic() +
  theme(axis.text.x = element_text(color = "grey20", size = 10, vjust = , hjust = .9, face = "plain", angle = 30),
        axis.text.y = element_text(color = "grey20", size = 10, vjust = ,face = "plain"),  
        axis.title.x = element_text(color = "grey20",vjust = , size = 12, face = "plain"),
        axis.title.y = element_text(color = "grey20",vjust = , size = 12, face = "plain"),
        legend.position = "none",
        legend.box = "vertical",
        legend.key.size = unit(5, "mm"),
        legend.title = element_text(size=12),
        legend.text = element_text(size = 10))

fig4c <- ggboxplot(vags2, x = "Phylogroup", y = "Total VAGs",
                   color = "Phylogroup",
                   order = c(sort(unique(args2$Phylogroup))),
                   ylab = "VAGs", xlab = "Phylogroup") +
  scale_y_continuous(name = "Total VAGs", expand = c(0, 0), limits = c(0,150), n.breaks = 7) +
  scale_color_manual(values = phylo_clrs) +
  theme_classic() +
  theme(axis.text.x = element_text(color = "grey20", size = 10, vjust = , hjust = .5, face = "plain", angle = ),
        axis.text.y = element_text(color = "grey20", size = 10, vjust = ,face = "plain"),  
        axis.title.x = element_text(color = "grey20",vjust = , size = 12, face = "plain"),
        axis.title.y = element_text(color = "grey20",vjust = , size = 12, face = "plain"),
        legend.position = "none",
        legend.box = "vertical",
        legend.key.size = unit(5, "mm"),
        legend.title = element_text(size= 12),
        legend.text = element_text(size = 10))

fig4d <- ggboxplot(vags2 %>% filter(ST != "Other"), x = "ST", y = "Total VAGs",
                   color = "ST", 
                   order = c(levels(fct_infreq(factor(args2$ST)))),
                   ylab = "", xlab = "ST") +
  scale_y_continuous(expand = c(0, 0), limits = c(0,150), n.breaks = 7) +
  scale_color_manual(values = ST_clrs[2:length(ST_clrs)]) +
  theme_classic() +
  theme(axis.text.x = element_text(color = "grey20", size = 10, vjust = , hjust = .9, face = "plain", angle = 30),
        axis.text.y = element_text(color = "grey20", size = 10, vjust = ,face = "plain"),  
        axis.title.x = element_text(color = "grey20",vjust = , size = 12, face = "plain"),
        axis.title.y = element_text(color = "grey20",vjust = , size = 12, face = "plain"),
        legend.position = "none",
        legend.box = "vertical",
        legend.key.size = unit(5, "mm"),
        legend.title = element_text(size=12),
        legend.text = element_text(size = 10))

fig4 <- ggarrange(plotlist = list(fig4a, fig4b, fig4c, fig4d),
                  ncol = 2, nrow =2,
                  labels =c("a)","b)", "c)","d)"), legend = "none", align = "h")

ggsave("fig4.argvag.pdf", 
       fig4, 
       path = "outputs/figures", 
       device = "pdf", 
       width = 297,
       height = 210,
       unit ="mm")

####------FIGURES S1-S3 - GENE HEATMAPS-----####
# Define the tree
pheno.tree.heatmap <- ggtree(mproot_core, branch.length = "none") %<+%
  phylo_tips +
  geom_tippoint(size = .6, aes(color = Phylogroup))

tree.heatmap <- ggtree(mproot_core, branch.length = "none") 

pheno_clrs <- c("Resistant" = "#000000", "Susceptible" = "#FFFFFF", "ND"= "grey")

heat_meta_clrs <- c(tree_vars, phylo_clrs)

## AMR HEATMAP ##
figS1 <- gheatmap(tree.heatmap,
                  offset = -5,
                  args_heat, 
                  width = 30,
                  font.size = 2,
                  colnames_offset_x = 0.00001,
                  colnames_offset_y = 2.2,
                  colnames_position = "top",
                  colnames_angle = 45,
                  hjust = 0,
                  color = NULL) +
  scale_fill_manual(name = "Data", values = c(pheno_clrs, phylo_clrs, ST_clrs, colv_clrs, puti89_clrs, heat_clrs)) +
  ylim(c(0,400)) +
  theme(legend.position = "none",
        legend.key.size = unit(3, "mm"),
        legend.title = element_text(),
        legend.text = element_text(size = 6))

ggsave("figS1.args.pdf", 
       figS1, 
       path = "outputs/figures", 
       device = "pdf", 
       width = 297,
       height = 210,
       unit ="mm")

## VIRULENCE HEATMAP ##
figS2 <- gheatmap(tree.heatmap,
                  offset = -3,
                  vags_heat, 
                  width = 30,
                  font.size = 2,
                  colnames_offset_x = 0.00001,
                  colnames_offset_y = 2.2,
                  colnames_position = "top",
                  colnames_angle = 60,
                  hjust = 0,
                  color = NULL) +
  scale_fill_manual(name = "Data", values = c(heat_clrs, heat_meta_clrs)) +
  ylim(c(0,400)) +
  theme(legend.position = "none",
        legend.key.size = unit(3, "mm"),
        legend.title = element_text(),
        legend.text = element_text(size = 6))

ggsave("figS2.vags.pdf", 
       figS2, 
       path = "outputs/figures", 
       device = "pdf", 
       width = 297,
       height = 210,
       unit ="mm")

## PLASMID REPLICON HEATMAP ##
figS3 <- gheatmap(tree.heatmap, 
                  offset = -13,
                  plas_heat, 
                  width = 30,
                  font.size = 2.5,
                  colnames_offset_x = 0.00001,
                  colnames_offset_y = 2.2,
                  colnames_position = "top",
                  colnames_angle = 45,
                  hjust = 0,
                  color = NULL) +
  scale_fill_manual(name = "Data", values = c(heat_clrs, heat_meta_clrs)) +
  ylim(c(0,400)) +
  theme(legend.position = "none",
        legend.key.size = unit(3, "mm"),
        legend.title = element_text(),
        legend.text = element_text(size = 6))

ggsave("figS3.plas.pdf", 
       figS3, 
       path = "outputs/figures", 
       device = "pdf", 
       width = 297,
       height = 210,
       unit ="mm")

####------FIGURES S4, S5 PLASMID MAPS-----####
## Fig S4
# Use the tree.heatmap tree with phylogroup tips
# Refine meta data for use in this figure
puti_meta <- meta %>% select(Phylogroup, ST, pUTI89)

# Recode pUTI89 values for clarity on the figure
puti_meta$pUTI89 <- puti_meta$pUTI89 %>% recode(Yes = "pUTI89-Pos", No = "pUTI89-Neg")

# Add rownames
rownames(puti_meta) <- nems

# Visualise the tree with metadata
figS4a <- gheatmap(tree.heatmap, 
           puti_meta,
           offset = ,
           width = .05,
           font.size = 2,
           colnames= TRUE,
           colnames_position = "top",
           colnames_angle = 70,
           colnames_offset_y = 5,
           hjust = 0,
           color = rgb(0, 0, 0, alpha = 0)) +
  scale_fill_manual(name = "Data", values = c(phylo_clrs, ST_clrs, "pUTI89-Pos" = "#df03fc","pUTI89-Neg" = "white"))+
  ylim(c(0,400))+
  theme(legend.key.size = unit(3, "mm"),
        legend.title = element_text(),
        legend.text = element_text(size = 6))

# Visualise the tree with hits to pUTI89
figS4b <- gheatmap(tree.heatmap, 
         pUTI89_binned_hits,
         offset = 35,
         width = 30,
         font.size = 2,
         colnames= FALSE,
         hjust = 0,
         color = rgb(0, 0, 0, alpha = 0)) +
  scale_fill_gradient(low = "white", high = "#c6ffb8", na.value = "white")+
  theme(legend.key.size = unit(3, "mm"),
        legend.title = element_text(),
        legend.text = element_text(size = 6))

ggsave("figS4a.puti89.meta.pdf", 
       figS4a, 
       path = "outputs/figures", 
       device = "pdf", 
       width = 297,
       height = 210,
       unit ="mm")       

ggsave("figS4b.puti89.heat.pdf", 
       figS4b, 
       path = "outputs/figures", 
       device = "pdf", 
       width = 297,
       height = 210,
       unit ="mm")

## Fig S5
# Refine meta data for use in this figure
pCERC4_meta <- meta %>% select(Phylogroup, ST, ColV)

# Recode ColV values for clarity on the figure
pCERC4_meta$ColV <- pCERC4_meta$ColV %>% recode(Yes = "ColV-Pos", No = "ColV-Neg")

# Add rownames
rownames(pCERC4_meta) <- nems

# Visualise the tree with metadata
figS5a <- gheatmap(tree.heatmap, 
         pCERC4_meta,
         offset = ,
         width = .05,
         font.size = 2,
         colnames= FALSE,
         colnames_position = "top",
         colnames_angle = 70,
         colnames_offset_y = 5,
         hjust = 0,
         color = rgb(0, 0, 0, alpha = 0)) +
  scale_fill_manual(name = "Data", values = c(phylo_clrs, ST_clrs, colv_clrs))+
  ylim(c(0,400))+
  theme(legend.key.size = unit(3, "mm"),
        legend.title = element_text(),
        legend.text = element_text(size = 6))

# Visualise the tree with hits to pUTI89
figS5b <- gheatmap(tree.heatmap, 
         pCERC4_binned_hits,
         offset = 35,
         width = 30,
         font.size = 2,
         colnames= FALSE,
         hjust = 0,
         color = rgb(0, 0, 0, alpha = 0)) +
  scale_fill_gradient(low = "white", high = "#8dd3c7", na.value = "white")+
  theme(legend.key.size = unit(3, "mm"),
        legend.title = element_text(),
        legend.text = element_text(size = 6))

ggsave("figS5a.pCERC4.meta.pdf", 
       figS5a, 
       path = "outputs/figures", 
       device = "pdf", 
       width = 297,
       height = 210,
       unit ="mm",
       dpi=300) 

ggsave("figS5b.pCERC4.heat.pdf", 
       figS5b, 
       path = "outputs/figures", 
       device = "pdf", 
       width = 297,
       height = 210,
       unit ="mm",
       dpi=300)

####-----TABLE S1-----####
# Generate Table S1 by replacing simplified ST, OH and F designations with full data
tableS1 <- geno_meta %>% select(Name, where(is.integer))
tableS1 <- left_join(tableS1, full_mlst, by ="Name")
tableS1 <- left_join(tableS1, full_OH, by ="Name")
tableS1 <- left_join(tableS1, full_F, by ="Name")
tableS1 <- left_join(tableS1, meta %>% select(Name, Year = year, Specimen, Detail = Info.Nota,
                              Phylogroup, ColV, pCERC4_ID, pUTI89, pUTI89_ID, 
                              fimH_allele = fimH, FQR, MDR, IS26_signature)) %>%
  select(Name, Year, Specimen, Detail,
         Phylogroup, ST, OH_type, `F Plasmid`, ColV, pCERC4_ID, pUTI89, pUTI89_ID, 
         fimH_allele, FQR, MDR, IS26_signature, where(is.integer))


# Write Table S1
write_csv(tableS1,
          "outputs/data/TableS1.csv")

####-----VIRULENCE GENE CULL-----####
# The gene heatmap is huge so I would like to reduce it, particularly in instances where operon
# presence can be approximated by a marker gene. Also to remove genes that aren't really 
# considered bonafide VFs
# gheatmap(tree.heatmap,
#          vags_heat %>% select(starts_with("insert gene here")), 
#          width = 30,
#          font.size = 7,
#          colnames_offset_x = 0.00001,
#          colnames_offset_y = 2.2,
#          colnames_position = "top",
#          colnames_angle = 45,
#          hjust = 0,
#          color = rgb(0, 0, 0, alpha = .2)) +
#   scale_fill_manual(name = "Data", values = c(heat_clrs, heat_meta_clrs)) +
#   ylim(c(0,400)) +
#   theme(legend.key.size = unit(3, "mm"),
#         legend.title = element_text(),
#         legend.text = element_text(size = 6))
# 
# vags %>% select(starts_with("clb")) %>% colSums ## Colibactin: Consistent carriage across the operon, use clbA as marker
# vags %>% select(starts_with("cps")) %>% colSums ## E. faecalis capsular polysaccharide: only in one isolate: Exclude?
# vags %>% select(starts_with("ent")) %>% colSums ## Enterobactin: Fairly consistent carriage: use entB as marker
# vags %>% select(starts_with("esc")) %>% colSums ## T3SS: consistent carriage: use escC
# vags %>% select(starts_with("esp")) %>% colSums ## T3SS: some very rare, espL1 and X1, X4 and X5 in about 1/3 though, espL4 in 50: Keep espL, espX and espY
# vags %>% select(starts_with("fae")) %>% colSums ## Fimbriae: rare so exclude
# vags %>% select(starts_with("fep")) %>% colSums ## Ferrienterobactin: almost same for all use fepA
# vags %>% select(starts_with("fim")) %>% colSums ## Consistent, use fimH
# vags %>% select(starts_with("gsp")) %>% colSums ## T2SS: Fairly consistent, use gspM. Missing from most ST372
# vags %>% select(starts_with("kps")) %>% colSums ## Variable and not too many genes, keep all
# vags %>% select(starts_with("nle")) %>% colSums ## T3SS: rare, exclude
# vags %>% select(starts_with("pap")) %>% colSums ## Quite variable so keep all
# vags %>% select(starts_with("sfa")) %>% colSums ## S-fimbriae: highly variable - keep all for now
# vags %>% select(starts_with("shu")) %>% colSums ## Shigella-associated outer membrane heme receptor: Just use shuA
# vags %>% select(starts_with("tss")) %>% colSums ## T6SS: consistent, use tssA 
# vags %>% select(starts_with("yag")) %>% colSums ## E. coli common pilus, near ubiquitous so probably not 'virulence gene' per sé: exclude
# vags %>% select(starts_with("ybt")) %>% colSums ## Yersinabactin: consistent and associated with phylogroup B2, use ybtA

# Genes to select out of operons
# clbA, entB, escC, espL*, espX*, espY*, fepA, fimH, gspM, shuA, tssA, ybtA

# Genes to exclude
# cps, fae, nle, yag
