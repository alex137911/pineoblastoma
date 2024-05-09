# Human-Zebrafish Orthologs

# Remove objects in workspace
rm(list = ls())

# Required packages
suppressMessages(library(readr))
suppressMessages(library(dplyr))

# -------------------------------------------------------------------
# LOAD DATA
inpath <- "C:/Users/acale/OneDrive/Documents/Waterloo BME/Co-op/4C/pineoblastoma/Data/Input/"
inDir  <- sprintf("%s/Input", dirname(inpath))
setwd(inDir)

# Import HCOP Data (HGNC Comparison of Orthology Predictions)
# Downloaded May 8, 2024 from: https://genenames.org/tools/hcop/
# Methodology described in: doi: 10.1016/j.devcel.2023.11.001
hgnc_orthologs <- read_delim("human_zebrafish_hcop_fifteen_column.txt.gz", 
                             delim = "\t", escape_double = FALSE, trim_ws = TRUE)

hgnc_orthologsDF <- as.data.frame(hgnc_orthologs)
message(sprintf("Loaded %i records", nrow(hgnc_orthologsDF)))            # 132 678 records

# Import ensembl BioMart Data
# Exported May 8, 2024 from: https://www.ensembl.org/biomart/martview 
# Methodology described in: https://useast.ensembl.org/info/data/biomart/how_to_use_biomart.html
ensembl_orthologs <- read_delim("human_zebrafish_biomart_export.txt", 
                                delim = "\t", escape_double = FALSE, trim_ws = TRUE)

ensembl_orthologsDF <- as.data.frame(ensembl_orthologs)
message(sprintf("Loaded %i records", nrow(ensembl_orthologsDF)))         # 72 462 records

# Import ZFIN Data
# Downloaded May 8, 2024 from: https://zfin.org/downloads
header_column <- c("ZFIN_ID", "ZFIN_Symbol", "ZFIN_Name", "Human_Symbol", 
                   "Human_Name", "OMIM_ID", "Gene_ID", "HGNC_ID", "Evidence", "Pub_ID")

ZFIN_orthologs <- read_delim("human_zebrafish_zfin.txt",
                             delim = "\t", col_names = header_column, escape_double = FALSE, trim_ws = TRUE)

ZFIN_orthologsDF <- as.data.frame(ZFIN_orthologs)
message(sprintf("Loaded %i records", nrow(ZFIN_orthologsDF)))            # 43 612 records 

# -------------------------------------------------------------------
# PRE-PROCESSING
# Quality Check
no_entries <- sum(hgnc_orthologsDF$zebrafish_entrez_gene == "-")  # 11 645 blank
no_entries <- sum(hgnc_orthologsDF$zebrafish_ensembl_gene == "-") # 2 081 blank
no_entries <- sum(hgnc_orthologsDF$zebrafish_name == "-")         # 5 253 blank
no_entries <- sum(hgnc_orthologsDF$zfin_id == "-")                # 16 005 blank

# Drop rows with no zebrafish_ensembl_gene 
hgnc_orthologsDF <- hgnc_orthologsDF %>% filter(zebrafish_ensembl_gene != "-")

# HGNC Data
# Count occurrences of each human gene
# HGNChuman_counts <- hgnc_orthologsDF %>%
#   group_by(human_ensembl_gene) %>%
#   summarize(zebrafish_count = n_distinct(zebrafish_ensembl_gene), .groups = 'drop')

HGNChuman_counts <- hgnc_orthologsDF %>%
  group_by(human_ensembl_gene) %>%
  summarize(zebrafish_count = n(), .groups = 'drop')

# Count occurrences of each zebrafish gene
HGNCzebrafish_counts <- hgnc_orthologsDF %>%
  group_by(zebrafish_ensembl_gene) %>%
  summarize(human_count = n(), .groups = 'drop')

# Classify human to zebrafish relationships
HGNChuman_zebrafish <- HGNChuman_counts %>%
  mutate(relation_type = case_when(
    zebrafish_count == 1 ~ "One to One",
    zebrafish_count > 1 ~ "One to Many",
    TRUE ~ "No Match"
  ))

# Classify zebrafish to human relationships
HGNCzebrafish_human <- HGNCzebrafish_counts %>%
  mutate(relation_type = case_when(
    human_count == 1 ~ "One to One",
    human_count > 1 ~ "Many to One",
    TRUE ~ "No Match"
  ))

# Output the summary of relationships
print(human_to_zebrafish)
print(zebrafish_to_human)

# Ensembl data --> go both directions and drop rows with "NA" for human gene stable ID column





