suppressMessages(library(rtracklayer))
suppressMessages(library(readr))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(ggplot2))


inpath <- "C:/Users/acale/OneDrive/Documents/Waterloo BME/Co-op/4C/pineoblastoma/Data/Input/"
inDir  <- sprintf("%s/Input", dirname(inpath))
setwd(inDir)

# Load and Clean Data ---------------------------------------------------------
# Import ZFIN Data
# Downloaded May 8, 2024 from: https://zfin.org/downloads
header_column <- c("ZFIN_ID", "ZFIN_Symbol", "ZFIN_Name", "Human_Symbol", 
                   "Human_Name", "OMIM_ID", "Gene_ID", "HGNC_ID", "Evidence", "Pub_ID")

ZFIN_orthologs <- read_delim("human_zebrafish_zfin.txt",
                             delim = "\t", col_names = header_column, escape_double = FALSE, trim_ws = TRUE)

ZFIN_orthologsDF <- as.data.frame(ZFIN_orthologs)
message(sprintf("Loaded %i records", nrow(ZFIN_orthologsDF)))            # 43 612 records 

HCOP_zfinOrthologs <- read_delim("HCOP_zfinOrthologs.tsv",
                                 delim = "\t", col_names = TRUE, escape_double = FALSE, trim_ws = TRUE)

# not_hcop <- unique(ZFIN_orthologsDF[!ZFIN_orthologsDF$ZFIN_Symbol %in% HCOP_zfinOrthologs$zebrafish_gene_id, ])
not_hcop <- ZFIN_orthologsDF[!ZFIN_orthologsDF$ZFIN_Symbol %in% HCOP_zfinOrthologs$zebrafish_symbol, ] %>%
  distinct(ZFIN_Symbol,  .keep_all = TRUE)

not_hcop <- not_hcop[!not_hcop$ZFIN_ID %in% HCOP_zfinOrthologs$zfin_id, ] %>%
  distinct(ZFIN_ID,  .keep_all = TRUE)


not_zfin <- HCOP_zfinOrthologs[!HCOP_zfinOrthologs$zebrafish_symbol %in% ZFIN_orthologsDF$ZFIN_Symbol, ] %>%
  distinct(zebrafish_symbol,  .keep_all = TRUE)

not_zfin <- not_zfin[!not_zfin$zfin_id %in% ZFIN_orthologsDF$ZFIN_ID, ] %>%
  distinct(zfin_id,  .keep_all = TRUE)


# VERIFY
# Remove duplicates
ZFIN_orthologsDF_unique <- ZFIN_orthologsDF %>% distinct(ZFIN_Symbol, .keep_all = TRUE)
HCOP_zfinOrthologs_unique <- HCOP_zfinOrthologs %>% distinct(zebrafish_symbol, .keep_all = TRUE)

# Verify the counts of unique symbols
n_distinct(ZFIN_orthologsDF_unique$ZFIN_Symbol)  # This should match the first count of 17029
n_distinct(HCOP_zfinOrthologs_unique$zebrafish_symbol)  # This should match the second count of 16977

# Subset and count the unique symbols again
verifyNot_hcop <- ZFIN_orthologsDF_unique[!ZFIN_orthologsDF_unique$ZFIN_Symbol %in% HCOP_zfinOrthologs_unique$zebrafish_symbol, ]

n_distinct(verifyNot_hcop$ZFIN_Symbol)  # This should give us the expected 52

which(verifyNot_hcop$ZFIN_Symbol %in% HCOP_zfinOrthologs$zebrafish_symbol)

which(HCOP_zfinOrthologs$zebrafish_symbol %in% verifyNot_hcop$ZFIN_Symbol)
