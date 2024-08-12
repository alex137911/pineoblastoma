# Human-Zebrafish Orthologs listed in HGNC Comparison of Orthology Predictions (HCOP)
# Note that HCOP (i.e., HGNC Comparison of Orthology Predictions) uses 
# Ensembl v109 as of May 8, 2024
# https://genenames.org/help/hcop/

# Remove objects in workspace
rm(list = ls())

# Required Packages -----------------------------------------------------------
suppressMessages(library(rtracklayer))
suppressMessages(library(readr))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(scales))
suppressMessages(library(tidyr))
suppressMessages(library(stringr))
suppressMessages(library(readxl))
suppressMessages(library(tidyverse))
suppressMessages(library(UpSetR))
suppressMessages(library(data.table))
suppressMessages(library(patchwork))
suppressMessages(library(grDevices))
suppressMessages(library(ggdist))
suppressMessages(library(igraph))
suppressMessages(library(e1071))
suppressMessages(library(openxlsx))
suppressMessages(library(devEMF))
#remotes::install_github("davidsjoberg/ggsankey")
suppressMessages(library(ggsankey))

# Load Data -------------------------------------------------------------------
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

# Farrell Data (https://doi.org/10.1016/j.devcel.2023.11.001)
# Categorized genes based on spatiotemporal variation (Table S3)
multipleSheets <- function(fname) { 
  # Getting info about all sheets 
  sheets <- readxl::excel_sheets(fname) 
  
  # Deal with type mismatches among Excel sheets
  convertData <- function(sheet, fname) {
    df <- readxl::read_excel(fname, sheet = sheet) %>% as.data.frame()
    df$CV.celltype.per.stage.mean.log <- as.numeric(df$CV.celltype.per.stage.mean.log)
    return(df)
  }
  
  # Read the first sheet
  combined_df <- convertData(sheets[1], fname)
  
  # Append subsequent sheets to the data frame
  for (i in 2:length(sheets)) {
    sheet_df <- convertData(sheets[i], fname)
    combined_df <- bind_rows(combined_df, sheet_df)
  }
  
  # Return the combined dataframe
  return(combined_df)
} 

# Path to Excel file
path <- "farrell_categorized-genes-spatiotemporal-variation.xlsx"

# Calling function and storing result
farrell_categorizedGenes <- multipleSheets(path)

farrell_categorizedGenes <- read_excel("farrell_categorized-genes-spatiotemporal-variation.xlsx")

farrell_categorizedGenesDF <- as.data.frame(farrell_categorizedGenes)
message(sprintf("Loaded %i records", nrow(farrell_categorizedGenesDF)))  # 863 records

# Retinoblastoma gene in cancer (WP2446) pathway
# Downloaded May 10, 2024 from https://www.wikipathways.org/pathways/WP2446.html
# Pathway upregulated in in PB tumours (RB subgroup) - from Oliva
WP2446_pathwayGenes <- read_delim("WP2446_pathwayGenes.tsv", delim = "\t", escape_double = FALSE,
                                  trim_ws = TRUE)

# Drop tags
WP2446_pathwayGenes <- WP2446_pathwayGenes %>% mutate(Ensembl = sub("ensembl:", "", Ensembl))
WP2446_pathwayGenes <- WP2446_pathwayGenes %>% mutate(HGNC = sub("hgnc.symbol:", "", HGNC))

# 89 genes
n_distinct(WP2446_pathwayGenes$HGNC)

# Ensembl v109 BioMart Export (Zebrafish Gene IDs)
# True ZF Ensembl IDs from Ensembl
zebrafish_genesMart <- read_delim("zebrafish-genes-mart-export.txt", delim = "\t", escape_double = FALSE, 
                                  trim_ws = TRUE)

zebrafish_genesMartDF <- as.data.frame(zebrafish_genesMart)


human_genesMart <- read_delim("human_ensembl-entrez-mart_export.txt", delim = "\t", escape_double = FALSE,
                              trim_ws = TRUE)

human_genesMartDF <- as.data.frame(human_genesMart)

# ZFIN-Ensembl 1:1 mapping
zfin_ensemblGene <- read_delim("zfin-ensembl_1_to_1_2024.07.09.txt", 
                               delim = "\t", escape_double = FALSE, trim_ws = TRUE)

zfin_ensemblGeneDF <- as.data.frame(zfin_ensemblGene)
zfin_ensemblGeneDF$...5 <- NULL

# Out Dir
outpath <- "C:/Users/acale/OneDrive/Documents/Waterloo BME/Co-op/4C/pineoblastoma/Data/Output/"
outDir  <- sprintf("%s/Output", dirname(outpath))
setwd(outDir)

# Read in Farrell Ortholog list
farrellHuman_zebrafish <- read_delim("farrellHuman_zebrafishOrthologs.tsv",
                                     delim = "\t", escape_double = FALSE, trim_ws = TRUE)

# Quality Check ---------------------------------------------------------------

# HGNC null entries 
no_entries <- sum(hgnc_orthologsDF$zebrafish_entrez_gene == "-")  # 11 645 blank
no_entries <- sum(hgnc_orthologsDF$zebrafish_ensembl_gene == "-") # 2 081  blank
no_entries <- sum(hgnc_orthologsDF$zebrafish_symbol == "-")       # 5 295  blank
no_entries <- sum(hgnc_orthologsDF$zfin_id == "-")                # 16 005 blank

no_entries <- sum(hgnc_orthologsDF$human_symbol == "-")           # 377    blank
no_entries <- sum(hgnc_orthologsDF$human_entrez_gene == "-")      # 416    blank
no_entries <- sum(hgnc_orthologsDF$human_ensembl_gene == "-")     # 0      blank

# Multiple null entries
total_null <- hgnc_orthologs %>%
  filter(grepl("-", zebrafish_symbol) & grepl("-", human_symbol)) # 2 040  blank (both species)

# total_null %>%
#   filter(grepl("-", zebrafish_symbol) & grepl("-", human_symbol)) %>%
#   summarise(total = n())

no_entries <- sum(total_null$human_entrez_gene == "-")            # 86     blank     
no_entries <- sum(total_null$human_ensembl_gene == "-")           # 0      blank

no_entries <- sum(total_null$zebrafish_entrez_gene == "-")        # 1 107  blank
no_entries <- sum(total_null$zfin_id == "-")                      # 507    blank
no_entries <- sum(total_null$zebrafish_ensembl_gene == "-")       # 15     blank

# Define the categories with prefixes
# cDNA clones (si:ch211-, si:dkey-, zgc:, etc.)
cDNA_prefixes <- c("^si:", "^sb:", "^sc:", "^wu:", "^gb:", "zgc:")

# Computationally predicted genes (LOC-, CABZ-, BX-, AL-, etc.)
# zmpste24 is a zebrafish gene (need "zmp:" not "zmp")
predicted_prefixes <- c("^LOC", "^im:", "^CABZ", "^BX", "^AL", "^Gm", "^OTTDARG", 
                        "^XP", "^ENS", "^zmp:")

# Create patterns from prefixes
cDNA_pattern <- paste(cDNA_prefixes, collapse = "|")
predicted_pattern <- paste(predicted_prefixes, collapse = "|")

cDNA <- hgnc_orthologs %>%                                       # 48 335 entries (NOT genes)
  filter(grepl(cDNA_pattern, zebrafish_symbol)) %>%
  summarise(total = n())

predicted <- hgnc_orthologs %>%                                  # 10 290 entries
  filter(grepl(predicted_pattern, zebrafish_symbol)) %>%
  summarise(total = n())

# Putative genes identified through de novo transcriptome assembly by the Lawson lab (XLOC-)
putative <- hgnc_orthologs %>%
  filter(grepl("^XLOC|^TCONS", zebrafish_symbol)) %>%
  summarise(total = n())

# Records with named genes = 132 678 - 48 727 - 10 290 - 5 295 = 68 366

# test <- hgnc_orthologs %>% filter(grepl("^zmp", zebrafish_symbol))
# test <- hgnc_orthologs %>% filter(grepl(":", zebrafish_symbol))
# test <- test %>% filter(!grepl("^si:ch211|^si:ch73|^si:busm1|^si:dkey|^si:hs|^zgc:", zebrafish_symbol))
# test <- test %>% filter(!grepl("^zmp:", zebrafish_symbol))
# 
# unique_tags <- unique(test$zebrafish_symbol)

# UpSet Plot ------------------------------------------------------------------
# Function to process the support column
# Remove instances of multiple assertions made for the same gene pair from 
# the same source (want to consider each unique source only once per ortholog pair)
process_support <- function(support) {
  unique_sources <- unique(trimws(strsplit(support, ",")[[1]]))
  paste(sort(unique_sources), collapse = "&")
}

# Apply the process_support function to each row
hgnc_orthologsDF <- hgnc_orthologsDF %>%
  mutate(support_combination = sapply(support, process_support))

# Count number of intersections per combination
support_counts <- hgnc_orthologsDF %>% count(support_combination)

# Create a named vector for the UpSet plot input
upset_input <- setNames(support_counts$n, support_counts$support_combination)

# Generate the UpSet plot
upset(fromExpression(upset_input),
      nintersects = 15,
      sets = c("EggNOG", "Ensembl", "HomoloGene", "Inparanoid",
               "OMA", "OrthoDB", "OrthoMCL", "NCBI", "Panther",
               "PhylomeDB", "Treefam", "ZFIN"),
      mb.ratio = c(0.65, 0.35),
      order.by = "freq",
      decreasing = TRUE,
      text.scale = 1.1,
      point.size = 2.8,
      set_size.show = FALSE,
      line.size = 1)

# Define a function to generate and save the UpSet plot
save_upsetPlot <- function(input, file, width, height, res) {
  plot <- upset(fromExpression(input),
                nintersects = 15,
                sets = c("EggNOG", "Ensembl", "HomoloGene", "Inparanoid",
                         "OMA", "OrthoDB", "OrthoMCL", "NCBI", "Panther",
                         "PhylomeDB", "Treefam", "ZFIN"),
                mb.ratio = c(0.65, 0.35),
                order.by = "freq",
                decreasing = TRUE,
                number.angles = 0,
                text.scale = c(2, 2, 2, 2, 2, 2),  # size title, size numbers, set size title,
                point.size = 2.8,                  # set size numbers, set names, main bar plot title              
                line.size = 1)
  png(file, width = width, height = height, units = "in", res = res,
      bg = "transparent")
  print(plot)
  dev.off()
}

# Call the function to save the UpSet plot
save_upsetPlot(upset_input, "HCOP Ortholog Source.png", 
               width = 17, height = 9, res = 600)

# Sum the counts for rows containing "ZFIN" or "NCBI" (manually curated)
manual_count <- support_counts %>%
  filter(grepl("ZFIN", support_combination) | grepl("NCBI", support_combination))

# UpSet plot for manually curated/suported entries
manualUpset_input <- setNames(manual_count$n, manual_count$support_combination)

upset(fromExpression(manualUpset_input),
      nintersects = 15,
      sets = c("EggNOG", "Ensembl", "HomoloGene", "Inparanoid",
               "OMA", "OrthoDB", "OrthoMCL", "NCBI", "Panther",
               "PhylomeDB", "Treefam", "ZFIN"),
      mb.ratio = c(0.65, 0.35),
      order.by = "freq",
      decreasing = TRUE,
      text.scale = 1.1,
      point.size = 2.8,
      set_size.show = FALSE,
      line.size = 1)

save_upsetPlot(manualUpset_input, "HCOP Manual Orthologs Source.png", 
               width = 17, height = 9, res = 600)

# Remove Duplicates ------------------------------------------------------------
# Create a new column with the number of databases supporting a particular assertion
hgnc_orthologsDF$support_count <- sapply(strsplit(hgnc_orthologsDF$support_combination, "&"), length)

# Set to data table for faster operations
setDT(hgnc_orthologsDF)

# Subset for testing
# hgnc_orthologsDF[, cluster := NULL]
# subset_hgnc_orthologsDF <- hgnc_orthologsDF[1:50]
# 
# subset_hgnc_orthologsDF <- hgnc_clusters_original[, c("hgnc_id", "human_name", "human_chr", "human_assert_ids",
#                                                        "zebrafish_name", "zebrafish_chr", "zebrafish_assert_ids",
#                                                        "support") := NULL]
# 
# subset_hgnc_orthologsDF <- subset_hgnc_orthologsDF[cluster == "Cluster 509", ]
# 
# subset_hgnc_orthologsDF <- subset_hgnc_orthologsDF %>%
#   group_by(cluster) %>%
#   mutate(cluster_count = n()) %>%
#   ungroup() %>%
#   filter(cluster_count > 1 & cluster_count <= 5) %>%
#   select(-cluster_count)  # Optionally remove the temporary count column
# 
# # Unique IDs to check
# subset_distinct <- subset_hgnc_orthologsDF %>%
#   dplyr::distinct(cluster, .keep_all = TRUE)
# 
# store_cluster <- subset_hgnc_orthologsDF[subset_hgnc_orthologsDF$cluster == "Cluster 3988", ]
# 
# write.xlsx(store_cluster, 
#            file = "cluster_3988.xlsx", colNames = TRUE, rowNames = FALSE)

# Add a column for the original row indices
hgnc_orthologsDF[, row_idx := .I]

hgnc_orthologsDF <- hgnc_orthologsDF %>%
  relocate(row_idx, .before = human_entrez_gene)

# Create copy of DF
hgnc_orthologsDF_backup <- copy(hgnc_orthologsDF)

# Create temporary columns with the transformed case and check for duplicates
# Different databases follow different nomenclature (e.g., ABR vs abr, MCAT vs mcat)
# Need to convert to a standard nomenclature so that ortholog counts are not falsely inflated
hgnc_orthologsDF[, zebrafish_symbol_lower := tolower(zebrafish_symbol)]
hgnc_orthologsDF[, human_symbol_upper := toupper(human_symbol)]

# Identify duplicates
dupes_zfSymbol_lower <- hgnc_orthologsDF[, .N, by = zebrafish_symbol_lower][N > 1]$zebrafish_symbol_lower
dupes_humanSymbol_upper <- hgnc_orthologsDF[, .N, by = human_symbol_upper][N > 1]$human_symbol_upper

# Combine the vectors, remove duplicates from the combined list
zebrafish_dupes <- unique(c(dupes_zfSymbol_lower))
human_dupes <- unique(c(dupes_humanSymbol_upper))

# Update the original columns based on the duplicates found
hgnc_orthologsDF[, zebrafish_symbol := ifelse(tolower(zebrafish_symbol) %in% zebrafish_dupes, tolower(zebrafish_symbol), zebrafish_symbol)]
hgnc_orthologsDF[, human_symbol := ifelse(toupper(human_symbol) %in% human_dupes, toupper(human_symbol), human_symbol)]

# Remove the temporary columns
hgnc_orthologsDF[, c("zebrafish_symbol_lower", "human_symbol_upper") := NULL]

# DOUBLE CHECK FOR ANO5 AND PINX1

# 5 Entrez IDs with conflicting human IDs (but these genes occur in different locations)
# Additional 15 Human Symbols with no Entrez ID mapping
# 21 Human Symbols with conflicting Entrez IDs
# 0 Entrez IDs with no Human Symbol mapping
hcopGene_conflicts <- hgnc_orthologsDF %>%
  group_by(human_symbol) %>%
  summarize(
    unique_ensembl_count = n_distinct(str_trim(human_ensembl_gene[human_ensembl_gene != "-" & str_trim(human_ensembl_gene) != ""])),
    unique_ensembl = paste(unique(str_trim(human_ensembl_gene[human_ensembl_gene != "-" & str_trim(human_ensembl_gene) != ""])), collapse = ", "),
    unique_symbol_count = n_distinct(str_trim(as.character(human_entrez_gene[human_entrez_gene != "-" & str_trim(as.character(human_entrez_gene)) != ""]))),
    unique_symbol = paste(unique(str_trim(as.character(human_entrez_gene[human_entrez_gene != "-" & str_trim(as.character(human_entrez_gene)) != ""]))), collapse = ", ")
  )


sum(hcopGene_conflicts$unique_symbol_count>1)

# 99 Entrez IDs with conflicting zebrafish IDs
# Additional 1017 Zebrafish Symbols with no Entrez ID mapping
# 88 zebrafish symbols with conflicting Entrez IDs
# 750 Entrez IDs with no Zebrafish Symbol mapping
hcopGene_conflicts <- hgnc_orthologsDF %>%
  group_by(zebrafish_symbol) %>%
  summarize(
    unique_ensembl_count = n_distinct(zebrafish_ensembl_gene[zebrafish_ensembl_gene != "-"]),
    unique_ensembl = paste(unique(zebrafish_ensembl_gene[zebrafish_ensembl_gene != "-"]), collapse = ", "),
    unique_symbol_count = n_distinct(zebrafish_entrez_gene[zebrafish_entrez_gene != "-"]),
    unique_symbol = paste(unique(zebrafish_entrez_gene[zebrafish_entrez_gene != "-"]), collapse = ", ")
  )

sum(hcopGene_conflicts$unique_symbol_count>1)

# Missing gene IDs
# 377 genes without Human Symbol
sum(hgnc_orthologsDF$human_symbol == "-")
# 5295 genes without Zebrafish Symbol (11645 without Entrez ID)
sum(hgnc_orthologsDF$zebrafish_symbol == "-")

# All 377 genes also don't have Entrez ID
missing_gene <- subset(hgnc_orthologsDF, human_symbol == "-")
sum(missing_gene$human_entrez_gene == "-")

# 4544/5295 genes also don't have Entrez ID (751 do)
missing_gene <- subset(hgnc_orthologsDF, zebrafish_symbol == "-")
sum(missing_gene$zebrafish_entrez_gene == "-")

# Create generic zebrafish_gene_id using zebrafish symbol + Entrez ID (for 751 genes)
hgnc_orthologsDF[, zebrafish_gene_id := fifelse(zebrafish_symbol != "-",
                                                zebrafish_symbol,
                                                zebrafish_entrez_gene)]

sum(hgnc_orthologsDF$zebrafish_gene_id == "-") # 4544 missing IDs

# Merge non-unique ortholog predictions occurring from different databases
# e.g., appl2-APPL2 (prediction made 3 separate times)
# Group identical (non-unique) ortholog predictions based on gene ID and chromosome of gene
# After grouping identical predictions, remove duplicate assertions made from the same database
# Recount true number of supporting databases
colnames <- colnames(hgnc_orthologsDF)

hgnc_orthologsDF <- hgnc_orthologsDF[, .(
  row_idx = paste(unique(row_idx), collapse = ","),
  human_entrez_gene = paste(setdiff(unique(human_entrez_gene), "-"), collapse = ", "),
  human_ensembl_gene = paste(setdiff(unique(human_ensembl_gene), "-"), collapse = ", "),
  hgnc_id = paste(setdiff(unique(hgnc_id), "-"), collapse = ", "),
  human_name = paste(setdiff(unique(human_name), "-"), collapse = ", "),
  human_assert_ids = paste(setdiff(unique(human_assert_ids), "-"), collapse = ", "),
  zebrafish_entrez_gene = paste(setdiff(unique(zebrafish_entrez_gene), "-"), collapse = ", "),
  zebrafish_ensembl_gene = paste(setdiff(unique(zebrafish_ensembl_gene), "-"), collapse = ", "),
  zebrafish_symbol = paste(setdiff(unique(zebrafish_symbol), "-"), collapse = ", "),
  zfin_id = paste(setdiff(unique(zfin_id), "-"), collapse = ", "),
  zebrafish_name = paste(setdiff(unique(zebrafish_name), "-"), collapse = ", "),
  zebrafish_assert_ids = paste(setdiff(unique(zebrafish_assert_ids), "-"), collapse = ", "),
  support = paste(unique(support), collapse = ","),
  support_combination = paste(unique(unlist(strsplit(support_combination, "&"))), collapse = "&"),
  support_count = length(unique(unlist(strsplit(support_combination, "&"))))
), by = .(human_symbol, zebrafish_gene_id, zebrafish_chr, human_chr)]

# Reorder the columns to match the original order
setcolorder(hgnc_orthologsDF, colnames)

# Identify rows with more than one ID (i.e., collapsed rows)
# 954 collapsed rows (non-unique predictions)
# 132,678 rows --> 129,593 rows (-3,085 rows)
collapsed_rows <- hgnc_orthologsDF[sapply(strsplit(row_idx, ","), length) > 1]

# 18,396 human genes
n_distinct(hgnc_orthologsDF$human_symbol)
# 23,900 zebrafish genes
n_distinct(hgnc_orthologsDF$zebrafish_gene_id)

# Copy
hgnc_orthologsDF_backup <- hgnc_orthologsDF

# Drop rows with no gene ID
hgnc_orthologsDF <- hgnc_orthologsDF[human_symbol != "-"]        # 129,240 rows (-353 rows)
n_distinct(hgnc_orthologsDF$human_symbol)                        # 18,396 human genes
hgnc_orthologsDF <- hgnc_orthologsDF[zebrafish_gene_id != "-"]   # 127,268 rows (-1972 rows)
n_distinct(hgnc_orthologsDF$zebrafish_gene_id)                   # 23, 898 (-2 genes)

# si:ch211-171h4.7, si:dkey-97l20.6
unique(hgnc_orthologsDF_backup$zebrafish_gene_id[!(hgnc_orthologsDF_backup$zebrafish_gene_id %in% hgnc_orthologsDF$zebrafish_gene_id)])
unique(hgnc_orthologsDF_backup$human_symbol[!(hgnc_orthologsDF_backup$human_symbol %in% hgnc_orthologsDF$human_symbol)])

# Consensus Ensembl gene ID ----------------------------------------------------

# Check Ensembl ID against Entrez ID (zebrafish)
entrez_ensemblConflict <- zebrafish_genesMartDF %>%
  group_by(`NCBI gene (formerly Entrezgene) ID`) %>%
  summarise(across(everything(), ~ paste(unique(na.omit(.)), collapse = ", ")),
            .groups = "drop")

# Rename columns for merge
entrez_ensemblConflict <- entrez_ensemblConflict %>%
  dplyr::rename(zebrafish_entrez_gene = `NCBI gene (formerly Entrezgene) ID`, 
                entrez_ensembl_gene = `Gene stable ID`)

# Merge
entrez_ensemblConflict <- entrez_ensemblConflict %>%
  mutate(zebrafish_entrez_gene = as.character(zebrafish_entrez_gene))

hgnc_orthologsDF <- hgnc_orthologsDF %>%
  left_join(entrez_ensemblConflict %>% select(zebrafish_entrez_gene, entrez_ensembl_gene), 
            by = "zebrafish_entrez_gene")

hgnc_orthologsDF <- hgnc_orthologsDF %>%
  relocate(entrez_ensembl_gene, .after = zebrafish_ensembl_gene)

# DELETE?
# zfin_entrezConflict <- zfin_entrezGeneDF %>%
#   group_by(`ZFIN ID`) %>%
#   summarise(across(everything(), ~ paste(unique(na.omit(.)), collapse = ", ")),
#             .groups = "drop")

# Check Ensembl ID against ZFIN
zfin_ensemblConflict <- zfin_ensemblGeneDF %>%
  group_by(`ZFIN ID`) %>%
  summarise(across(everything(), ~ paste(unique(.), collapse = ", ")),
            .groups = "drop")

# Rename columns for merge
zfin_ensemblConflict <- zfin_ensemblConflict %>%
  dplyr::rename(zfin_id = `ZFIN ID`, zfin_ensembl_gene = `Ensembl ID`)

# Merge
hgnc_orthologsDF <- hgnc_orthologsDF %>%
  left_join(zfin_ensemblConflict %>% select(zfin_id, zfin_ensembl_gene), by = "zfin_id")

hgnc_orthologsDF <- hgnc_orthologsDF %>%
  relocate(zfin_ensembl_gene, .after = zebrafish_ensembl_gene)

# Find the most frequent gene IDs across three zebrafish columns
# When the maximum frequency is n = 1, consider keeping only zebrafish_ensembl_gene
find_consensusEnsembl <- function(row) {
  genes <- c(str_split(row[["zebrafish_ensembl_gene"]], ",")[[1]], 
             str_split(row[["zfin_ensembl_gene"]], ",")[[1]], 
             str_split(row[["entrez_ensembl_gene"]], ",")[[1]])
  
  genes <- genes[!is.na(genes) & genes != ""]
  
  if (length(genes) == 0) {
    return(NA)
  }
  
  gene_table <- sort(table(genes), decreasing = TRUE)
  most_frequent <- names(gene_table)[gene_table == max(gene_table)]
  paste(most_frequent, collapse = ", ")
}

# Apply function
hgnc_orthologsDF[, consensus_ensembl_gene := apply(.SD, 1, find_consensusEnsembl), 
                 .SDcols = c("zebrafish_ensembl_gene", "zfin_ensembl_gene", "entrez_ensembl_gene")]

hgnc_orthologsDF <- hgnc_orthologsDF %>%
  relocate(consensus_ensembl_gene, .before = zebrafish_ensembl_gene)

# Replace "-" (representing NULL entries) with NA
hgnc_orthologsDF[hgnc_orthologsDF == "-"] <- NA

# Clustering -------------------------------------------------------------------
# Columns to check for overlaps
columns <- c("human_entrez_gene", "human_ensembl_gene", "human_symbol",
             "zebrafish_entrez_gene", "consensus_ensembl_gene", "zfin_id", "zebrafish_symbol")

# Melt the data.table to long format
long_hgnc <- melt(hgnc_orthologsDF, id.vars = "row_idx", measure.vars = columns, 
                  na.rm = TRUE, variable.name = "gene_type", value.name = "gene_id")

# Remove rows with NA gene_id (if any)
long_hgnc <- long_hgnc[!is.na(gene_id)]

# Create edges for the graph by connecting rows that share the same gene_id
# Row indices are nodes, edges are created between nodes that share a gene ID in any of the columns
edges <- long_hgnc[, if (.N > 1) .(from = combn(row_idx, 2)[1,], 
                                   to = combn(row_idx, 2)[2,]), by = gene_id]

# Ensure all rows are included in the graph (even those without edges)
all_nodes  <- unique(long_hgnc$row_idx)
edge_nodes <- unique(c(edges$from, edges$to))
isolated_nodes <- setdiff(all_nodes, edge_nodes)

# Using the edges, create graph where each unique gene id is connected to its corresponding row indices
cluster_graph <- graph_from_data_frame(edges[, .(from, to)], directed = FALSE, vertices = all_nodes)

# Find the clusters
clusters <- components(cluster_graph)

# Create a mapping from row_idx to cluster membership
cluster_membership <- data.table(row_idx = all_nodes, cluster = clusters$membership)

# Map the cluster membership back to the original data.table
hgnc_orthologsDF <- merge(hgnc_orthologsDF, cluster_membership, by = "row_idx", all.x = TRUE)
hgnc_orthologsDF[, cluster := paste0("Cluster ", cluster)]
hgnc_orthologsDF[, cluster := ifelse(is.na(cluster), "Cluster NA", cluster)]

# Create a summary table of unique clusters and their counts
cluster_summary <- hgnc_orthologsDF[, .N, by = cluster][order(cluster)]

# Sort by numerical value
cluster_summary[, cluster_num := as.integer(sub("Cluster ", "", cluster))]
cluster_summary <- cluster_summary[order(cluster_num)]
cluster_summary[, cluster_num := NULL]
setnames(cluster_summary, "N", "count")

# Function to verify a specific cluster is correct
examine_cluster <- function(cluster_num, df) {
  cluster_label <- paste0("Cluster ", cluster_num)
  rows_in_cluster <- hgnc_orthologsDF[cluster == cluster_label]
  
  result <- data.table(cluster = character(), column = character(), shared_gene_ids = character())
  
  for (col in columns) {
    gene_ids <- rows_in_cluster[[col]]
    gene_ids <- gene_ids[!is.na(gene_ids)]  # Exclude NA values
    shared_gene_ids <- gene_ids[gene_ids %in% gene_ids[duplicated(gene_ids)]]
    if (length(shared_gene_ids) > 0) {
      result <- rbind(result, data.table(cluster = cluster_label, column = col, 
                                         shared_gene_ids = toString(unique(shared_gene_ids))))
    }
  }
  
  return(result)
}

# Check
result <- examine_cluster(509)

# Calculate kurtosis for the count column
kurtosis_value <- kurtosis(cluster_summaryDF$count)

# Density Plot
cluster_summaryDF <- as.data.frame(cluster_summary)

plot<-cluster_summaryDF %>%
  ggplot(aes(x=count)) +
  geom_density(fill = "#F6A9CA", color = "#EE5496", alpha = 0.8) +
  labs(x = "Node Count per Ortholog Cluster", y = "Density") +
  # coord_cartesian(xlim = c(0, 5000)) +
  xlim(0, 50) +
  theme_minimal() +
  theme(axis.title.x = element_text(size = 18, margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size = 18, margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 12))
  # + annotate("text", x = Inf, y = Inf, label = paste("Kurtosis:", round(kurtosis_value, 2)),
  #          hjust = 1.2, vjust = 5, size = 5, color = "black")

# Histogram
# Create a frequency table for the filtered counts
frequency_table <- cluster_summaryDF %>%
  #filter(count <= 50) %>%
  group_by(count) %>%
  summarize(frequency = n())

plot<-ggplot(cluster_summaryDF, aes(x=count)) +
  geom_histogram(binwidth = 1, fill = "#FFE69F", color = "#FFC625", alpha = 0.8) +
  labs(x = "Node Count per Ortholog Cluster", y = "Frequency") +
  coord_cartesian(xlim = c(0, 50)) +
  theme_minimal() +
  theme(axis.title.x = element_text(size = 18, margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size = 18, margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12))

ggsave("Histogram - Clusters.png", plot = plot, device = "png", 
       width = 9.50, height = 8.19, dpi = 600)

# # Store clustered data (June 26, 2024)
# write.table(hgnc_orthologsDF, file = "clustered_hgncOrthologs.tsv",
#             sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

# write.table(cluster_summary, file = "clusterRanks.tsv",
#             sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

# Rows where consensus Ensembl ID is different from HCOP Ensembl ID
# 516 rows
consensus_conflict <- hgnc_orthologsDF[hgnc_orthologsDF$consensus_ensembl_gene != hgnc_orthologsDF$zebrafish_ensembl_gene, ]

subset_df <- hgnc_orthologsDF[grepl(",", hgnc_orthologsDF$consensus_ensembl_gene), ]

# Filter HGNC Data -------------------------------------------------------------
# Create generic human_gene_id column
# hgnc_orthologsDF[, human_gene_id := fifelse(!is.na(human_symbol),
#                                             human_symbol,
#                                             fifelse(!is.na(human_entrez_gene),
#                                                     human_entrez_gene,
#                                                     human_ensembl_gene))]
# 
# # hgnc_orthologsDF[, human_gene_id := fifelse(human_symbol != "-",
# #                                             human_symbol,
# #                                             fifelse(human_entrez_gene != "-",
# #                                                     human_entrez_gene,
# #                                                     human_ensembl_gene))]
# 
# # Check that only 377 entries have Ensembl ID (number of NULL entries for HGNC)
# hgnc_orthologsDF %>%
#   filter(grepl("ENSG00000", human_gene_id)) %>%
#   summarise(total = n())
# 
# # Check if zebrafish_ensembl_gene is in consensus_ensembl_gene
# check_ensemblConsensus <- function(zebrafish_ensembl, consensus_ensembl) {
#   zebrafish_genes <- unlist(str_split(zebrafish_ensembl, ","))
#   consensus_genes <- unlist(str_split(consensus_ensembl, ","))
#   
#   # Check if any zebrafish_ensembl_gene is not in consensus_ensembl_gene
#   any(!zebrafish_genes %in% consensus_genes)
# }
# 
# # Apply function
# hgnc_orthologsDF[, not_in_consensus := check_ensemblConsensus(zebrafish_ensembl_gene, consensus_ensembl_gene), 
#                  by = 1:nrow(hgnc_orthologsDF)]
# 
# # 5822 rows where HCOP ZF Ensembl Gene is not the consensus Ensembl Gene
# sum(hgnc_orthologsDF$not_in_consensus, na.rm = TRUE)
# 
# hgnc_orthologsDF <- hgnc_orthologsDF %>%
#   relocate(not_in_consensus, .after = zebrafish_ensembl_gene)
# 
# # Create generic zebrafish_gene_id column
# hgnc_orthologsDF[, zebrafish_gene_id := fifelse(!is.na(zebrafish_ensembl_gene),
#                                                 zebrafish_ensembl_gene,
#                                                 fifelse(!is.na(zfin_id),
#                                                         zfin_id,
#                                                         fifelse(!is.na(zebrafish_entrez_gene),
#                                                                 zebrafish_entrez_gene,
#                                                                 zebrafish_symbol)))]
# 
# # hgnc_orthologsDF[, zebrafish_gene_id := fifelse(zebrafish_symbol != "-",
# #                                                 zebrafish_symbol,
# #                                                 fifelse(zebrafish_entrez_gene != "-",
# #                                                         zebrafish_entrez_gene,
# #                                                         fifelse(zfin_id != "-",
# #                                                                 zfin_id,
# #                                                                 zebrafish_ensembl_gene)))]
# 
# # Check how many entries have Ensembl ID (4529)
# hgnc_orthologsDF %>%
#   filter(grepl("ENSDARG000", zebrafish_gene_id)) %>%
#   summarise(total = n())
# 
# # Create a copy of the original data
# hgnc_orthologsDF_original <- copy(hgnc_orthologsDF)
# 
# # Merge non-unique ortholog predictions occurring from different databases
# # e.g., appl2-APPL2 (prediction made 3 separate times)
# # Group identical (non-unique) ortholog predictions based on gene ID and chromosome of gene
# # After grouping identical predictions, remove duplicate assertions made from the same database
# # Recount true number of supporting databases
# hgnc_orthologsDF <- hgnc_orthologsDF[, .(
#   row_id = paste(unique(row_id), collapse = ", "),
#   human_entrez_gene = paste(na.omit(unique(human_entrez_gene)), collapse = ", "),
#   human_ensembl_gene = paste(na.omit(unique(human_ensembl_gene)), collapse = ", "),
#   hgnc_id = paste(na.omit(unique(hgnc_id)), collapse = ", "),
#   human_name = paste(na.omit(unique(human_name)), collapse = ", "),
#   human_symbol = paste(na.omit(unique(human_symbol)), collapse = ", "),
#   human_assert_ids = paste(na.omit(unique(human_assert_ids)), collapse = ", "),
#   zebrafish_entrez_gene = paste(na.omit(unique(zebrafish_entrez_gene)), collapse = ", "),
#   zebrafish_ensembl_gene = paste(na.omit(unique(zebrafish_ensembl_gene)), collapse = ", "),
#   zfin_id = paste(na.omit(unique(zfin_id)), collapse = ", "),
#   zebrafish_name = paste(na.omit(unique(zebrafish_name)), collapse = ", "),
#   zebrafish_symbol = paste(na.omit(unique(zebrafish_symbol)), collapse = ", "),
#   zebrafish_assert_ids = paste(na.omit(unique(zebrafish_assert_ids)), collapse = ", "),
#   support = paste(unique(support), collapse = ", "),
#   support_combination = paste(unique(unlist(strsplit(support_combination, "&"))), collapse = "&"),
#   support_count = length(unique(unlist(strsplit(support_combination, "&")))),
#   cluster = paste(na.omit(unique(cluster)), collapse = ", ")
# ), by = .(human_gene_id, zebrafish_gene_id, zebrafish_chr, human_chr)]
# 
# 
# # # Store as .tsv (July 5, 2024)
# # write.table(hgnc_orthologsDF, file = "hgnc_orthologsFiltered.tsv",
# #             sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

# HGNC Data -------------------------------------------------------------------

# CHECK THAT NA.OMIT IS WORKING FOR BLANK ENTRIES (E.G., ENSDARG00000096473)

# Count occurrences and list all matching zebrafish genes for each human gene
# HGNChuman_zebrafish <- hgnc_orthologsDF %>%
#   group_by(human_symbol) %>%
#   summarise(
#     # Concatenate all unique human symbols associated with each human gene ID, excluding NA
#     human_symbol = paste(na.omit(unique(human_symbol)), collapse = ", "),
#     human_entrez_gene = paste(na.omit(unique(human_entrez_gene)), collapse = ", "),
#     human_ensembl_gene = paste(na.omit(unique(human_ensembl_gene)), collapse = ", "),
#     zebrafish_count = n_distinct(zebrafish_gene_id),
#     # Concatenate all unique zebrafish gene IDs
#     zebrafish_gene_id = paste(na.omit(unique(zebrafish_gene_id)), collapse = ", "),
#     zebrafish_symbol = paste(na.omit(unique(zebrafish_symbol)), collapse = ", "),
#     zebrafish_entrez_gene = paste(na.omit(unique(zebrafish_entrez_gene)), collapse = ", "),
#     zfin_id = paste(na.omit(unique(zfin_id)), collapse = ", "),
#     zebrafish_ensembl_gene = paste(na.omit(unique(zebrafish_ensembl_gene)), collapse = ", "),
#     cluster = paste(na.omit(unique(cluster)), collapse = ", "),
#     average_support = round(mean(support_count, na.rm = TRUE), 1),
#     median_support = round(median(support_count, na.rm = TRUE), 1),
#     # Determine weighted average support for number of zebrafish orthologs per human gene
#     weighted_average = round(sum(support_count * (support_count / 12), na.rm = TRUE) / sum(support_count / 12, na.rm = TRUE), 1),
#     .groups = 'drop'
#   ) %>%
#   mutate(
#     composite_score = (weighted_average + median_support) / 2,
#     relation_type = case_when(
#       zebrafish_count == 1 ~ "One to One",
#       zebrafish_count > 1 ~ "One to Many",
#       TRUE ~ "No Match"
#     )
#   )

HGNChuman_zebrafish <- hgnc_orthologsDF %>%
  group_by(human_symbol) %>%
  # 18,210 humana genes
  summarise(
    # Concatenate all unique human symbols associated with each human gene ID, excluding "-"
    human_symbol = paste(setdiff(unique(human_symbol), "-"), collapse = ", "),
    human_entrez_gene = paste(setdiff(unique(human_entrez_gene), "-"), collapse = ", "),
    human_ensembl_gene = paste(setdiff(unique(human_ensembl_gene), "-"), collapse = ", "),
    zebrafish_count = n_distinct(zebrafish_gene_id),
    # Concatenate all unique zebrafish gene IDs
    zebrafish_gene_id = paste(setdiff(unique(zebrafish_gene_id), "-"), collapse = ", "),
    # zebrafish_symbol = paste(setdiff(unique(zebrafish_symbol), "-"), collapse = ", "),
    zebrafish_entrez_gene = paste(setdiff(unique(zebrafish_entrez_gene), "-"), collapse = ", "),
    zfin_id = paste(setdiff(unique(zfin_id), "-"), collapse = ", "),
    consensus_ensembl_gene = paste(setdiff(unique(consensus_ensembl_gene), "-"), collapse = ", "),
    average_support = round(mean(support_count, na.rm = TRUE), 1),
    median_support = round(median(support_count, na.rm = TRUE), 1),
    # Determine weighted average support for number of zebrafish orthologs per human gene
    weighted_average = round(sum(support_count * (support_count / 12), na.rm = TRUE) / sum(support_count / 12, na.rm = TRUE), 1),
    .groups = 'drop'
  ) %>%
  mutate(
    composite_score = (weighted_average + median_support) / 2,
    relation_type = case_when(
    zebrafish_count == 1 ~ "One to One",
    zebrafish_count > 1 ~ "One to Many",
    TRUE ~ "No Match"
    )
  )

# Identify rows with more than one entry in the cluster column
# multiple_clusters <- HGNChuman_zebrafish %>%
#   filter(str_detect(cluster, ","))

# Count occurrences and list all matching human genes for each zebrafish gene
HGNCzebrafish_human <- hgnc_orthologsDF %>%
  group_by(zebrafish_gene_id) %>%
  summarise(
    zebrafish_entrez_gene = paste(setdiff(unique(zebrafish_entrez_gene), "-"), collapse = ", "),
    zfin_id = paste(setdiff(unique(zfin_id), "-"), collapse = ", "),
    consensus_ensembl_gene = paste(setdiff(unique(consensus_ensembl_gene), "-"), collapse = ", "),
    human_count = n_distinct(human_symbol),
    human_symbol = paste(setdiff(unique(human_symbol), "-"), collapse = ", "),
    human_entrez_gene = paste(setdiff(unique(human_entrez_gene), "-"), collapse = ", "),
    human_ensembl_gene = paste(setdiff(unique(human_ensembl_gene), "-"), collapse = ", "),
    average_support = round(mean(support_count, na.rm = TRUE), 1),
    median_support = round(median(support_count, na.rm = TRUE), 1),
    weighted_average = round(sum(support_count * (support_count / 12), na.rm = TRUE) / sum(support_count / 12, na.rm = TRUE), 1),
    .groups = 'drop'
  ) %>%
  mutate(
    composite_score = (weighted_average + median_support) / 2,
    relation_type = case_when(
      human_count == 1 ~ "One to One",
      human_count > 1 ~ "Many to One",
      TRUE ~ "No Match"
    )
  )

# Single Gene Ortholog analysis
# Filter to keep only unique, one-to-one mappings in both datasets
one_to_one_human <- HGNChuman_zebrafish %>%
  filter(zebrafish_count == 1) %>%
  dplyr::select(human_symbol, zebrafish_gene_id)

one_to_one_zebrafish <- HGNCzebrafish_human %>%
  filter(human_count == 1) %>%
  dplyr::select(zebrafish_gene_id, human_symbol)

# Merge the datasets on the gene IDs, ensuring only those that have a one-to-one relationship in both are kept
SGO_analysis <- inner_join(one_to_one_human, one_to_one_zebrafish, 
                           by = c("human_symbol" = "human_symbol", 
                                  "zebrafish_gene_id" = "zebrafish_gene_id"))

# Mark these as True SGOs
SGO_status <- SGO_analysis %>%
  mutate(SGO_status = "True SGO") %>%
  dplyr::select(human_symbol, zebrafish_gene_id, SGO_status)

# Join SGO status to ortholog data
HGNChuman_zebrafish <- HGNChuman_zebrafish %>% 
  left_join(SGO_status[, c("human_symbol", "SGO_status")], by = "human_symbol") %>%
  mutate(SGO_status = ifelse(is.na(SGO_status), "Not SGO", SGO_status))

HGNCzebrafish_human <- HGNCzebrafish_human %>% 
  left_join(SGO_status[, c("zebrafish_gene_id", "SGO_status")], by = "zebrafish_gene_id") %>%
  mutate(SGO_status = ifelse(is.na(SGO_status), "Not SGO", SGO_status))

# Average orthologs
average_one_to_many <- HGNChuman_zebrafish %>%
  filter(relation_type == "One to Many") %>%        # Filter to keep only 'One to Many' rows
  summarise(average_count = mean(zebrafish_count),
            median_count = median(zebrafish_count))  

average_one_to_many <- HGNCzebrafish_human %>%
  filter(relation_type == "Many to One") %>%        # Filter to keep only 'Many to One' rows
  summarise(average_count = mean(human_count),
            median_count = median(human_count))    

print(average_one_to_many) # 10.9 (3) zebrafish/9.8 (3) humans

HGNCzebrafish_human %>% filter(average_support > 6) %>% summarise(count = n())

# Store as .tsv (August 5, 2024)
# write.table(HGNChuman_zebrafish, file = "HGNChuman_zebrafishOrthologs-final.tsv",
#             sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
# 
# write.table(HGNCzebrafish_human, file = "HGNCzebrafish_humanOrthologs-final.tsv",
#             sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

# HGNC Filter "Many to Many" --------------------------------------------------
# Filter for predictions with manual support
manualHGNC_orthologs <- hgnc_orthologsDF %>% 
  filter(grepl("ZFIN", hgnc_orthologsDF$support) | grepl("NCBI", hgnc_orthologsDF$support))

setDT(manualHGNC_orthologs)
setDT(HGNChuman_zebrafish)
setDT(HGNCzebrafish_human)

# Create a subsets where relation_type is "Many to Many"
subset_HGNChuman_zebrafish <- HGNChuman_zebrafish[relation_type == "One to Many"]
subset_HGNCzebrafish_human <- HGNCzebrafish_human[relation_type == "Many to One"]

# Drop the "SGO_status" column if present
subset_HGNChuman_zebrafish[, SGO_status := NULL]
subset_HGNCzebrafish_human[, SGO_status := NULL]

# Ensure human_gene_id columns are of the same type (character)
manualHGNC_orthologs[, human_symbol := as.character(human_symbol)]
subset_HGNChuman_zebrafish[, human_symbol := as.character(human_symbol)]

# Create a subset of manualHGNC_orthologs where human_gene_id appears in subset
# Where multiple orthology assertions are made, only assertions with manual support are kept
subset_manualHGNC <- manualHGNC_orthologs[human_symbol %in% subset_HGNChuman_zebrafish$human_symbol]
subset_manualHGNC <- subset_manualHGNC[zebrafish_gene_id != "-"]        # 8,129 genes

# Subset of orthologs without manual support
# i.e., assertions from the subset of "Many to Many" but no manual support
subsetNot_manualHGNC <- subset_HGNChuman_zebrafish[!(human_symbol %in% subset_manualHGNC$human_symbol)]
subsetNot_manualHGNC <- hgnc_orthologsDF[human_symbol %in% subsetNot_manualHGNC$human_symbol]
subsetNot_manualHGNC <- subsetNot_manualHGNC[zebrafish_gene_id != "-"]  # 2,887 genes (8,129 + 2,887 = 11,016, correct)

# Append data
manualHGNC_bind <- rbind(subset_manualHGNC, subsetNot_manualHGNC)

# Count occurrences and list all matching zebrafish genes for each human gene
subset_HGNChuman_zebrafish <- manualHGNC_bind %>%
  # 11 016 unique human gene IDs
  group_by(human_symbol) %>%
  summarise(
    # Concatenate all unique human symbols associated with each human gene ID, excluding "-"
    human_entrez_gene = paste(setdiff(unique(human_entrez_gene), "-"), collapse = ", "),
    human_ensembl_gene = paste(setdiff(unique(human_ensembl_gene), "-"), collapse = ", "),
    zebrafish_count = n_distinct(zebrafish_gene_id),
    # Concatenate all unique zebrafish gene IDs
    zebrafish_gene_id = paste(setdiff(unique(zebrafish_gene_id), "-"), collapse = ", "),
    zebrafish_entrez_gene = paste(setdiff(unique(zebrafish_entrez_gene), "-"), collapse = ", "),
    zfin_id = paste(setdiff(unique(zfin_id), "-"), collapse = ", "),
    consensus_ensembl_gene = paste(setdiff(unique(consensus_ensembl_gene), "-"), collapse = ", "),
    average_support = round(mean(support_count, na.rm = TRUE), 1),
    median_support = round(median(support_count, na.rm = TRUE), 1),
    # Determine weighted average support for number of zebrafish orthologs per human gene
    weighted_average = round(sum(support_count * (support_count / 12), na.rm = TRUE) / sum(support_count / 12, na.rm = TRUE), 1),
    .groups = 'drop'
  ) %>%
  mutate(
    composite_score = (weighted_average + median_support) / 2,
    relation_type = case_when(
      zebrafish_count == 1 ~ "One to One",
      zebrafish_count > 1 ~ "One to Many",
      TRUE ~ "No Match"
    )
  )

# Identify dropped zebrafish orthologs
setDT(subset_HGNChuman_zebrafish)

# Merge to compare zebrafish counts
comparison_df <- merge(
  subset_HGNChuman_zebrafish[, .(human_symbol, new_zebrafish_count = zebrafish_count, new_zebrafish_gene_id = zebrafish_gene_id)],
  HGNChuman_zebrafish[, .(human_symbol, original_zebrafish_count = zebrafish_count, original_zebrafish_gene_id = zebrafish_gene_id)],
  by = "human_symbol",
  all = TRUE
)

# Drop rows where filtering was not applied (i.e., one-to-one relationship)
comparison_df <- comparison_df[!is.na(new_zebrafish_count)]

# Identify differences in zebrafish counts
comparison_df[, zebrafish_count_difference := new_zebrafish_count - original_zebrafish_count]

# List dropped zebrafish genes
comparison_df[, dropped_zebrafish_orthologs := mapply(function(new_genes, original_genes) {
  new_genes <- strsplit(new_genes, ", ")[[1]]
  original_genes <- strsplit(original_genes, ", ")[[1]]
  paste(setdiff(original_genes, new_genes), collapse = ", ")
}, new_zebrafish_gene_id, original_zebrafish_gene_id, SIMPLIFY = TRUE)]

# Create a subset where there has been a change in zebrafish count
filtered_changesHZ <- comparison_df[
  zebrafish_count_difference != 0 | dropped_zebrafish_orthologs != "",
  .(human_symbol, zebrafish_count_difference, dropped_zebrafish_orthologs)
]

# Store as .tsv (August 5, 2024)
# write.table(filtered_changesHZ, file = "droppedHuman_zebrafishOrthologs.tsv",
#             sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

# Ensure zebrafish_gene_id columns are of the same type (character)
manualHGNC_orthologs[, zebrafish_gene_id := as.character(zebrafish_gene_id)]
subset_HGNCzebrafish_human[, zebrafish_gene_id := as.character(zebrafish_gene_id)]

# Create a subset of manualHGNC_orthologs where zebrafish_gene_id appears in subset
subset_manualHGNC <- manualHGNC_orthologs[zebrafish_gene_id %in% subset_HGNCzebrafish_human$zebrafish_gene_id]
subset_manualHGNC <- subset_manualHGNC[zebrafish_gene_id != "-"]        # 8,234 genes

# Subset of orthologs without manual support
# i.e., assertions from the subset of "Many to Many" but no manual support
subsetNot_manualHGNC <- subset_HGNCzebrafish_human[!(zebrafish_gene_id %in% subset_manualHGNC$zebrafish_gene_id)]
subsetNot_manualHGNC <- hgnc_orthologsDF[zebrafish_gene_id %in% subsetNot_manualHGNC$zebrafish_gene_id]
subsetNot_manualHGNC <- subsetNot_manualHGNC[zebrafish_gene_id != "-"]  # 3,509 genes (8,234 + 3,509 = 11,743, correct)

# Append data
manualHGNC_bind <- rbind(subset_manualHGNC, subsetNot_manualHGNC)

# Count occurrences and list all matching human genes for each zebrafish gene
subset_HGNCzebrafish_human <- manualHGNC_bind %>%
  # 11 743 unique zebrafish gene IDs
  group_by(zebrafish_gene_id) %>%
  summarise(
    # Concatenate all unique zebrafish symbols associated with each zebrafish gene ID, excluding "-"
    zebrafish_entrez_gene = paste(setdiff(unique(zebrafish_entrez_gene), "-"), collapse = ", "),
    zfin_id = paste(setdiff(unique(zfin_id), "-"), collapse = ", "),
    consensus_ensembl_gene = paste(setdiff(unique(consensus_ensembl_gene), "-"), collapse = ", "),
    human_count = n_distinct(human_symbol),
    # Concatenate all unique human IDs
    human_symbol = paste(setdiff(unique(human_symbol), "-"), collapse = ", "),
    human_entrez_gene = paste(setdiff(unique(human_entrez_gene), "-"), collapse = ", "),
    human_ensembl_gene = paste(setdiff(unique(human_ensembl_gene), "-"), collapse = ", "),
    average_support = round(mean(support_count, na.rm = TRUE), 1),
    median_support = round(median(support_count, na.rm = TRUE), 1),
    # Determine weighted average support for number of human orthologs per zebrafish gene
    weighted_average = round(sum(support_count * (support_count / 12), na.rm = TRUE) / sum(support_count / 12, na.rm = TRUE), 1),
    .groups = 'drop'
  ) %>%
  mutate(
    composite_score = (weighted_average + median_support) / 2,
    relation_type = case_when(
      human_count == 1 ~ "One to One",
      human_count > 1 ~ "Many to One",
      TRUE ~ "No Match"
    )
  )

# Identify dropped human orthologs
setDT(subset_HGNCzebrafish_human)

# Merge to compare human counts
comparison_df <- merge(
  subset_HGNCzebrafish_human[, .(zebrafish_gene_id, new_human_count = human_count, new_human_gene_id = human_symbol)],
  HGNCzebrafish_human[, .(zebrafish_gene_id, original_human_count = human_count, original_human_gene_id = human_symbol)],
  by = "zebrafish_gene_id",
  all = TRUE
)

# Drop rows where filtering was not applied (i.e., one-to-one relationship)
comparison_df <- comparison_df[!is.na(new_human_count)]

# Identify differences in human counts
comparison_df[, human_count_difference := new_human_count - original_human_count]

# List dropped human genes
comparison_df[, dropped_human_orthologs := mapply(function(new_genes, original_genes) {
  new_genes <- strsplit(new_genes, ", ")[[1]]
  original_genes <- strsplit(original_genes, ", ")[[1]]
  paste(setdiff(original_genes, new_genes), collapse = ", ")
}, new_human_gene_id, original_human_gene_id, SIMPLIFY = TRUE)]

# Create a subset where there has been a change in human count
filtered_changesZH <- comparison_df[
  human_count_difference != 0 | dropped_human_orthologs != "",
  .(zebrafish_gene_id, human_count_difference, dropped_human_orthologs)
]

# Store as .tsv (August 5, 2024)
# write.table(filtered_changesZH, file = "droppedZebrafish_humanOrthologs.tsv",
#             sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

# Check if there are any new SGOs
one_to_one_human <- subset_HGNChuman_zebrafish %>%
  filter(zebrafish_count == 1) %>%
  dplyr::select(human_symbol, zebrafish_gene_id)

one_to_one_zebrafish <- subset_HGNCzebrafish_human %>%
  filter(human_count == 1) %>%
  dplyr::select(zebrafish_gene_id, human_symbol)

# Merge the datasets on the gene IDs, ensuring only those that have a one-to-one relationship in both are kept
SGO_analysis <- inner_join(one_to_one_human, one_to_one_zebrafish, 
                           by = c("human_symbol" = "human_symbol", 
                                  "zebrafish_gene_id" = "zebrafish_gene_id"))

# Mark these as True SGOs
SGO_status <- SGO_analysis %>%
  mutate(SGO_status = "True SGO") %>%
  dplyr::select(human_symbol, zebrafish_gene_id, SGO_status)

# Join SGO status to ortholog data
subset_HGNChuman_zebrafish <- subset_HGNChuman_zebrafish %>% 
  left_join(SGO_status[, c("human_symbol", "SGO_status")], by = "human_symbol") %>%
  mutate(SGO_status = ifelse(is.na(SGO_status), "Not SGO", SGO_status))

subset_HGNCzebrafish_human <- subset_HGNCzebrafish_human %>% 
  left_join(SGO_status[, c("zebrafish_gene_id", "SGO_status")], by = "zebrafish_gene_id") %>%
  mutate(SGO_status = ifelse(is.na(SGO_status), "Not SGO", SGO_status))

# Drop rows where relation_type is "Many to Many"
filteredHGNChuman_zebrafish <- HGNChuman_zebrafish[relation_type != "One to Many"]
filteredHGNCzebrafish_human <- HGNCzebrafish_human[relation_type != "Many to One"]

# Append data
filteredHGNChuman_zebrafish <- rbind(filteredHGNChuman_zebrafish, subset_HGNChuman_zebrafish)
filteredHGNCzebrafish_human <- rbind(filteredHGNCzebrafish_human, subset_HGNCzebrafish_human)

# Average orthologs
average_one_to_many <- filteredHGNChuman_zebrafish %>%
  filter(relation_type == "One to Many") %>%        # Filter to keep only 'One to Many' rows
  summarise(average_count = mean(zebrafish_count),
            median_count = median(zebrafish_count))  

average_one_to_many <- filteredHGNCzebrafish_human %>%
  filter(relation_type == "Many to One") %>%        # Filter to keep only 'Many to One' rows
  summarise(average_count = mean(human_count),
            median_count = median(human_count))    

print(average_one_to_many) # 15.7 (2) zebrafish/25.2 (4) humans

# # Store as .tsv (August 5, 2024)
# write.table(filteredHGNChuman_zebrafish, file = "filteredHGNChuman_zebrafishOrthologs-final.tsv",
#             sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
# 
# write.table(filteredHGNCzebrafish_human, file = "filteredHGNCzebrafish_humanOrthologs-final.tsv",
#             sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

# Visualize HGNC Data ---------------------------------------------------------

# Bar plot
colours <- c("One to One" = "#EE5496", "One to Many" = "#FFC000", "Many to One" = "#78DAD5")
lighter_colours <- alpha(colours, 0.7)

# Convert back to dataframe if throwing errors
# HGNChuman_zebrafish <- as.data.frame(HGNChuman_zebrafish)

## Human to Zebrafish Orthologs including SGO status
plot<-ggplot(HGNChuman_zebrafish, aes(x = relation_type, fill = relation_type)) +
  geom_bar(stat = "count", position = position_dodge(), color = colours[2:1], size = 1, alpha = 0.5) +
  geom_text(stat = 'count', aes(label = ..count..), 
            vjust = -0.5, position = position_dodge(width = 0.9), size = 7) +
  geom_bar(data = subset(HGNChuman_zebrafish, SGO_status == "True SGO"), 
           aes(x = relation_type, fill = relation_type),
           stat = "count", position = position_dodge(), color = colours[1], size = 1) +
  geom_text(data = subset(HGNChuman_zebrafish, SGO_status == "True SGO"), aes(label = ..count..), 
            stat = 'count', vjust = -0.5, position = position_dodge(width = 0.9), size = 7) +
  labs(title = "Human to Zebrafish Orthology Counts By Type",
       x = NULL, y = "Count") +
  theme_minimal() +
  theme(plot.title = element_text(size = 20, hjust = 0.5, 
                                  margin = margin(t = 0, r = 0, b = 10, l = 0)),
        axis.title.y = element_text(size = 18, margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 15),
        legend.position = "none") +
  scale_fill_manual(values = lighter_colours) +
  scale_color_manual(values = colours[2:1]) +  # Apply darker colours for outlines 
  scale_y_continuous(
    limits = c(0, 17500),       
    breaks = seq(0, 17500, by = 2500)
  )

ggsave("HCOP HGNC Symbol Human-Zebrafish.png", plot = plot, device = "png", 
       width = 9.50, height = 8.19, dpi = 600)

## Zebrafish to Human Orthologs including SGO status
plot<-ggplot(HGNCzebrafish_human, aes(x = relation_type, fill = relation_type)) +
  geom_bar(stat = "count", position = position_dodge(), color = colours[c(3, 1)], size = 1, alpha = 0.5) +
  geom_text(stat = 'count', aes(label = ..count..), 
            vjust = -0.5, position = position_dodge(width = 0.9), size = 7) +
  geom_bar(data = subset(HGNCzebrafish_human, SGO_status == "True SGO"), 
           aes(x = relation_type, fill = relation_type),
           stat = "count", position = position_dodge(), color = colours[1], size = 1) +
  geom_text(data = subset(HGNCzebrafish_human, SGO_status == "True SGO"), aes(label = ..count..), 
            stat = 'count', vjust = -0.5, position = position_dodge(width = 0.9), size = 7) +
  labs(title = "Zebrafish to Human Orthology Counts By Type",
       x = NULL, y = "Count") +
  theme_minimal() +
  theme(plot.title = element_text(size = 20, hjust = 0.5, 
                                  margin = margin(t = 0, r = 0, b = 10, l = 0)),
        axis.title.y = element_text(size = 18, margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 15),
        legend.position = "none") +
  scale_fill_manual(values = lighter_colours) +
  scale_color_manual(values = colours[c(3, 1)]) +  # Apply darker colours for outlines 
  scale_y_continuous(
    limits = c(0, 17500),       
    breaks = seq(0, 17500, by = 2500)
  )

ggsave("HCOP HGNC Zebrafish-Human.png", plot = plot, device = "png", 
       width = 9.50, height = 8.19, dpi = 600)

# Boxplot
# Extract the composite_score columns
zebraOrtho_score <- HGNChuman_zebrafish %>% select(composite_score)
humanOrtho_score <- HGNCzebrafish_human %>% select(composite_score)

# Determine the maximum length
max_length <- max(nrow(zebraOrtho_score), nrow(humanOrtho_score))

# Ensure both columns are the same length by filling with NA
zebraOrtho_score <- zebraOrtho_score %>% slice(1:max_length)
humanOrtho_score <- humanOrtho_score %>% slice(1:max_length)

# Combine the composite_score columns into a new dataframe
composite_scores <- data.frame(
  zebraOrtho_score = c(zebraOrtho_score$composite_score, rep(NA, max_length - nrow(zebraOrtho_score))),
  humanOrtho_score = c(humanOrtho_score$composite_score, rep(NA, max_length - nrow(humanOrtho_score)))
)

# Ensure the columns are numeric
composite_scores$zebraOrtho_score <- as.numeric(composite_scores$zebraOrtho_score)
composite_scores$humanOrtho_score <- as.numeric(composite_scores$humanOrtho_score)

# Convert to long format for plot
composite_scores <- composite_scores %>%
  pivot_longer(cols = everything(), names_to = "group", values_to = "score") %>%
  filter(!is.na(score))  # Remove NA values

# Sample 250 points from each group
set.seed(321) # For reproducibility
samples_compositeScores <- composite_scores %>%
  group_by(group) %>%
  sample_n(250)

ggplot(samples_compositeScores, aes(x = group, y = score, fill = group)) +
  stat_halfeye(adjust = 0.4, justification = -0.2, .width = 0, point_color = NA) +
  geom_boxplot(width = 0.25, outlier.color = NA, alpha = 0.5) +
  geom_point(aes(color = group), fill = "white", shape = 21, stroke = .4, size = 2, 
             position = position_jitter(seed = 1, width = .12)) +
  geom_point(aes(fill = group), color = "transparent", shape = 21, stroke = .4, size = 2, 
             alpha = .3, position = position_jitter(seed = 1, width = .12)) +
  stat_summary(geom = "text", fun = "median", aes(label = round(..y.., 2), color = group), 
               family = "Roboto Mono", fontface = "bold", size = 3, vjust = -3.5) +
  stat_summary(geom = "text", fun.data = function(x) return(c(y = max(x) + .025, label = length(x))), 
               aes(label = paste("n =", ..label..), color = group), 
               family = "Roboto Condensed", size = 3, hjust = 0) +
  theme_minimal() +
  labs(x = "Group", y = "Composite Score") +
  theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(), 
        axis.ticks = element_blank(), axis.text.x = element_text(family = "Roboto Mono", size = 10), 
        axis.text.y = element_text(size = 12), 
        axis.title.x = element_text(margin = margin(t = 10), size = 12), 
        plot.subtitle = element_text(color = "grey40", hjust = 0, margin = margin(0, 0, 20, 0)), 
        plot.title.position = "plot", 
        plot.caption = element_text(color = "grey40", lineheight = 1.2, margin = margin(20, 0, 0, 0)), 
        plot.margin = margin(15, 15, 10, 15), 
        legend.position = "none") +
  scale_color_manual(values = c("#EE5496", "#FFC000")) +
  scale_fill_manual(values = c("#EE5496", "#FFC000"))

# WP2446 (RB) Pathway Analysis ------------------------------------------------
# 89 genes in RB pathway (y-axis)
WP2446_humanSymbol <- unique(WP2446_pathwayGenes$HGNC)

subset_WP2446 <- c(
  "RB1", "TP53", "CDK4", "CDK6", "CCND1", "CCNE1", "CDK2", "E2F1", "CCNA2",
  "CDKN1A", "CDKN1B", "MYC", "MCM3", "MCM4", "MCM6", "MCM7", "PCNA", "CDT1",
  "SMARCA2", "CHEK1"
)

# Create a data frame for WP2446_humanSymbol
# RBhuman_genes <- data.frame(
#   human_symbol = WP2446_humanSymbol,
#   stringsAsFactors = FALSE
# )

RBhuman_genes <- data.frame(
  human_symbol = WP2446_humanSymbol[WP2446_humanSymbol %in% subset_WP2446],
  stringsAsFactors = FALSE
)

# Farrell data
matches <- sum(farrellHuman_zebrafish$human_symbol %in% RBhuman_genes$human_symbol)

farrell_HGNChuman_zebrafish <- HGNChuman_zebrafish %>%
  semi_join(farrellHuman_zebrafish, by = c("human_symbol", "zebrafish_symbol"))

# Need number of databases (%) supporting orthology assertions (take average for "Many to Many" relationships)
# Need zebrafish_count - i.e., number of zebrafish orthologs for each respective human gene in pathway
prepare_data <- function(data){
  filtered_data <- RBhuman_genes %>%
    left_join(data, by = "human_symbol") %>%
    mutate(
      human_symbol = factor(human_symbol, levels = WP2446_humanSymbol), # Maintain order
      composite_score = ifelse(is.na(composite_score), 0, composite_score),
      # compositeScore_percent = ifelse(is.na(composite_score), 0, composite_score / 12 * 100),
      zebrafish_count = ifelse(is.na(zebrafish_count), 0, zebrafish_count),
      size_category = case_when(
        zebrafish_count == 0 ~ "0",
        zebrafish_count == 1 ~ "1",
        zebrafish_count == 2 ~ "2",
        zebrafish_count == 3 ~ "3",
        zebrafish_count == 4 ~ "4",
        zebrafish_count == 5 ~ "5",
        zebrafish_count >= 6 ~ "6 or more"
      ),
      # color_group = ifelse(compositeScore_percent == 0, "zero", as.character(seq_along(human_symbol)))
      color_group = ifelse(composite_score == 0, "zero", as.character(seq_along(human_symbol)))
    ) %>%
    mutate(
      size_category = factor(size_category, levels = c("0", "1", "2", "3", "4", "5", "6 or more"))
    )
  return(filtered_data)
}

# No filtering, gene = generic ID
filtered_data1 <- prepare_data(HGNChuman_zebrafish)

# Standard filtering (only filtered "one-to-many" orthologs)
filtered_data2 <- prepare_data(filteredHGNChuman_zebrafish)

# Farrell ortholog genes
filtered_data3 <- prepare_data(farrell_HGNChuman_zebrafish)

# Arbitrarily drop any row where composite score < 5
filtered_data4 <- HGNChuman_zebrafish %>%
  mutate(composite_score = ifelse(composite_score > 5, composite_score, 0),
         zebrafish_count = ifelse(composite_score > 5, zebrafish_count, 0))

filtered_data4 <- prepare_data(filtered_data4)

# Strip plot data
prepare_stripData <- function(data){
  filtered_data <- RBhuman_genes %>%
    left_join(data, by = "human_symbol") %>%
    mutate(
      human_symbol = factor(human_symbol, levels = WP2446_humanSymbol), # Maintain order
      support_count = ifelse(is.na(support_count), 0, support_count)
    )
  return(filtered_data)
}

strip_plotAll <- prepare_stripData(hgnc_orthologsDF)

# Pathway Visualization -------------------------------------------------------

# Define the color palette
palette_colors <- c("#886889", "#B46D78", "#F89771", "#FCBF68", "#FAE079", "#ECDB66", "#E0D567", "#9EAB7C", "#81908B", "#7999B2")
palette_colors <- c("#0070C0", "#B0DE7E")

# Generate 89 colors from the given palette
color_palette <- colorRampPalette(palette_colors)(89)
names(color_palette) <- as.character(1:89)

color_palette <- colorRampPalette(palette_colors)(20)
names(color_palette) <- as.character(1:20)
color_palette <- c(color_palette, "zero" = "#9400D3") 

# Ensure the color palette matches the levels of color_group
color_levels <- c(as.character(1:89), "zero")
color_levels <- c(as.character(1:20), "zero")

# Gray color for the size category legend
gray_color <- "#252525"

# Create the bubble plot (shape = 21 is a filled circle with border)
create_bubblePlot <- function(data){
  median_score <- median(data$composite_score, na.rm = TRUE)
  ggplot(data, aes(x = human_symbol, y = composite_score)) +
    geom_hline(yintercept = median_score, color = "#797979", linetype = "dashed") +
    geom_point(aes(size = size_category, fill = color_group, color = color_group),
               shape = 21, stroke = 1, alpha = 0.6) +
    scale_size_manual(values = c("0" = 2, "1" = 2, "2" = 4,
                                 "3" = 6, "4" = 8, "5" = 10, "6 or more" = 12),
                      drop = FALSE) +
    scale_fill_manual(values = color_palette, breaks = color_levels) +
    scale_color_manual(values = color_palette, breaks = color_levels) +
    scale_y_continuous(limits = c(0, 12), breaks = seq(0, 12, by = 2)) +
    labs(x = "Human Gene Symbol",
         y = "Composite Score Percent (%)",
         size = "Number of Zebrafish Orthologs") +
    theme_minimal() +
    theme(axis.text.x = element_text(size = 8, angle = 90, vjust = 0.2, hjust = 1),
          axis.text.y = element_text(size = 10),
          plot.margin = margin(25, 25, 10, 25),
          legend.position = "top",
          legend.justification = "left",
          legend.title = element_text(size = 10),
          legend.text = element_text(size = 8),
          legend.key.size = unit(0.5, "lines"),
          # Remove axis titles for combined plot
          axis.title = element_blank()) +
    guides(color = "none", fill = "none",
           size = guide_legend(override.aes = list(fill = gray_color), nrow = 1)) +
    coord_cartesian(clip = "off")
}

# No legend
bubble_noLegend <- function(data){
  median_score <- median(data$composite_score, na.rm = TRUE)
  ggplot(data, aes(x = human_symbol, y = composite_score)) +
    geom_hline(yintercept = median_score, color = "#797979", linetype = "dashed") +
    geom_point(aes(size = size_category, fill = color_group, color = color_group),
               shape = 21, stroke = 1, alpha = 0.6) +
    scale_size_manual(values = c("0" = 2, "1" = 2, "2" = 4,
                                 "3" = 6, "4" = 8, "5" = 10, "6 or more" = 12),
                      drop = FALSE) +
    scale_fill_manual(values = color_palette, breaks = color_levels) +
    scale_color_manual(values = color_palette, breaks = color_levels) +
    scale_y_continuous(limits = c(0, 12), breaks = seq(0, 12, by = 2)) +
    labs(x = "Human Gene Symbol",
         y = "Composite Score Percent (%)",
         size = "Number of Zebrafish Orthologs") +
    theme_minimal() +
    theme(axis.text.x = element_text(size = 8, angle = 90, vjust = 0.2, hjust = 1),
          axis.text.y = element_text(size = 10),
          plot.margin = margin(25, 25, 10, 25),
          legend.position = "none",
          # Remove axis titles for combined plot
          axis.title = element_blank()) +
    guides(color = "none", fill = "none",
           size = guide_legend(override.aes = list(fill = gray_color), nrow = 1)) +
    coord_cartesian(clip = "off")
}

plot1<-create_bubblePlot(filtered_data1)
plot2<-bubble_noLegend(filtered_data2)
# plot2a<-create_bubblePlot(filtered_data2)
# plot3<-bubble_noLegend(filtered_data3)
# plot4<-bubble_noLegend(filtered_data4)

# Create a common y-axis title plot
y_title <- ggplot() + 
  annotate("text", x = 1, y = 0.5, label = "Supporting Databases", 
           angle = 90, size = 4) +
  theme_void() +
  coord_cartesian(clip = "off")

# Create a common x-axis title plot
x_title <- ggplot() + 
  annotate("text", x = 0.5, y = 0, label = "Human Gene Symbol", size = 4, vjust = 0.5) +
  theme_void() +
  coord_cartesian(clip = "off")

# Combine the plots with common axis titles
combined_plot <- (y_title | wrap_plots(plot1, plot2, ncol = 1)) +
  plot_layout(widths = c(0.05, 1))

# combined_plot <- (y_title | wrap_plots(plot2a, plot4, plot3, ncol = 1)) +
#   plot_layout(widths = c(0.05, 1))

final_plot <- combined_plot / x_title +
  plot_layout(heights = c(1, 0.05))

# Print the final plot
print(final_plot)

#6.96
ggsave("Overall Confidence in Orthology per Gene.png", plot = final_plot, device = "png", 
       width = 13.34, height = 7.5, dpi = 900)

# Pathway Visualization (Poster Size) ------------------------------------------
create_bubblePlot <- function(data){
  median_score <- median(data$composite_score, na.rm = TRUE)
  ggplot(data, aes(x = human_symbol, y = composite_score)) +
    geom_hline(yintercept = median_score, color = "#797979", linetype = "dashed") +
    geom_point(aes(size = size_category, fill = color_group, color = color_group), 
               shape = 21, stroke = 1, alpha = 0.6) +
    scale_size_manual(values = c("0" = 4, "1" = 4, "2" = 6,
                                 "3" = 8, "4" = 10, "5" = 12, "6 or more" = 14), 
                      drop = FALSE) +
    scale_fill_manual(values = color_palette, breaks = color_levels) +
    scale_color_manual(values = color_palette, breaks = color_levels) +
    scale_y_continuous(limits = c(0, 12), breaks = seq(0, 12, by = 2)) +
    labs(x = NULL, y = NULL, size = "Number of Zebrafish Orthologs") +
    theme_minimal() +
    theme(axis.text.x = element_text(size = 12, angle = 45, vjust = 1, hjust = 1), 
          axis.text.y = element_text(size = 10),
          plot.margin = margin(5, 5, 5, 5),
          legend.position = "top",
          legend.justification = "left",
          legend.title = element_text(size = 10),
          legend.text = element_text(size = 8),
          legend.key.size = unit(0.5, "lines")) +
    guides(color = "none", fill = "none", 
           size = guide_legend(override.aes = list(fill = gray_color), nrow = 1)) +
    coord_cartesian(clip = "off")
}

# No legend
bubble_noLegend <- function(data){
  median_score <- median(data$composite_score, na.rm = TRUE)
  ggplot(data, aes(x = human_symbol, y = composite_score)) +
    geom_hline(yintercept = median_score, color = "#797979", linetype = "dashed") +
    geom_point(aes(size = size_category, fill = color_group, color = color_group), 
               shape = 21, stroke = 1, alpha = 0.6) +
    scale_size_manual(values = c("0" = 4, "1" = 4, "2" = 6,
                                 "3" = 8, "4" = 10, "5" = 12, "6 or more" = 14), 
                      drop = FALSE) +
    scale_fill_manual(values = color_palette, breaks = color_levels) +
    scale_color_manual(values = color_palette, breaks = color_levels) +
    scale_y_continuous(limits = c(0, 12), breaks = seq(0, 12, by = 2)) +
    labs(x = NULL, y = NULL, size = "Number of Zebrafish Orthologs") +
    theme_minimal() +
    theme(axis.text.x = element_text(size = 12, angle = 45, vjust = 1, hjust = 1), 
          axis.text.y = element_text(size = 10),
          plot.margin = margin(5, 5, 5, 5),
          legend.position = "none") +
    guides(color = "none", fill = "none", 
           size = guide_legend(override.aes = list(fill = gray_color), nrow = 1)) +
    coord_cartesian(clip = "off")
}

plot1 <- create_bubblePlot(filtered_data1)
plot2 <- bubble_noLegend(filtered_data2)

# Combine the plots without axis titles
combined_plot <- wrap_plots(plot1, plot2, ncol = 1) +
  plot_layout(widths = c(1))

final_plot <- combined_plot +
  plot_layout(heights = c(1, 1))

# Print the final plot
print(final_plot)

# Save the plot with the specified size in cm
ggsave("Overall_Confidence_in_Orthology_per_Gene_x3.png", plot = final_plot, device = "png", 
       width = 28, height = 13.03, units = "cm", dpi = 900)

# Rotated Bubble Plot ----------------------------------------------------------
create_bubblePlot <- function(data){
  median_score <- median(data$composite_score, na.rm = TRUE)
  ggplot(data, aes(y = human_symbol, x = composite_score)) +
    geom_vline(xintercept = median_score, color = "#797979", linetype = "dashed") +
    geom_point(aes(size = size_category, fill = color_group, color = color_group), 
               shape = 21, stroke = 1, alpha = 0.6) +
    scale_size_manual(values = c("0" = 4, "1" = 4, "2" = 6,
                                 "3" = 8, "4" = 10, "5" = 12, "6 or more" = 14), 
                      drop = FALSE) +
    scale_fill_manual(values = color_palette, breaks = color_levels) +
    scale_color_manual(values = color_palette, breaks = color_levels) +
    scale_x_continuous(limits = c(0, 12), breaks = seq(0, 12, by = 2)) +
    labs(x = "Number of Zebrafish Orthologs", y = NULL, size = "Number of Zebrafish Orthologs") +
    theme_minimal() +
    theme(axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 10, angle = 45, vjust = 1, hjust = 1),
          plot.margin = margin(5, 5, 5, 5),
          legend.position = "top",
          legend.justification = "left",
          legend.title = element_text(size = 10),
          legend.text = element_text(size = 8),
          legend.key.size = unit(0.5, "lines")) +
    guides(color = "none", fill = "none", 
           size = guide_legend(override.aes = list(fill = gray_color), nrow = 1)) +
    coord_cartesian(clip = "off")
}

bubble_noLegend <- function(data){
  median_score <- median(data$composite_score, na.rm = TRUE)
  ggplot(data, aes(y = human_symbol, x = composite_score)) +
    geom_vline(xintercept = median_score, color = "#797979", linetype = "dashed") +
    geom_point(aes(size = size_category, fill = color_group, color = color_group), 
               shape = 21, stroke = 1, alpha = 0.6) +
    scale_size_manual(values = c("0" = 4, "1" = 4, "2" = 6,
                                 "3" = 8, "4" = 10, "5" = 12, "6 or more" = 14), 
                      drop = FALSE) +
    scale_fill_manual(values = color_palette, breaks = color_levels) +
    scale_color_manual(values = color_palette, breaks = color_levels) +
    scale_x_continuous(limits = c(0, 12), breaks = seq(0, 12, by = 2)) +
    labs(x = "Number of Zebrafish Orthologs", y = NULL, size = "Number of Zebrafish Orthologs") +
    theme_minimal() +
    theme(axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 10, angle = 45, vjust = 1, hjust = 1),
          plot.margin = margin(5, 5, 5, 5),
          legend.position = "none") +
    guides(color = "none", fill = "none", 
           size = guide_legend(override.aes = list(fill = gray_color), nrow = 1)) +
    coord_cartesian(clip = "off")
}

plot1 <- create_bubblePlot(filtered_data1)
plot2 <- bubble_noLegend(filtered_data2)

# Combine the plots with axis titles and 2 columns
combined_plot <- wrap_plots(plot1, plot2, nrow = 1) +
  plot_layout(heights = c(1))

final_plot <- combined_plot +
  plot_layout(widths = c(1, 1))

# Print the final plot
print(final_plot)

# Save the plot with the specified size in cm
ggsave("Column Overall Confidence in Orthology per Gene.png", plot = final_plot, device = "png", 
       width = 28, height = 13.03, units = "cm", dpi = 900)

# Strip plot poster ------------------------------------------------------------
g <- ggplot(strip_plotAll, aes(x = support_count, y = human_symbol, color = human_symbol)) +
  coord_flip() +
  scale_y_discrete(limits = rev(levels(strip_plotAll$human_symbol))) + # Ensures the y-axis labels are ordered as per WP2446_humanSymbol
  scale_x_continuous(limits = c(0, max(strip_plotAll$support_count, na.rm = TRUE)), expand = c(0.02, 0.02)) +
  labs(x = "Number of Supporting Databases", y = "Human Gene Symbol") +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 16),
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 12),
    panel.grid = element_blank()
  ) +
  geom_hline(aes(yintercept = median(support_count, na.rm = TRUE)), color = "gray70", size = 0.6) +
  stat_summary(fun = mean, geom = "point", size = 5) +
  geom_jitter(size = 2, alpha = 0.25, width = 0.2)

# Print the plot
print(g)

# Internal Validation Work -----------------------------------------------------
# Liming Email 07/09/2024
# Ensembl v109 BioMart Export (Zebrafish Gene IDs)
zebrafish_genesMart <- read_delim("zebrafish-genes-mart-export.txt", delim = "\t", escape_double = FALSE, 
                                  trim_ws = TRUE)

zebrafish_genesMartDF <- as.data.frame(zebrafish_genesMart)


human_genesMart <- read_delim("human_ensembl-entrez-mart_export.txt", delim = "\t", escape_double = FALSE,
                              trim_ws = TRUE)

human_genesMartDF <- as.data.frame(human_genesMart)

# Download the zebrafish id conversion table from ZFIN database and validate 
# the Ensembl ID for each gene based on ZFIN ID
# ZFIN-Ensembl 1:1 mapping
zfin_ensemblGene <- read_delim("zfin-ensembl_1_to_1_2024.07.09.txt", 
                               delim = "\t", escape_double = FALSE, trim_ws = TRUE)

zfin_ensemblGeneDF <- as.data.frame(zfin_ensemblGene)
zfin_ensemblGeneDF$...5 <- NULL

zfin_ensemblConflict <- zfin_ensemblGeneDF %>%
  group_by(`ZFIN ID`) %>%
  summarise(across(everything(), ~ paste(unique(.), collapse = ", ")),
            .groups = "drop")

zfin_entrezGene <- read_delim("zfin-entrezgene_2024.07.09.txt", delim = "\t", 
                              escape_double = FALSE, trim_ws = TRUE)

zfin_entrezGeneDF <- as.data.frame(zfin_entrezGene)
zfin_entrezGeneDF$...5 <- NULL

zfin_entrezConflict <- zfin_entrezGeneDF %>%
  group_by(`ZFIN ID`) %>%
  summarise(across(everything(), ~ paste(unique(na.omit(.)), collapse = ", ")),
            .groups = "drop")

# Add the new Ensembl ID to your table with cluster number and record how  
# many conflicts has been found
# 108 ZFIN IDs with > 1 corresponding Ensembl ID
ensembl_conflict <- zfin_ensemblConflict %>%
  filter(str_count(`Ensembl ID`, ",") >= 1)

# 2 ZFIN IDs with > 2 corresponding Entrez ID
zfin_conflict <- zfin_entrezConflict %>%
  filter(str_count(`NCBI Gene ID`, ",") >= 1)

# Rename columns for merge
zfin_ensemblConflict <- zfin_ensemblConflict %>%
  dplyr::rename(zfin_id = `ZFIN ID`, zfin_ensembl_gene = `Ensembl ID`)

# Merge
hgnc_orthologsDF_copy <- copy(hgnc_orthologsDF)

hgnc_orthologsDF <- hgnc_orthologsDF %>%
  left_join(zfin_ensemblConflict %>% select(zfin_id, zfin_ensembl_gene), by = "zfin_id")

hgnc_orthologsDF <- hgnc_orthologsDF %>%
  relocate(row_id, .before = human_entrez_gene)

hgnc_orthologsDF <- hgnc_orthologsDF %>%
  relocate(zfin_ensembl_gene, .after = zebrafish_ensembl_gene)

# Compare Ensembl IDs as listed in ZFIN to Ensembl IDs from HCOP
compare_ensemblID <- function(col1, col2) {
  if (is.na(col1) | is.na(col2)) {
    return(NA)
  }
  
  values1 <- unlist(str_split(col1, ","))
  values2 <- unlist(str_split(col2, ","))
  
  if (length(setdiff(values1, values2)) == 0 & length(setdiff(values2, values1)) == 0) {
    return("Ensembl Match")
  } else {
    # Mismatches in both directions
    # E.g., if genes1 = ["A", "B", "C"] and genes2 = ["B", "C", "D"]
    # Returns A + D (2 mismatches)
    mismatches <- length(setdiff(values1, values2)) + length(setdiff(values2, values1))
    return(paste(mismatches, "mismatches"))
    # return(paste(mismatches))
  }
}

# Apply function to each row to create the zfin_ensembl_conflict column
hgnc_orthologsDF <- hgnc_orthologsDF %>%
  rowwise() %>%
  mutate(zfin_ensembl_conflict = compare_ensemblID(zebrafish_ensembl_gene, zfin_ensembl_gene)) %>%
  ungroup()

subset_count <- hgnc_orthologsDF %>%
  # Pull mismatch count
  filter(zfin_ensembl_conflict != "Ensembl Match") %>%
  mutate(mismatches = as.numeric(str_extract(zfin_ensembl_conflict, "\\d+"))) %>%
  pull(mismatches)

mean(subset_count) # 1 mismatches

# 42/132,195 rows do not match
hgnc_orthologsDF %>% filter(str_detect(zfin_ensembl_conflict, "mismatches")) %>% nrow()

# 109,164/132,195 rows match
hgnc_orthologsDF %>% filter(zfin_ensembl_conflict == "Ensembl Match") %>% nrow()

# 22,989/132,195 NA (16.1%)
hgnc_orthologsDF %>% filter(is.na(zfin_ensembl_conflict)) %>% nrow()

# Check Ensembl ID against Entrez ID for zebrafish
entrez_ensemblConflict <- zebrafish_genesMartDF %>%
  group_by(`NCBI gene (formerly Entrezgene) ID`) %>%
  summarise(across(everything(), ~ paste(unique(na.omit(.)), collapse = ", ")),
            .groups = "drop")

# 3700 Entrez IDs with > 1 corresponding Ensembl ID
ensembl_conflict <- entrez_ensemblConflict %>%
  filter(str_count(`Gene stable ID`, ",") >= 1)

# 236 Entrez IDs with > 1 corresponding ZFIN ID
zfin_conflict <- entrez_ensemblConflict %>%
  filter(str_count(`ZFIN ID`, ",") >= 1)

# Rename columns for merge
entrez_ensemblConflict <- entrez_ensemblConflict %>%
  dplyr::rename(zebrafish_entrez_gene = `NCBI gene (formerly Entrezgene) ID`, 
                entrez_ensembl_gene = `Gene stable ID`)

# Merge
hgnc_orthologsDF_copy <- copy(hgnc_orthologsDF)

entrez_ensemblConflict <- entrez_ensemblConflict %>%
  mutate(zebrafish_entrez_gene = as.character(zebrafish_entrez_gene))

hgnc_orthologsDF <- hgnc_orthologsDF %>%
  left_join(entrez_ensemblConflict %>% select(zebrafish_entrez_gene, entrez_ensembl_gene), 
            by = "zebrafish_entrez_gene")

hgnc_orthologsDF <- hgnc_orthologsDF %>%
  relocate(entrez_ensembl_gene, .after = zfin_ensembl_gene)

# Compare Ensembl ID from NCBI to Ensembl IDs from HCOP
hgnc_orthologsDF <- hgnc_orthologsDF %>%
  rowwise() %>%
  mutate(entrez_ensembl_conflict = compare_ensemblID(zebrafish_ensembl_gene, entrez_ensembl_gene)) %>%
  ungroup()

subset_count <- hgnc_orthologsDF %>%
  # Pull mismatch count
  filter(entrez_ensembl_conflict != "Ensembl Match") %>%
  mutate(mismatches = as.numeric(str_extract(entrez_ensembl_conflict, "\\d+"))) %>%
  pull(mismatches)

mean(subset_count) # 1 mismatches

# 34/132,195 rows do not match
hgnc_orthologsDF %>% filter(str_detect(entrez_ensembl_conflict, "mismatches")) %>% nrow()

# 110,051/132,195 rows match
hgnc_orthologsDF %>% filter(entrez_ensembl_conflict == "Ensembl Match") %>% nrow()

# 22,110/132,195 NA
hgnc_orthologsDF %>% filter(is.na(entrez_ensembl_conflict)) %>% nrow()

hgnc_orthologsDF <- hgnc_orthologsDF %>%
  relocate(zfin_ensembl_conflict, .after = zfin_ensembl_gene)

hgnc_orthologsDF <- hgnc_orthologsDF %>%
  relocate(entrez_ensembl_conflict, .after = entrez_ensembl_gene)


# Check Human Ensembl ID against Entrez ID
humanEntrez_ensemblConflict <- human_genesMartDF %>%
  group_by(`NCBI gene (formerly Entrezgene) ID`) %>%
  summarise(across(everything(), ~ paste(unique(na.omit(.)), collapse = ", ")),
            .groups = "drop")

# 2134 Entrez IDs with > 1 corresponding Ensembl ID
ensembl_conflict <- humanEntrez_ensemblConflict %>%
  filter(str_count(`Gene stable ID`, ",") >= 1)

# Rename columns for merge
humanEntrez_ensemblConflict <- humanEntrez_ensemblConflict %>%
  dplyr::rename(human_entrez_gene = `NCBI gene (formerly Entrezgene) ID`, 
                humanEntrez_ensembl_gene = `Gene stable ID`)

# Merge
hgnc_orthologsDF_copy <- copy(hgnc_orthologsDF)

humanEntrez_ensemblConflict <- humanEntrez_ensemblConflict %>%
  mutate(human_entrez_gene = as.character(human_entrez_gene))

hgnc_orthologsDF <- hgnc_orthologsDF %>%
  left_join(humanEntrez_ensemblConflict %>% select(human_entrez_gene, humanEntrez_ensembl_gene), 
            by = "human_entrez_gene")

hgnc_orthologsDF <- hgnc_orthologsDF %>%
  relocate(humanEntrez_ensembl_gene, .after = human_ensembl_gene)

# Compare Ensembl ID from NCBI to Ensembl IDs from HCOP
hgnc_orthologsDF <- hgnc_orthologsDF %>%
  rowwise() %>%
  mutate(human_ensembl_conflict = compare_ensemblID(human_ensembl_gene, humanEntrez_ensembl_gene)) %>%
  ungroup()

subset_count <- hgnc_orthologsDF %>%
  # Pull mismatch count
  filter(human_ensembl_conflict != "Ensembl Match") %>%
  mutate(mismatches = as.numeric(str_extract(human_ensembl_conflict, "\\d+"))) %>%
  pull(mismatches)

mean(subset_count) # 1 mismatches

# 9/132,195 rows do not match
hgnc_orthologsDF %>% filter(str_detect(human_ensembl_conflict, "mismatches")) %>% nrow()

# 117,295/132,195 rows match
hgnc_orthologsDF %>% filter(human_ensembl_conflict == "Ensembl Match") %>% nrow()

# 14,891/132,195 NA
hgnc_orthologsDF %>% filter(is.na(human_ensembl_conflict)) %>% nrow()

hgnc_orthologsDF <- hgnc_orthologsDF %>%
  relocate(human_ensembl_conflict, .after = humanEntrez_ensembl_gene)


write.table(hgnc_orthologsDF, file = "hgnc_ensemblConflict_update.tsv",
            sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

# Liming email: "Ensembl has some additional genes from scaffold genome, so if you 
# have multiple matches for one gene, only keep the one that exists in your original table"
# Function to filter entries within each cell
filter_entries <- function(gene_list, reference_list) {
  if (is.na(gene_list) | is.na(reference_list)) {
    return(gene_list)
  }
  
  genes <- unlist(str_split(gene_list, ","))
  reference_genes <- unlist(str_split(reference_list, ","))
  
  filtered_genes <- genes[genes %in% reference_genes]
  
  if (length(filtered_genes) == 0) {
    return(NA)
  } else {
    return(paste(filtered_genes, collapse = ","))
  }
}

hgnc_orthologsDF <- hgnc_orthologsDF %>%
  rowwise() %>%
  mutate(
    humanEntrez_ensembl_gene = filter_entries(humanEntrez_ensembl_gene, human_ensembl_gene),
    zfin_ensembl_gene = filter_entries(zfin_ensembl_gene, zebrafish_ensembl_gene),
    entrez_ensembl_gene = filter_entries(entrez_ensembl_gene, zebrafish_ensembl_gene)
  ) %>%
  ungroup()

# Recalculate mismatches
hgnc_orthologsDF <- hgnc_orthologsDF %>%
  rowwise() %>%
  mutate(zfin_ensembl_conflict = compare_ensemblID(zebrafish_ensembl_gene, zfin_ensembl_gene)) %>%
  ungroup()

hgnc_orthologsDF <- hgnc_orthologsDF %>%
  rowwise() %>%
  mutate(entrez_ensembl_conflict = compare_ensemblID(zebrafish_ensembl_gene, entrez_ensembl_gene)) %>%
  ungroup()

hgnc_orthologsDF <- hgnc_orthologsDF %>%
  rowwise() %>%
  mutate(human_ensembl_conflict = compare_ensemblID(human_ensembl_gene, humanEntrez_ensembl_gene)) %>%
  ungroup()

# Ensembl ID from ZFIN vs Ensembl ID (using Entrez ID to match)
hgnc_orthologsDF <- hgnc_orthologsDF %>%
  rowwise() %>%
  mutate(zfin_entrezEnsembl_conflict = compare_ensemblID(zfin_ensembl_gene, entrez_ensembl_gene)) %>%
  ungroup()

# Subset of mismatches
subset_df <- hgnc_orthologsDF %>%
  filter(zfin_entrezEnsembl_conflict != "Ensembl Match") %>%
  mutate(mismatches = as.numeric(str_extract(zfin_entrezEnsembl_conflict, "\\d+")))

subset_count <- hgnc_orthologsDF %>%
  # Pull mismatch count
  filter(zfin_entrezEnsembl_conflict != "Ensembl Match") %>%
  mutate(mismatches = as.numeric(str_extract(zfin_entrezEnsembl_conflict, "\\d+"))) %>%
  pull(mismatches)

mean(subset_count) # 1 mismatches

# 10/132,195 rows do not match
hgnc_orthologsDF %>% filter(str_detect(zfin_entrezEnsembl_conflict, "mismatches")) %>% nrow()

# 95,084/132,195 rows match
hgnc_orthologsDF %>% filter(zfin_entrezEnsembl_conflict == "Ensembl Match") %>% nrow()

# 37,101/132,195 NA
hgnc_orthologsDF %>% filter(is.na(zfin_entrezEnsembl_conflict)) %>% nrow()

# ID correction
# Find the most frequent ZF Ensembl ID (not including NA)
consensus_ensemblGene <- function(..., na.rm = TRUE) {
  values <- c(...)
  values <- values[!is.na(values)]
  if (length(values) == 0) {
    return(NA)
  }
  most_frequent <- sort(table(values), decreasing = TRUE)[1]
  names(most_frequent)
}

# Create the consensus_zebrafish_ensembl column
hgnc_orthologsDF <- hgnc_orthologsDF %>%
  rowwise() %>%
  mutate(consensus_zebrafish_ensembl = consensus_ensemblGene(zebrafish_ensembl_gene, zfin_ensembl_gene, entrez_ensembl_gene)) %>%
  ungroup()

hgnc_orthologsDF <- hgnc_orthologsDF %>%
  relocate(consensus_zebrafish_ensembl, .before = zebrafish_ensembl_gene)

# Liming's email: "For each zebrafish gene, it should have a unique Entrez ID + 
# Ensembl ID + ZFIN ID + gene symbol combination."

unique_combinations <- function(df, group_col, check_cols) {
  # Check for unique combinations within each group
  unique_combinations <- df %>%
    group_by(!!sym(group_col)) %>%
    summarise(
      across(all_of(check_cols), ~ n_distinct(.), .names = "unique_{col}"),
      .groups = "drop"
    )
  
  # Filter to find groups with inconsistencies
  inconsistent_groups <- unique_combinations %>%
    filter(if_any(starts_with("unique_"), ~ . > 1))
  
  return(list(unique_combinations = unique_combinations, 
              inconsistent_groups = inconsistent_groups))
}

# Columns to check for unique combinations
group_col <- "zebrafish_symbol"
check_cols <- c("consensus_zebrafish_ensembl", "zfin_id", "zebrafish_entrez_gene")

# Call the function and store the results
unique_results <- unique_combinations(hgnc_orthologsDF, group_col, check_cols)

# Extract the results
unique_groups <- unique_results$unique_combinations
inconsistent_groups <- unique_results$inconsistent_groups

# Store table of mismatched IDs
write.table(inconsistent_groups, file = "zfSymbol_inconsistentID.tsv",
            sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)







# Check for unique combinations within each group of zebrafish genes
unique_combinations <- hgnc_orthologsDF %>%
  group_by(zebrafish_entrez_gene) %>%
  summarise(
    unique_ensembl = n_distinct(consensus_zebrafish_ensembl),
    unique_zfin = n_distinct(zfin_id),
    unique_symbol = n_distinct(zebrafish_symbol),
    .groups = "drop"
  )

# 11 475
sum(hgnc_orthologsDF$zebrafish_entrez_gene == "")

# 22 995
n_distinct(hgnc_orthologsDF$zebrafish_entrez_gene)

# Filter to find groups with inconsistencies
inconsistent_groups <- unique_combinations %>%
  filter(unique_ensembl > 1 | unique_zfin > 1 | unique_symbol > 1)




# # Combine unique entries and ensure consistency for each group
# consistent_hgnc_orthologsDF <- hgnc_orthologsDF %>%
#   group_by(zebrafish_entrez_gene) %>%
#   summarise(
#     consensus_zebrafish_ensembl = paste(unique(consensus_zebrafish_ensembl), collapse = ","),
#     zfin_id = paste(unique(zfin_id), collapse = ","),
#     zebrafish_symbol = paste(unique(zebrafish_symbol), collapse = ","),
#     .groups = "drop"
#   )



# OrthoDB only subset (should be 70732 rows)
hgnc_orthoDB <- subset(hgnc_orthologsDF, support_count == 1 & support == "OrthoDB")

# Opposite
notHgnc_orthoDB <- subset(hgnc_orthologsDF, !(support_count == 1 & support == "OrthoDB"))
n_distinct(notHgnc_orthoDB$human_gene_id) # 18 835 (260 genes less than full data)

not_orthoDB <- unique(notHgnc_orthoDB$zebrafish_gene_id)
only_orthoDB <- as.data.frame(unique(hgnc_orthoDB$zebrafish_gene_id[which(!hgnc_orthoDB$zebrafish_gene_id %in% not_orthoDB)]))
only_orthoDB <- only_orthoDB[[1]]

housekeeping <- subset(farrell_categorizedGenes, class == "Housekeeping")

class_onlyOrtho <- as.data.frame(only_orthoDB[which(only_orthoDB %in% farrell_categorizedGenes$gene)])
class_onlyOrtho <- subset(farrell_categorizedGenes, gene %in% only_orthoDB)


hgnc_housekeeping <- subset(hgnc_orthoDB, zebrafish_gene_id %in% housekeeping$gene)
hgnc_housekeeping <- subset(HGNCzebrafish_human, zebrafish_gene_id %in% housekeeping$gene)

hgnc_housekeeping <- subset(housekeeping, gene %in% only_orthoDB)



zfin <- subset(hgnc_orthologsDF, grepl("ZFIN", support))
# write.table(zfin, file = "HCOP_zfinOrthologs.tsv",
#             sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)


single_support <- subset(hgnc_orthologsDF, support_count == 1 & !support %in% c("NCBI", "ZFIN"))


# Number of intersections per combination
manualSupport_counts <- subset(support_counts, grepl("NCBI|ZFIN", support_combination))

# Create a new column with the number of databases supporting a particular assertion
manualSupport_counts$support_count <- sapply(strsplit(manualSupport_counts$support_combination, "&"), length)

manualSupport_counts <- subset(manualSupport_counts, support_count > 5)
sum(manualSupport_counts$n)

test_counts <- subset(support_counts, !grepl("NCBI|ZFIN", support_combination))

test_counts <- support_counts
test_counts$support_count <- sapply(strsplit(test_counts$support_combination, "&"), length)
test_counts <- subset(test_counts, support_count > 5)
sum(test_counts$n)

appl2 <- subset(hgnc_orthologsDF, zebrafish_gene_id == "appl2") %>%
  mutate(row_id = row_number())

U3 <- subset(hgnc_orthologsDF, zebrafish_gene_id == "U3") %>% 
  mutate(row_id = row_number())

mcat <- subset(hgnc_orthologsDF, grepl("mcat|MCAT", zebrafish_gene_id)) %>%
  mutate(row_id = row_number())

treefam <- subset(hgnc_orthologsDF, grepl("Treefam", support)) %>%
  mutate(row_id = row_number())

SNORA75 <- subset(hgnc_orthologsDF, human_gene_id == "SNORA75") %>%
  mutate(row_id = row_number())

# Group identical (non-unique) ortholog predictions based on gene ID and chromosome of gene
# After grouping identical predictions, remove duplicate assertions made from the same database
# Recount true number of supporting databases
collapsed_dt <- treefam[, .(
  row_id = paste(unique(row_id), collapse = ","),
  human_entrez_gene = paste(setdiff(unique(human_entrez_gene), "-"), collapse = ", "),
  human_ensembl_gene = paste(setdiff(unique(human_ensembl_gene), "-"), collapse = ", "),
  hgnc_id = paste(setdiff(unique(hgnc_id), "-"), collapse = ", "),
  human_name = paste(setdiff(unique(human_name), "-"), collapse = ", "),
  human_symbol = paste(setdiff(unique(human_symbol), "-"), collapse = ", "),
  human_assert_ids = paste(setdiff(unique(human_assert_ids), "-"), collapse = ", "),
  zebrafish_entrez_gene = paste(setdiff(unique(zebrafish_entrez_gene), "-"), collapse = ", "),
  zebrafish_ensembl_gene = paste(setdiff(unique(zebrafish_ensembl_gene), "-"), collapse = ", "),
  zfin_id = paste(setdiff(unique(zfin_id), "-"), collapse = ", "),
  zebrafish_name = paste(setdiff(unique(zebrafish_name), "-"), collapse = ", "),
  zebrafish_symbol = paste(setdiff(unique(zebrafish_symbol), "-"), collapse = ", "),
  zebrafish_assert_ids = paste(setdiff(unique(zebrafish_assert_ids), "-"), collapse = ", "),
  support = paste(unique(support), collapse = ","),
  support_combination = paste(unique(unlist(strsplit(support_combination, "&"))), collapse = "&"),
  support_count = length(unique(unlist(strsplit(support_combination, "&"))))
), by = .(human_gene_id, zebrafish_gene_id, zebrafish_chr, human_chr)]

# Reorder the columns to match the original order
setcolorder(collapsed_dt, colnames(treefam))

# Identify rows with more than one ID (i.e., collapsed rows)
collapsed_rows <- collapsed_dt[sapply(strsplit(row_id, ","), length) > 1]

test <- hgnc_orthologsDF %>%
  mutate(zebrafish_gene_id = tolower(zebrafish_gene_id))

TESTzebrafish_human <- hgnc_orthologsDF %>%
  # 25 158 unique zebrafish gene IDs
  group_by(zebrafish_gene_id) %>%
  summarise(
    # Concatenate all unique zebrafish symbols associated with each zebrafish gene ID, excluding "-"
    zebrafish_symbol = paste(setdiff(unique(zebrafish_symbol), "-"), collapse = ", "),
    zebrafish_entrez_gene = paste(setdiff(unique(zebrafish_entrez_gene), "-"), collapse = ", "),
    zfin_id = paste(setdiff(unique(zfin_id), "-"), collapse = ", "),
    zebrafish_ensembl_gene = paste(setdiff(unique(zebrafish_ensembl_gene), "-"), collapse = ", "),
    human_count = n_distinct(human_gene_id),
    # Concatenate all unique human IDs
    human_gene_id = paste(setdiff(unique(human_gene_id), "-"), collapse = ", "),
    human_symbol = paste(setdiff(unique(human_symbol), "-"), collapse = ", "),
    human_entrez_gene = paste(setdiff(unique(human_entrez_gene), "-"), collapse = ", "),
    human_ensembl_gene = paste(setdiff(unique(human_ensembl_gene), "-"), collapse = ", "),
    average_support = mean(support_count, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  mutate(relation_type = case_when(
    human_count == 1 ~ "One to One",
    human_count > 1 ~ "Many to One",
    TRUE ~ "No Match"
  ))

case_sensitiveGenes <- HGNCzebrafish_human %>%
  mutate(zebrafish_gene_id = tolower(zebrafish_gene_id))

case_sensitiveGenes <- subset(HGNCzebrafish_human,
                              !(zebrafish_gene_id %in% TESTzebrafish_human$zebrafish_gene_id))