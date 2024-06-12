# Human-Zebrafish Orthologs listed in HGNC Comparison of Orthology Predictions (HCOP)
# Note that HCOP (i.e., HGNC Comparison of Orthology Predictions) uses 
# Ensembl v109 as of May 15, 2024
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

# Path to the Excel file
path <- "farrell_categorized-genes-spatiotemporal-variation.xlsx"

# Calling the function and storing the result
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

outpath <- "C:/Users/acale/OneDrive/Documents/Waterloo BME/Co-op/4C/pineoblastoma/Data/Output/"
outDir  <- sprintf("%s/Output", dirname(outpath))
setwd(outDir)

# Read in Farrell Ortholog ist
farrellHuman_zebrafish <- read_delim("farrellHuman_zebrafishOrthologs.tsv",
                                     delim = "\t", escape_double = FALSE, trim_ws = TRUE)

# Quality Check ---------------------------------------------------------------

# HGNC null entries 
no_entries <- sum(hgnc_orthologsDF$zebrafish_entrez_gene == "-")  # 11 645 blank
no_entries <- sum(hgnc_orthologsDF$zebrafish_ensembl_gene == "-") # 2 081  blank
no_entries <- sum(hgnc_orthologsDF$zebrafish_symbol == "-")       # 5 295  blank
no_entries <- sum(hgnc_orthologsDF$zfin_id == "-")                # 16 005 blank

no_entries <- sum(hgnc_orthologsDF$human_symbol == "-")           # 337    blank

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

# Filter HGNC Data ------------------------------------------------------------
# Create a new column with the number of databases supporting a particular assertion
hgnc_orthologsDF$support_count <- sapply(strsplit(hgnc_orthologsDF$support_combination, "&"), length)

# Set to data table for faster operations
setDT(hgnc_orthologsDF)

# Create generic human_gene_id column
hgnc_orthologsDF[, human_gene_id := fifelse(human_symbol != "-", 
                                            human_symbol, 
                                            fifelse(human_entrez_gene != "-", 
                                                    human_entrez_gene, 
                                                    human_ensembl_gene))]

# Check that only 377 entries have Ensembl ID (number of NULL entries for HGNC)
hgnc_orthologsDF %>%
  filter(grepl("ENSG00000", human_gene_id)) %>%
  summarise(total = n())

# Create generic zebrafish_gene_id column
hgnc_orthologsDF[, zebrafish_gene_id := fifelse(zebrafish_symbol != "-", 
                                                zebrafish_symbol, 
                                                fifelse(zebrafish_entrez_gene != "-", 
                                                        zebrafish_entrez_gene, 
                                                        fifelse(zfin_id != "-", 
                                                                zfin_id, 
                                                                zebrafish_ensembl_gene)))]

# Check how many entries have Ensembl ID (4529)
hgnc_orthologsDF %>%
  filter(grepl("ENSDARG000", zebrafish_gene_id)) %>%
  summarise(total = n())

# Create a copy of the original data
hgnc_orthologsDF_original <- copy(hgnc_orthologsDF)

# Create temporary columns with the transformed case and check for duplicates
# Different databases follow different nomenclature (e.g., ABR vs abr, MCAT vs mcat)
# Need to convert to a standard nomenclature so that gene counts are not falsely inflated
hgnc_orthologsDF[, zebrafish_gene_id_lower := tolower(zebrafish_gene_id)]
hgnc_orthologsDF[, human_gene_id_upper := toupper(human_gene_id)]
zebrafish_dupes <- hgnc_orthologsDF[, .N, by = zebrafish_gene_id_lower][N > 1]$zebrafish_gene_id_lower
human_dupes <- hgnc_orthologsDF[, .N, by = human_gene_id_upper][N > 1]$human_gene_id_upper

# Update the original columns based on the duplicates found
hgnc_orthologsDF[, zebrafish_gene_id := ifelse(tolower(zebrafish_gene_id) %in% zebrafish_dupes, tolower(zebrafish_gene_id), zebrafish_gene_id)]
hgnc_orthologsDF[, human_gene_id := ifelse(toupper(human_gene_id) %in% human_dupes, toupper(human_gene_id), human_gene_id)]

# Remove the temporary columns
hgnc_orthologsDF[, c("zebrafish_gene_id_lower", "human_gene_id_upper") := NULL]

hgnc_orthologsDF[grepl("ensdarg0000", zebrafish_gene_id, ignore.case = TRUE), 
                 zebrafish_gene_id := toupper(zebrafish_gene_id)]

# Create a copy of the original data (CHECK: should now have 25 158 zebrafish genes)
hgnc_orthologsDF_original <- copy(hgnc_orthologsDF)

# Keep track of collapsed/merged predictions
hgnc_orthologsDF <- hgnc_orthologsDF %>% mutate(row_id = row_number())

colnames <- colnames(hgnc_orthologsDF)

# Merge non-unique ortholog predictions occurring from different databases
# e.g., appl2-APPL2 (prediction made 3 separate times)
# Group identical (non-unique) ortholog predictions based on gene ID and chromosome of gene
# After grouping identical predictions, remove duplicate assertions made from the same database
# Recount true number of supporting databases
hgnc_orthologsDF <- hgnc_orthologsDF[, .(
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
setcolorder(hgnc_orthologsDF, colnames)

# Identify rows with more than one ID (i.e., collapsed rows)
collapsed_rows <- hgnc_orthologsDF[sapply(strsplit(row_id, ","), length) > 1]

# Drop rows where zebrafish_gene_id is "-" (i.e., no entry)
hgnc_orthologsDF <- hgnc_orthologsDF[zebrafish_gene_id != "-"]        # 132 149 entries (529 less entries, no genes lost)

# HGNC Data -------------------------------------------------------------------

# Count occurrences and list all matching zebrafish genes for each human gene
HGNChuman_zebrafish <- hgnc_orthologsDF %>%
  # 18 595 unique human gene IDs (no change after standardizing nomenclature or merge)
  group_by(human_gene_id) %>%
  summarise(
    # Concatenate all unique human symbols associated with each human gene ID, excluding "-"
    human_symbol = paste(setdiff(unique(human_symbol), "-"), collapse = ", "),
    human_entrez_gene = paste(setdiff(unique(human_entrez_gene), "-"), collapse = ", "),
    human_ensembl_gene = paste(setdiff(unique(human_ensembl_gene), "-"), collapse = ", "),
    zebrafish_count = n_distinct(zebrafish_gene_id),
    # Concatenate all unique zebrafish gene IDs
    zebrafish_gene_id = paste(setdiff(unique(zebrafish_gene_id), "-"), collapse = ", "),
    zebrafish_symbol = paste(setdiff(unique(zebrafish_symbol), "-"), collapse = ", "),
    zebrafish_entrez_gene = paste(setdiff(unique(zebrafish_entrez_gene), "-"), collapse = ", "),
    zfin_id = paste(setdiff(unique(zfin_id), "-"), collapse = ", "),
    zebrafish_ensembl_gene = paste(setdiff(unique(zebrafish_ensembl_gene), "-"), collapse = ", "),
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

# Count occurrences and list all matching human genes for each zebrafish gene
HGNCzebrafish_human <- hgnc_orthologsDF %>%
  # 25 158 unique zebrafish gene IDs (down from 25 209 after standardizing, no loss from merge)
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

# # Store as .tsv (June 6, 2024)
# write.table(HGNChuman_zebrafish, file = "HGNChuman_zebrafishOrthologs.tsv",
#             sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
# 
# write.table(HGNCzebrafish_human, file = "HGNCzebrafish_humanOrthologs.tsv",
#             sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

# Single Gene Ortholog analysis
# Filter to keep only unique, one-to-one mappings in both datasets
one_to_one_human <- HGNChuman_zebrafish %>%
  filter(zebrafish_count == 1) %>%
  dplyr::select(human_gene_id, zebrafish_gene_id)

one_to_one_zebrafish <- HGNCzebrafish_human %>%
  filter(human_count == 1) %>%
  dplyr::select(zebrafish_gene_id, human_gene_id)

# Merge the datasets on the gene IDs, ensuring only those that have a one-to-one relationship in both are kept
SGO_analysis <- inner_join(one_to_one_human, one_to_one_zebrafish, 
                           by = c("human_gene_id" = "human_gene_id", 
                                  "zebrafish_gene_id" = "zebrafish_gene_id"))

# Mark these as True SGOs
SGO_status <- SGO_analysis %>%
  mutate(SGO_status = "True SGO") %>%
  dplyr::select(human_gene_id, zebrafish_gene_id, SGO_status)

# Join SGO status to ortholog data
HGNChuman_zebrafish <- HGNChuman_zebrafish %>% 
  left_join(SGO_status[, c("human_gene_id", "SGO_status")], by = "human_gene_id") %>%
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

print(average_one_to_many) # 11.0 zebrafish/9.7 humans

HGNCzebrafish_human %>% filter(average_support > 6) %>% summarise(count = n())

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
manualHGNC_orthologs[, human_gene_id := as.character(human_gene_id)]
subset_HGNChuman_zebrafish[, human_gene_id := as.character(human_gene_id)]

# Create a subset of manualHGNC_orthologs where human_gene_id appears in subset
# Where multiple orthology assertions are made, only assertions with manual support are kept
subset_manualHGNC <- manualHGNC_orthologs[human_gene_id %in% subset_HGNChuman_zebrafish$human_gene_id]
subset_manualHGNC <- subset_manualHGNC[zebrafish_gene_id != "-"]        # 8 315 genes

# Subset of orthologs without manual support
# i.e., assertions from the subset of "Many to Many" but no manual support
subsetNot_manualHGNC <- subset_HGNChuman_zebrafish[!(human_gene_id %in% subset_manualHGNC$human_gene_id)]
subsetNot_manualHGNC <- hgnc_orthologsDF[human_gene_id %in% subsetNot_manualHGNC$human_gene_id]
subsetNot_manualHGNC <- subsetNot_manualHGNC[zebrafish_gene_id != "-"]  # 3 062 genes (8315 + 3062 = 11 377, correct)

# Append data
subset_manualHGNC <- rbind(subset_manualHGNC, subsetNot_manualHGNC)

# Count occurrences and list all matching zebrafish genes for each human gene
subset_HGNChuman_zebrafish <- subset_manualHGNC %>%
  # 11 377 unique human gene IDs
  group_by(human_gene_id) %>%
  summarise(
    # Concatenate all unique human symbols associated with each human gene ID, excluding "-"
    human_symbol = paste(setdiff(unique(human_symbol), "-"), collapse = ", "),
    human_entrez_gene = paste(setdiff(unique(human_entrez_gene), "-"), collapse = ", "),
    human_ensembl_gene = paste(setdiff(unique(human_ensembl_gene), "-"), collapse = ", "),
    zebrafish_count = n_distinct(zebrafish_gene_id),
    # Concatenate all unique zebrafish gene IDs
    zebrafish_gene_id = paste(setdiff(unique(zebrafish_gene_id), "-"), collapse = ", "),
    zebrafish_symbol = paste(setdiff(unique(zebrafish_symbol), "-"), collapse = ", "),
    zebrafish_entrez_gene = paste(setdiff(unique(zebrafish_entrez_gene), "-"), collapse = ", "),
    zfin_id = paste(setdiff(unique(zfin_id), "-"), collapse = ", "),
    zebrafish_ensembl_gene = paste(setdiff(unique(zebrafish_ensembl_gene), "-"), collapse = ", "),
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

# Ensure zebrafish_gene_id columns are of the same type (character)
manualHGNC_orthologs[, zebrafish_gene_id := as.character(zebrafish_gene_id)]
subset_HGNCzebrafish_human[, zebrafish_gene_id := as.character(zebrafish_gene_id)]

# Create a subset of manualHGNC_orthologs where zebrafish_gene_id appears in subset
subset_manualHGNC <- manualHGNC_orthologs[zebrafish_gene_id %in% subset_HGNCzebrafish_human$zebrafish_gene_id]
subset_manualHGNC <- subset_manualHGNC[zebrafish_gene_id != "-"]        # 8 342 genes

# Subset of orthologs without manual support
# i.e., assertions from the subset of "Many to Many" but no manual support
subsetNot_manualHGNC <- subset_HGNCzebrafish_human[!(zebrafish_gene_id %in% subset_manualHGNC$zebrafish_gene_id)]
subsetNot_manualHGNC <- hgnc_orthologsDF[zebrafish_gene_id %in% subsetNot_manualHGNC$zebrafish_gene_id]
subsetNot_manualHGNC <- subsetNot_manualHGNC[zebrafish_gene_id != "-"]  # 3 868 genes (8342 + 3906 = 12 248, correct)

# Append data
subset_manualHGNC <- rbind(subset_manualHGNC, subsetNot_manualHGNC)

# Count occurrences and list all matching human genes for each zebrafish gene
subset_HGNCzebrafish_human <- subset_manualHGNC %>%
  # 12 248 unique zebrafish gene IDs
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

# Check if there are any new SGOs
one_to_one_human <- subset_HGNChuman_zebrafish %>%
  filter(zebrafish_count == 1) %>%
  dplyr::select(human_gene_id, zebrafish_gene_id)

one_to_one_zebrafish <- subset_HGNCzebrafish_human %>%
  filter(human_count == 1) %>%
  dplyr::select(zebrafish_gene_id, human_gene_id)

# Merge the datasets on the gene IDs, ensuring only those that have a one-to-one relationship in both are kept
SGO_analysis <- inner_join(one_to_one_human, one_to_one_zebrafish, 
                           by = c("human_gene_id" = "human_gene_id", 
                                  "zebrafish_gene_id" = "zebrafish_gene_id"))

# Mark these as True SGOs
SGO_status <- SGO_analysis %>%
  mutate(SGO_status = "True SGO") %>%
  dplyr::select(human_gene_id, zebrafish_gene_id, SGO_status)

# Join SGO status to ortholog data
subset_HGNChuman_zebrafish <- subset_HGNChuman_zebrafish %>% 
  left_join(SGO_status[, c("human_gene_id", "SGO_status")], by = "human_gene_id") %>%
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

print(average_one_to_many) # 15.8 zebrafish/23.7 humans

# # Store as .tsv (June 10, 2024)
# write.table(filteredHGNChuman_zebrafish, file = "filteredHGNChuman_zebrafishOrthologs.tsv",
#             sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
# 
# write.table(filteredHGNCzebrafish_human, file = "filteredHGNCzebrafish_humanOrthologs.tsv",
#             sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

# Visualize HGNC Data ---------------------------------------------------------

# Bar plot
colours <- c("One to One" = "#EE5496", "One to Many" = "#FFC000", "Many to One" = "#78DAD5")
lighter_colours <- alpha(colours, 0.7)

# Convert back to datafram if throwing errors
# HGNChuman_zebrafish <- as.data.frame(HGNChuman_zebrafish)

## Human to Zebrafish Orthologs including SGO status
ggplot(filteredHGNChuman_zebrafish, aes(x = relation_type, fill = relation_type)) +
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
    limits = c(0, 25000),       
    breaks = seq(0, 25000, by = 2500)
  )

ggsave("Filtered HCOP HGNC Symbol Human-Zebrafish.png", plot = plot, device = "png", 
       width = 9.50, height = 8.19, dpi = 600)

## Zebrafish to Human Orthologs including SGO status
ggplot(filteredHGNCzebrafish_human, aes(x = relation_type, fill = relation_type)) +
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
    limits = c(0, 25000),       
    breaks = seq(0, 25000, by = 2500)
  )

ggsave("Filtered HCOP HGNC Zebrafish-Human.png", plot = plot, device = "png", 
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

# Create a data frame for WP2446_humanSymbol
RBhuman_genes <- data.frame(
  human_symbol = WP2446_humanSymbol,
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

# Pathway Visualization -------------------------------------------------------

# Define the color palette
palette_colors <- c("#886889", "#B46D78", "#F89771", "#FCBF68", "#FAE079", "#ECDB66", "#E0D567", "#9EAB7C", "#81908B", "#7999B2")

# Generate 89 colors from the given palette
color_palette <- colorRampPalette(palette_colors)(89)
names(color_palette) <- as.character(1:89)
color_palette <- c(color_palette, "zero" = "#9400D3") 

# Ensure the color palette matches the levels of color_group
color_levels <- c(as.character(1:89), "zero")

# Gray color for the size category legend
gray_color <- "#252525"

# Create the bubble plot (shape = 21 is a filled circle with border)
create_bubblePlot <- function(data){
  ggplot(data, aes(x = human_symbol, y = composite_score)) +
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
  ggplot(data, aes(x = human_symbol, y = composite_score)) +
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
plot2a<-create_bubblePlot(filtered_data2)
plot3<-bubble_noLegend(filtered_data3)
plot4<-bubble_noLegend(filtered_data4)

# Create a common y-axis title plot
y_title <- ggplot() + 
  annotate("text", x = 1, y = 0.5, label = "Composite Score", angle = 90, size = 4) +
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

combined_plot <- (y_title | wrap_plots(plot2a, plot4, plot3, ncol = 1)) +
  plot_layout(widths = c(0.05, 1))

final_plot <- combined_plot / x_title +
  plot_layout(heights = c(1, 0.05))

# Print the final plot
print(final_plot)

ggsave("Overall Confidence in Orthology per Gene x3.png", plot = final_plot, device = "png", 
       width = 13.34, height = 8, dpi = 600)

# 6.96

# Internal Validation Work -----------------------------------------------------
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