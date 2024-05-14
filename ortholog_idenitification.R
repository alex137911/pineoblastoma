# Human-Zebrafish Orthologs

# Remove objects in workspace
rm(list = ls())

# Required packages
suppressMessages(library(readr))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(scales))
suppressMessages(library(data.table))
# suppressMessages(library(Matrix))
suppressMessages(library(tidyr))
suppressMessages(library(stringr))
# suppressMessages(library(pheatmap))
# suppressMessages(library(igraph))
# suppressMessages(library(ggraph))
suppressMessages(library(plotly))

# -------------------------------------------------------------------
# LOAD DATA
inpath <- "C:/Users/acale/OneDrive/Documents/Waterloo BME/Co-op/4C/pineoblastoma/Data/Input/"
inDir  <- sprintf("%s/Input", dirname(inpath))
setwd(inDir)

# Gene Conversion Table
# From: 10.1016/j.devcel.2023.11.001
# Jeffrey Farrell paper
HGNChuman_ZF <- read_delim("HGNC_human_ZF_genenamesonly.tsv.gz", 
                           delim = "\t", escape_double = FALSE, trim_ws = TRUE)

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
ensemblHZ_orthologs <- read_delim("human-zebrafish_biomart_export.txt", 
                                     delim = "\t", escape_double = FALSE, trim_ws = TRUE)

ensemblHZ_orthologsDF <- as.data.frame(ensemblHZ_orthologs)
message(sprintf("Loaded %i records", nrow(ensemblHZ_orthologsDF)))       # 330 460 records

ensemblZH_orthologs <- read_delim("zebrafish-human_biomart_export.txt",
                                  delim = "\t", escape_double = FALSE, trim_ws = TRUE)

ensemblZH_orthologsDF <- as.data.frame(ensemblZH_orthologs)
message(sprintf("Loaded %i records", nrow(ensemblZH_orthologsDF)))       # 72 462 records

# Import ZFIN Data
# Downloaded May 8, 2024 from: https://zfin.org/downloads
header_column <- c("ZFIN_ID", "ZFIN_Symbol", "ZFIN_Name", "Human_Symbol", 
                   "Human_Name", "OMIM_ID", "Gene_ID", "HGNC_ID", "Evidence", "Pub_ID")

ZFIN_orthologs <- read_delim("human_zebrafish_zfin.txt",
                             delim = "\t", col_names = header_column, escape_double = FALSE, trim_ws = TRUE)

ZFIN_orthologsDF <- as.data.frame(ZFIN_orthologs)
message(sprintf("Loaded %i records", nrow(ZFIN_orthologsDF)))            # 43 612 records 

# Import pineoblastoma data (from Laurie)
pineo_primaryTumours <- read_delim("pineo.primary_tumours.20240223.fpkm.txt", 
                                   delim = "\t", escape_double = FALSE, trim_ws = TRUE)

# Retinoblastoma gene in cancer (WP2446) pathway
# Downloaded May 10, 2024 from https://www.wikipathways.org/pathways/WP2446.html
# Pathway upregulated in in PB tumours (RB subgroup) - from Oliva
WP2446_pathwayGenes <- read_delim("WP2446_pathwayGenes.tsv", delim = "\t", escape_double = FALSE,
                                  trim_ws = TRUE)

# Drop ensembl tag
WP2446_pathwayGenes <- WP2446_pathwayGenes %>% mutate(Ensembl = sub("ensembl:", "", Ensembl))

# Quality Check --------------------------------------------

# HGNC 
no_entries <- sum(hgnc_orthologsDF$zebrafish_entrez_gene == "-")  # 11 645 blank
no_entries <- sum(hgnc_orthologsDF$zebrafish_ensembl_gene == "-") # 2 081 blank
no_entries <- sum(hgnc_orthologsDF$zebrafish_name == "-")         # 5 253 blank
no_entries <- sum(hgnc_orthologsDF$zfin_id == "-")                # 16 005 blank

# Drop rows with no zebrafish_ensembl_gene 
hgnc_orthologsDF <- hgnc_orthologsDF %>% filter(zebrafish_ensembl_gene != "-")

# HZ BioMart
no_entries <- sum(is.na(ensemblHZ_orthologsDF$`Zebrafish gene stable ID`)) # 138 697 blank
no_entries <- sum(is.na(ensemblHZ_orthologsDF$`Zebrafish gene name`))      # 167 370 blank

# Drop rows with no zebrafish_ensembl_gene (Zebrafish gene stable ID)
ensemblHZ_orthologsDF <- ensemblHZ_orthologsDF %>% 
  filter(!is.na(ensemblHZ_orthologsDF$`Zebrafish gene stable ID`))

# ZH BioMart
no_entries <- sum(is.na(ensemblZH_orthologsDF$`Human gene stable ID`))     # 28 670 blank
no_entries <- sum(is.na(ensemblZH_orthologsDF$`Human gene name`))          # 28 847 blank

# Drop rows with no human_ensembl_gene (Human gene stable ID)
ensemblZH_orthologsDF <- ensemblZH_orthologsDF %>% 
  filter(!is.na(ensemblZH_orthologsDF$`Human gene stable ID`))

# HGNC Data -----------------------------------------------

# Count occurrences and list all matching zebrafish genes for each human gene
HGNChuman_zebrafish <- hgnc_orthologsDF %>%
  group_by(human_ensembl_gene) %>%
  summarise(
    human_symbol = first(human_symbol),         # Assumes there is a unique human_symbol
    zebrafish_count = n(),
    zebrafish_ensembl_gene = paste(zebrafish_ensembl_gene, collapse = ", "),
    .groups = 'drop'
  ) %>%
  mutate(relation_type = case_when(
    zebrafish_count == 1 ~ "One to One",
    zebrafish_count > 1 ~ "One to Many",
    TRUE ~ "No Match"
  ))

# Count occurrences and list all matching human genes for each zebrafish gene
HGNCzebrafish_human <- hgnc_orthologsDF %>%
  group_by(zebrafish_ensembl_gene) %>%
  summarise(
    zebrafish_symbol = first(zebrafish_symbol),  # Assumes there is a unique zebrafish_symbol (inaccurate)
    human_count = n(),
    human_ensembl_gene = paste(human_ensembl_gene, collapse = ", "),
    .groups = 'drop'
  ) %>%
  mutate(relation_type = case_when(
    human_count == 1 ~ "One to One",
    human_count > 1 ~ "Many to One",
    TRUE ~ "No Match"
  ))

# Single Gene Ortholog analysis
# Filter to keep only unique, one-to-one mappings in both datasets
one_to_one_human <- HGNChuman_zebrafish %>%
  filter(zebrafish_count == 1) %>%
  select(human_ensembl_gene, zebrafish_ensembl_gene)

one_to_one_zebrafish <- HGNCzebrafish_human %>%
  filter(human_count == 1) %>%
  select(zebrafish_ensembl_gene, human_ensembl_gene)

# Merge the datasets on the gene IDs, ensuring only those that have a one-to-one relationship in both are kept
SGO_analysis <- inner_join(one_to_one_human, one_to_one_zebrafish, 
                           by = c("human_ensembl_gene" = "human_ensembl_gene", 
                                  "zebrafish_ensembl_gene" = "zebrafish_ensembl_gene"))

# Mark these as True SGOs
SGO_status <- SGO_analysis %>%
  mutate(SGO_status = "True SGO") %>%
  select(human_ensembl_gene, zebrafish_ensembl_gene, SGO_status)

# Join SGO status to ortholog data
HGNChuman_zebrafish <- HGNChuman_zebrafish %>% 
  left_join(SGO_status[, c("human_ensembl_gene", "SGO_status")], by = "human_ensembl_gene") %>%
  mutate(SGO_status = ifelse(is.na(SGO_status), "Not SGO", SGO_status))

HGNCzebrafish_human <- HGNCzebrafish_human %>% 
  left_join(SGO_status[, c("zebrafish_ensembl_gene", "SGO_status")], by = "zebrafish_ensembl_gene") %>%
  mutate(SGO_status = ifelse(is.na(SGO_status), "Not SGO", SGO_status))

# Average orthologs
average_one_to_many <- HGNChuman_zebrafish %>%
  filter(relation_type == "One to Many") %>%  # Filter to keep only 'One to Many' rows
  summarise(average_count = mean(zebrafish_count))  # Calculate the average of zebrafish_count

average_one_to_many <- HGNCzebrafish_human %>%
  filter(relation_type == "Many to One") %>%  # Filter to keep only 'One to Many' rows
  summarise(average_count = mean(human_count))  # Calculate the average of zebrafish_count

print(average_one_to_many) # 11 zebrafish/9.7 humans

# ----------------------------------------------------------
# HGNC Pathway Analysis
# Keep only human-zebrafish orthologs with pineoblastoma genes
no_entries <- sum(HGNChuman_zebrafish$human_symbol == "-")                # 199 genes with no HGNC Id

HGNChuman_zebrafishPB <- HGNChuman_zebrafish %>% 
  filter(human_symbol %in% pineo_primaryTumours$gene)

# Keep only human-zebrafish orthologs in WP2446 pathway (RB)
HGNChuman_zebrafishWP2446 <- HGNChuman_zebrafish %>%
  filter(human_ensembl_gene %in% WP2446_pathwayGenes$Ensembl)

HGNChuman_zebrafishWP2446 <- subset(HGNChuman_zebrafishWP2446, human_symbol != "ZNF655")
HGNChuman_zebrafishWP2446 <- subset(HGNChuman_zebrafishWP2446, relation_type != "One to One")

# Visualize HGNC Data -----------------------------------------------

# Bar plot
colours <- c("One to One" = "#EE5496", "One to Many" = "#FFC000", "Many to One" = "#78DAD5")
lighter_colours <- alpha(colours, 0.7)

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
    limits = c(0, 15000),       
    breaks = seq(0, 15000, by = 2500)
  )

ggsave("HGNC Human-Zebrafish.png", plot = plot, device = "png", 
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
    limits = c(0, 15000),       
    breaks = seq(0, 15000, by = 2500)
  )

ggsave("HGNC Zebrafish-Human.png", plot = plot, device = "png", 
       width = 9.50, height = 8.19, dpi = 600)

# Pathway plots
ggplot(HGNChuman_zebrafishPB, aes(x=relation_type, fill=relation_type)) +
  geom_bar() +
  geom_text(stat='count', aes(label=..count..), vjust=-0.5, position=position_stack(vjust=1.0)) +
  labs(title="Distribution of PB Gene Relationship Types", x="Relationship Type", y="Count") +
  theme_minimal()

ggplot(HGNChuman_zebrafishWP2446, aes(x=relation_type, fill=relation_type)) +
  geom_bar() +
  geom_text(stat='count', aes(label=..count..), vjust=-0.5, position=position_stack(vjust=1.0)) +
  labs(title="Distribution of WP2446 Gene Relationship Types", x="Relationship Type", y="Count") +
  theme_minimal()

# Dot Plot (WP2446 Pathway)
ggplot(HGNChuman_zebrafishWP2446, aes(x=human_symbol, y=zebrafish_count, color=relation_type)) +
  geom_point() +
  labs(title="Human to Zebrafish Gene Relationships", x="Human Gene", y="Number of Zebrafish Genes") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle=90, hjust=1))

# Visualization
# Expand the zebrafish genes into separate rows
expanded_WP2446orthologs <- HGNChuman_zebrafishWP2446 %>%
  mutate(zebrafish_ensembl_genes = str_split(zebrafish_ensembl_genes, ", ")) %>%
  unnest(zebrafish_ensembl_genes)

# Creating a summary of counts
gene_counts <- expanded_WP2446orthologs %>%
  group_by(human_symbol) %>%
  summarise(count = n(), .groups = 'drop') %>%  # This line aggregates and counts the occurrences
  arrange(count)

# Reorder the factor levels based on the new ordering
gene_counts$human_symbol <- factor(gene_counts$human_symbol, levels = gene_counts$human_symbol)
  
# Bar Plot
ggplot(gene_counts, aes(x = human_symbol, y = count, fill = human_symbol)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), 
        legend.position="none", 
        panel.background = element_rect(fill="white"),
        panel.grid = element_line(colour="lightgrey"))+
  labs(x = "Human Gene", y = "Number of Zebrafish Orthologs", title = "Distribution of Orthologs per Human Gene")

# BioMart Data -----------------------------------------------

# Count occurrences and list all matching zebrafish genes for each human gene
ensemblHuman_zebrafish <- ensemblHZ_orthologsDF %>%
  group_by(`Gene stable ID`) %>%
  summarise(
    # human_symbol = first(human_symbol),         # Assumes there is a unique human_symbol
    zebrafish_count = n(),
    zebrafish_ensembl_gene = paste(`Zebrafish gene stable ID`, collapse = ", "),
    transcript_id = paste(`Transcript stable ID`, collapse = ", "),
    .groups = 'drop'
  ) %>%
  mutate(relation_type = case_when(
    zebrafish_count == 1 ~ "One to One",
    zebrafish_count > 1 ~ "One to Many",
    TRUE ~ "No Match"
  ))

# Count occurrences and list all matching human genes for each zebrafish gene
ensemblZebrafish_human <- ensemblZH_orthologsDF %>%
  group_by(`Gene stable ID`) %>%
  summarise(
    human_count = n(),
    human_ensembl_gene = paste(`Human gene stable ID`, collapse = ", "),
    transcript_id = paste(`Transcript stable ID`, collapse = ", "),
    .groups = 'drop'
  ) %>%
  mutate(relation_type = case_when(
    human_count == 1 ~ "One to One",
    human_count > 1 ~ "Many to One",
    TRUE ~ "No Match"
  ))

# Single Gene Ortholog analysis
# Filter to keep only unique, one-to-one mappings in both datasets
one_to_one_human <- ensemblHuman_zebrafish %>%
  filter(zebrafish_count == 1) %>%
  select(`Gene stable ID`, zebrafish_ensembl_gene)

one_to_one_zebrafish <- ensemblZebrafish_human %>%
  filter(human_count == 1) %>%
  select(`Gene stable ID`, human_ensembl_gene)

# Rename columns for the join
one_to_one_human <- rename(one_to_one_human, human_ensembl_gene = `Gene stable ID`)
one_to_one_zebrafish <- rename(one_to_one_zebrafish, zebrafish_ensembl_gene = `Gene stable ID`)
ensemblHuman_zebrafish <- rename(ensemblHuman_zebrafish, human_ensembl_gene = `Gene stable ID`)
ensemblZebrafish_human <- rename(ensemblZebrafish_human, zebrafish_ensembl_gene = `Gene stable ID`)

# Merge the datasets on the gene IDs, ensuring only those that have a one-to-one relationship in both are kept
SGO_analysis <- inner_join(one_to_one_human, one_to_one_zebrafish, 
                           by = c("human_ensembl_gene" = "human_ensembl_gene", 
                                  "zebrafish_ensembl_gene" = "zebrafish_ensembl_gene"))

# Mark these as True SGOs
SGO_status <- SGO_analysis %>%
  mutate(SGO_status = "True SGO") %>%
  select(human_ensembl_gene, zebrafish_ensembl_gene, SGO_status)

# Join SGO status to ortholog data
ensemblHuman_zebrafish <- ensemblHuman_zebrafish %>% 
  left_join(SGO_status[, c("human_ensembl_gene", "SGO_status")], by = "human_ensembl_gene") %>%
  mutate(SGO_status = ifelse(is.na(SGO_status), "Not SGO", SGO_status))

ensemblZebrafish_human <- ensemblZebrafish_human %>% 
  left_join(SGO_status[, c("zebrafish_ensembl_gene", "SGO_status")], by = "zebrafish_ensembl_gene") %>%
  mutate(SGO_status = ifelse(is.na(SGO_status), "Not SGO", SGO_status))

# Average orthologs
average_one_to_many <- ensemblHuman_zebrafish %>%
  filter(relation_type == "One to Many") %>%  # Filter to keep only 'One to Many' rows
  summarise(average_count = mean(zebrafish_count))  # Calculate the average of zebrafish_count

average_one_to_many <- ensemblZebrafish_human %>%
  filter(relation_type == "Many to One") %>%  # Filter to keep only 'One to Many' rows
  summarise(average_count = mean(human_count))  # Calculate the average of zebrafish_count

print(average_one_to_many) # 14.1 zebrafish/3.5 humans

# Visualize BioMart Data -----------------------------------------------

# Bar Plots
## Human to Zebrafish Orthologs including SGO status
plot<-ggplot(ensemblHuman_zebrafish, aes(x = relation_type, fill = relation_type)) +
  geom_bar(stat = "count", position = position_dodge(), color = colours[2:1], size = 1, alpha = 0.5) +
  geom_text(stat = 'count', aes(label = ..count..), 
            vjust = -0.5, position = position_dodge(width = 0.9), size = 7) +
  geom_bar(data = subset(ensemblHuman_zebrafish, SGO_status == "True SGO"), 
           aes(x = relation_type, fill = relation_type),
           stat = "count", position = position_dodge(), color = colours[1], size = 1) +
  geom_text(data = subset(ensemblHuman_zebrafish, SGO_status == "True SGO"), aes(label = ..count..), 
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
    limits = c(0, 15000),       
    breaks = seq(0, 15000, by = 2500)
  )

ggsave("Ensembl Human-Zebrafish.png", plot = plot, device = "png", 
       width = 9.50, height = 8.19, dpi = 600)

## Zebrafish to Human Orthologs including SGO status
plot<-ggplot(ensemblZebrafish_human, aes(x = relation_type, fill = relation_type)) +
  geom_bar(stat = "count", position = position_dodge(), color = colours[c(3, 1)], size = 1, alpha = 0.5) +
  geom_text(stat = 'count', aes(label = ..count..), 
            vjust = -0.5, position = position_dodge(width = 0.9), size = 7) +
  geom_bar(data = subset(ensemblZebrafish_human, SGO_status == "True SGO"), 
           aes(x = relation_type, fill = relation_type),
           stat = "count", position = position_dodge(), color = colours[1], size = 1) +
  geom_text(data = subset(ensemblZebrafish_human, SGO_status == "True SGO"), aes(label = ..count..), 
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
    limits = c(0, 15000),       
    breaks = seq(0, 15000, by = 2500)
  )

ggsave("Ensembl Zebrafish-Human.png", plot = plot, device = "png", 
       width = 9.50, height = 8.19, dpi = 600)
