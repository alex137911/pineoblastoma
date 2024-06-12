suppressMessages(library(rtracklayer))
suppressMessages(library(readr))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(ggplot2))


inpath <- "C:/Users/acale/OneDrive/Documents/Waterloo BME/Co-op/4C/pineoblastoma/Data/Input/"
inDir  <- sprintf("%s/Input", dirname(inpath))
setwd(inDir)

# Load and Clean Data ---------------------------------------------------------
HGNChuman_ZF <- read_delim("HGNC_human_ZF_genenamesonly.tsv.gz", 
                           delim = "\t", escape_double = FALSE, trim_ws = TRUE)

HGNChuman_ZFdf <- as.data.frame(HGNChuman_ZF)


# Separate zebrafish and human IDs
# Warning from rows without corresponding HGNC Human ID
HGNChuman_ZFdf <- HGNChuman_ZFdf %>%
  separate(col = ensHS_Gene.name, into = c("zebrafish_symbol", "human_symbol"), sep = "\t")

# Clean data (remove characters external to ID)
HGNChuman_ZFdf <- HGNChuman_ZFdf %>%
  mutate(zebrafish_symbol = gsub('\"', '', zebrafish_symbol), 
         human_symbol = gsub('\"', '', human_symbol))

# Clean data (remove characters external to ID)
HGNChuman_ZFdf <- HGNChuman_ZFdf %>%
  mutate(zebrafish_symbol = gsub('\"', '', zebrafish_symbol), 
         human_symbol = gsub('\"', '', human_symbol))

# Drop NAs (rows without a corresponding Human ortholog)
HGNChuman_ZFdf <- HGNChuman_ZFdf %>% filter(!is.na(human_symbol))

# Drop rows with blank entries (errors in data entry?)
HGNChuman_ZFdf <- HGNChuman_ZFdf %>% filter(human_symbol != "")

# Remove Lawson Gene Symbol column (probably corrupted when emailed over?)
HGNChuman_ZFdf <- HGNChuman_ZFdf %>% 
  dplyr::select(zebrafish_symbol, human_symbol)

# Add a unique ID to each row for verification
# Curate list of duplicates for Oliva to validate
HGNChuman_ZFdf <- HGNChuman_ZFdf %>% 
  mutate(id = row_number()) %>% dplyr::select(id, everything())

outpath <- "C:/Users/acale/OneDrive/Documents/Waterloo BME/Co-op/4C/pineoblastoma/Data/Output/"
outDir  <- sprintf("%s/Output", dirname(outpath))
setwd(outDir)

# # Write data (May 29, 2024)
# write.table(HGNChuman_ZFdf, file = "farrellHuman_zebrafishOrthologs.tsv",
#             sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

# Count orthologs  ------------------------------------------------------------
# Count occurrences and list all matching zebrafish genes for each human gene
farrellHuman_zebrafish <- HGNChuman_ZFdf %>%
  # 2 654 unique gene (HGNC) IDs
  group_by(human_symbol) %>%
  summarise(
    zebrafish_count = n_distinct(zebrafish_symbol),
    # Concatenate all unique zebrafish gene IDs
    zebrafish_symbol = paste(unique(zebrafish_symbol), collapse = ", "),
    .groups = 'drop'
  ) %>%
  mutate(relation_type = case_when(
    zebrafish_count == 1 ~ "One to One",
    zebrafish_count > 1 ~ "One to Many",
    TRUE ~ "No Match"
  ))

# Count occurrences and list all matching human genes for each zebrafish gene
farrellZebrafish_human <- HGNChuman_ZFdf %>%
  # 2 110 unique gene IDs
  group_by(zebrafish_symbol) %>%
  summarise(
    human_count = n_distinct(human_symbol),
    # Concatenate all unique zebrafish gene IDs
    human_symbol = paste(unique(human_symbol), collapse = ", "),
    .groups = 'drop'
  ) %>%
  mutate(relation_type = case_when(
    human_count == 1 ~ "One to One",
    human_count > 1 ~ "Many to One",
    TRUE ~ "No Match"
  ))

# Single Gene Ortholog analysis
# Filter to keep only unique, one-to-one mappings in both datasets
one_to_one_human <- farrellHuman_zebrafish %>%
  filter(zebrafish_count == 1) %>%
  dplyr::select(human_symbol, zebrafish_symbol)

one_to_one_zebrafish <- farrellZebrafish_human %>%
  filter(human_count == 1) %>%
  dplyr::select(zebrafish_symbol, human_symbol)

# Merge the datasets on the gene IDs, ensuring only those that have a one-to-one relationship in both are kept
SGO_analysis <- inner_join(one_to_one_human, one_to_one_zebrafish, 
                           by = c("human_symbol" = "human_symbol", 
                                  "zebrafish_symbol" = "zebrafish_symbol"))

# Mark these as True SGOs
SGO_status <- SGO_analysis %>%
  mutate(SGO_status = "True SGO") %>%
  dplyr::select(human_symbol, zebrafish_symbol, SGO_status)

# Join SGO status to ortholog data
farrellHuman_zebrafish <- farrellHuman_zebrafish %>% 
  left_join(SGO_status[, c("human_symbol", "SGO_status")], by = "human_symbol") %>%
  mutate(SGO_status = ifelse(is.na(SGO_status), "Not SGO", SGO_status))

farrellZebrafish_human <- farrellZebrafish_human %>% 
  left_join(SGO_status[, c("zebrafish_symbol", "SGO_status")], by = "zebrafish_symbol") %>%
  mutate(SGO_status = ifelse(is.na(SGO_status), "Not SGO", SGO_status))

# Average orthologs
average_one_to_many <- farrellHuman_zebrafish %>%
  filter(relation_type == "One to Many") %>%        # Filter to keep only 'One to Many' rows
  summarise(average_count = mean(zebrafish_count))  

average_one_to_many <- farrellZebrafish_human %>%
  filter(relation_type == "Many to One") %>%        # Filter to keep only 'One to Many' rows
  summarise(average_count = mean(human_count))  

print(average_one_to_many) # 5.7 zebrafish/7.7 humans

# Visualize Data --------------------------------------------------------------
colours <- c("One to One" = "#EE5496", "One to Many" = "#FFC000", "Many to One" = "#78DAD5")
lighter_colours <- alpha(colours, 0.7)

# Human to Zebrafish Orthologs including SGO status
plot<-ggplot(farrellHuman_zebrafish, aes(x = relation_type, fill = relation_type)) +
  geom_bar(stat = "count", position = position_dodge(), color = colours[2:1], size = 1, alpha = 0.5) +
  geom_text(stat = 'count', aes(label = ..count..), 
            vjust = -0.5, position = position_dodge(width = 0.9), size = 7) +
  geom_bar(data = subset(farrellHuman_zebrafish, SGO_status == "True SGO"), 
           aes(x = relation_type, fill = relation_type),
           stat = "count", position = position_dodge(), color = colours[1], size = 1) +
  geom_text(data = subset(farrellHuman_zebrafish, SGO_status == "True SGO"), aes(label = ..count..), 
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
    limits = c(0, 2000),       
    breaks = seq(0, 2000, by = 500)
  )

ggsave("Farrell Lab Human-Zebrafish.png", plot = plot, device = "png", 
       width = 9.50, height = 8.19, dpi = 600)

# Zebrafish to Human Orthologs including SGO status
plot<-ggplot(farrellZebrafish_human, aes(x = relation_type, fill = relation_type)) +
  geom_bar(stat = "count", position = position_dodge(), color = colours[c(3, 1)], size = 1, alpha = 0.5) +
  geom_text(stat = 'count', aes(label = ..count..), 
            vjust = -0.5, position = position_dodge(width = 0.9), size = 7) +
  geom_bar(data = subset(farrellZebrafish_human, SGO_status == "True SGO"), 
           aes(x = relation_type, fill = relation_type),
           stat = "count", position = position_dodge(), color = colours[1], size = 1) +
  geom_text(data = subset(farrellZebrafish_human, SGO_status == "True SGO"), aes(label = ..count..), 
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
    limits = c(0, 2000),       
    breaks = seq(0, 2000, by = 500)
  )

ggsave("Farrell Lab Zebrafish-Human.png", plot = plot, device = "png", 
       width = 9.50, height = 8.19, dpi = 600)