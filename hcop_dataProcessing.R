# Remove objects in workspace
rm(list = ls())

# Required Packages -----------------------------------------------------------
suppressMessages(library(readr))
suppressMessages(library(dplyr))
suppressMessages(library(data.table))
suppressMessages(library(UpSetR))
suppressMessages(library(ggplot2))
suppressMessages(library(ggdist))
suppressMessages(library(gghalves))
suppressMessages(library(forcats))
suppressMessages(library(colorspace))
suppressMessages(library(patchwork))
suppressMessages(library(tidyr))
suppressMessages(library(irr))
suppressMessages(library(ggalluvial))

# Load Data -------------------------------------------------------------------
inpath <- "C:/Users/acale/OneDrive/Documents/Waterloo BME/Co-op/4C/pineoblastoma/Data/Input/"
inDir  <- sprintf("%s/Input", dirname(inpath))
setwd(inDir)

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

outpath <- "C:/Users/acale/OneDrive/Documents/Waterloo BME/Co-op/4C/pineoblastoma/Data/Output/"
outDir  <- sprintf("%s/Output", dirname(outpath))
setwd(outDir)

# hgnc_orthologsDF <- read_delim("hgnc_orthologsFiltered.tsv",
#                                delim = "\t", escape_double = FALSE, trim_ws = TRUE)

hgnc_orthologsDF <- read_delim("clustered_hgncOrthologs.tsv",
                               delim = "\t", escape_double = FALSE, trim_ws = TRUE)

# HCOP Clustering Analysis ----------------------------------------------------

sum(hgnc_orthologsDF$zebrafish_symbol %in% farrell_categorizedGenes$gene)

# Define the columns to select from hgnc_orthologsDF
columns_hgnc <- c("human_entrez_gene", "human_ensembl_gene", "human_symbol",
                  "zebrafish_entrez_gene", "zebrafish_ensembl_gene", "zfin_id", "zebrafish_symbol", 
                  "support", "cluster")

# Filter hgnc_orthologsDF based on the condition and select the specified columns
filtered_hgnc_orthologsDF <- hgnc_orthologsDF %>%
  filter(zebrafish_symbol %in% farrell_categorizedGenes$gene) %>%
  select(all_of(columns_hgnc)) %>%
  group_by(cluster) %>%
  mutate(cluster_count = n()) %>%
  ungroup()

# Subset by housekeeping genes
housekeeping_genes <- farrell_categorizedGenes %>%
  filter(class == "Housekeeping")

housekeeping_genes <- hgnc_orthologsDF %>%
  filter(zebrafish_symbol %in% housekeeping_genes$gene) %>%
  select(all_of(columns_hgnc)) %>%
  group_by(cluster) %>%
  mutate(cluster_count = n()) %>%
  ungroup()

housekeeping_distinct <- housekeeping_genes %>% 
  dplyr::distinct(cluster, .keep_all = TRUE)

# Subset by ubiquitous-varying genes (genes expressed in many cell types
# but with significant variability in expression levels)
ubiquitous_varying <- farrell_categorizedGenes %>%
  filter(class == "Ubiquitous varying")

ubiquitous_varying <- hgnc_orthologsDF %>%
  filter(zebrafish_symbol %in% ubiquitous_varying$gene) %>%
  select(all_of(columns_hgnc)) %>%
  group_by(cluster) %>%
  mutate(cluster_count = n()) %>%
  ungroup()

ubiquitous_distinct <- ubiquitous_varying %>% 
  dplyr::distinct(cluster, .keep_all = TRUE)

# Subset by tissue-specific/restricted genes
tissue_specRestricted <- farrell_categorizedGenes %>%
  filter(class.devcell %in% c("Tissue restricted", "Tissue-specific"))

tissue_specRestricted <- hgnc_orthologsDF %>%
  filter(zebrafish_symbol %in% tissue_specRestricted$gene) %>%
  select(all_of(columns_hgnc)) %>%
  group_by(cluster) %>%
  mutate(cluster_count = n()) %>%
  ungroup()

tissue_distinct <- tissue_specRestricted %>% 
  dplyr::distinct(cluster, .keep_all = TRUE)

# Filter by cell-specific/restricted genes
cell_specRestricted <- farrell_categorizedGenes %>%
  filter(class.devcell %in% c("Cell type-specific", "Cell type-restricted"))

cell_specRestricted <- hgnc_orthologsDF %>%
  filter(zebrafish_symbol %in% cell_specRestricted$gene) %>%
  select(all_of(columns_hgnc)) %>%
  group_by(cluster) %>%
  mutate(cluster_count = n()) %>%
  ungroup()

cell_distinct <- cell_specRestricted %>% 
  dplyr::distinct(cluster, .keep_all = TRUE)

# Visualize number of nodes per cluster ---------------------------------------
# Add a category column to each data frame
housekeeping_genes$category <- "Housekeeping"
ubiquitous_varying$category <- "Ubiquitous varying"
tissue_specRestricted$category <- "Tissue specific/restricted"
cell_specRestricted$category <- "Cell type specific/restricted"

# All genes
hgnc_distinct <- hgnc_orthologsDF %>%
  select(all_of(columns_hgnc)) %>%
  group_by(cluster) %>%
  mutate(cluster_count = n()) %>%
  ungroup() %>%
  dplyr::distinct(cluster, .keep_all = TRUE)

hgnc_distinct$category <- "All"

# Combine the data frames
combined_data <- bind_rows(housekeeping_genes, ubiquitous_varying, 
                           tissue_specRestricted, cell_specRestricted)

combined_distinct <- bind_rows(housekeeping_distinct, ubiquitous_distinct,
                               tissue_distinct, cell_distinct)

# Enhanced
# Define the function to sample points
sample_points <- function(data, proportion) {
  data %>%
    group_by(category) %>%
    sample_frac(proportion)
}

# Sample a subset of points to reduce overcrowding
sampled_data <- sample_points(combined_data, proportion = 0.1)
sample_distinct <- sample_points(combined_distinct, proportion = 0.1)
sample_all <- sample_points(hgnc_distinct, proportion = 0.005)

# Calculate position and label for stat_summary
add_sample <- function(x){
  return(c(y = max(x) + .025, 
           label = length(x)))
}

colours <- c("#EE5496", "#FFC000", "#78DAD5", "#A6D96D", "#C197E1")

# Function to create Raincloud plot for each category
create_raincloudPlot <- function(graph_data, sampled_data, category_name, graph_colour) {
  data_filtered <- graph_data %>% filter(category == category_name)
  y_range <- max(data_filtered$cluster_count) - min(data_filtered$cluster_count)
  y_position_median <- min(data_filtered$cluster_count) + 0.4 * y_range
  y_position_sample <- min(data_filtered$cluster_count) + 0.9 * y_range
  median_value <- median(data_filtered$cluster_count, na.rm = TRUE)
  
  ggplot(graph_data %>% filter(category == category_name), 
         aes(x = category, y = cluster_count, fill = graph_colour, color = graph_colour)) + 
    ggdist::stat_halfeye(aes(fill = after_scale(lighten(graph_colour, .5))),
                         adjust = 5, width = 1.2, .width = 0, 
                         justification = -.4, point_color = NA, alpha = 0.7) +
    geom_boxplot(aes(fill = after_scale(desaturate(lighten(graph_colour, .8), .4))),
                 width = .42, outlier.shape = NA, 
                 alpha = 0.5, color = graph_colour) +
    geom_jitter(data = sampled_data %>% filter(category == category_name), 
                color = graph_colour, shape = 16,
                size = 2, alpha = .25, 
                position = position_jitter(seed = 1, width = .12)) +
    # geom_jitter(data = sampled_data %>% filter(category == category_name), 
    #             fill = graph_colour, color = "transparent", shape = 21,
    #             stroke = 0, size = 2, alpha = .6, 
    #             position = position_jitter(seed = 1, width = .12)) +
    # stat_summary(geom = "text", fun = "median",
    #   aes(label = round(..y.., 2)), 
    #   color = after_scale(darken(graph_colour, .1, space = "HLS")),
    #   fontface = "bold", size = 4.5, hjust = -1) +
    # annotate("text", x = 1, y = y_position_median, label = round(median_value), 
    #          size = 4.5, fontface = "bold", 
    #          color = colorspace::darken(graph_colour, 0.3)) +
    stat_summary(geom = "text", 
                 fun.data = function(y) 
                   data.frame(y = y_position_sample, 
                              label = paste("n =", length(y))),
                 aes(label = ..label..), size = 4, color = "grey40") +
    #coord_cartesian(ylim = c(0, 100)) +
    labs(x = "ZF Gene Category", y = "Nodes per Cluster") +
    theme_minimal() +
    #theme_blank() +
    theme(
      axis.text.x = element_text(size = 10, vjust = 0.2, color = "grey40"), 
      axis.text.y = element_text(size = 10, color = "grey40"),
      plot.margin = margin(5, 5, 10, 5),
      legend.position = "none",
      panel.background = element_blank(),
      plot.background = element_rect(fill = "transparent", color = NA),
      # Remove axis titles for combined plot
      axis.title = element_blank()
    )
}

# Create individual plots for each category
plot_housekeeping <- create_raincloudPlot(combined_data, sampled_data, "Housekeeping", colours[1])
plot_ubiquitous_varying <- create_raincloudPlot(combined_data, sampled_data, "Ubiquitous varying", colours[2])
plot_tissue_specRestricted <- create_raincloudPlot(combined_data, sampled_data, "Tissue specific/restricted", colours[3])
plot_cell_specRestricted <- create_raincloudPlot(combined_data, sampled_data, "Cell type specific/restricted", colours[4])
plot_full <- create_raincloudPlot(hgnc_distinct, sample_all, "All", colours[5])

plot_housekeeping <- create_raincloudPlot(combined_distinct, sample_distinct, "Housekeeping", colours[1])
plot_ubiquitous_varying <- create_raincloudPlot(combined_distinct, sample_distinct, "Ubiquitous varying", colours[2])
plot_tissue_specRestricted <- create_raincloudPlot(combined_distinct, sample_distinct, "Tissue specific/restricted", colours[3])
plot_cell_specRestricted <- create_raincloudPlot(combined_distinct, sample_distinct, "Cell type specific/restricted", colours[4])

# Create a common y-axis title plot
y_title <- ggplot() + 
  annotate("text", x = 1, y = 0.5, label = "Nodes per Cluster", 
           angle = 90, size = 5) +
  theme_void() +
  coord_cartesian(clip = "off")

# Create a common x-axis title plot
x_title <- ggplot() + 
  annotate("text", x = 0.5, y = 0, label = "ZF Gene Category", 
           size = 5, vjust = 0.5) +
  theme_void() +
  coord_cartesian(clip = "off")

# Combine the individual plots into a single layout with common axis titles
combined_plot <- (y_title | (plot_housekeeping + plot_ubiquitous_varying + plot_tissue_specRestricted + plot_cell_specRestricted + plot_full) + 
                    plot_layout(nrow = 1)) +
  plot_layout(widths = c(0.05, 1))

final_plot <- combined_plot / x_title +
  plot_layout(heights = c(1, 0.05))

# Display the combined plot
print(final_plot)

ggsave("Raincloud No Border - Nodes per Cluster.png", plot = final_plot,
       width = 13.34, height = 7.5, dpi = 900, bg = "transparent")


# Raincloud Plot Human Genes ---------------------------------------------------
cancerGenes <- read_tsv(file.path(inpath, "Cosmic_CancerGeneCensus_v100_GRCh38.tsv.gz"))
# cancerGenes <- read_tsv(file.path(inpath, "Cosmic_CancerGeneCensus_v98_GRCh38.tsv.gz"))

# Filter for oncogenes (includes genes that are classified as oncogenes + TSG)
# Hyperactivity of the gene drives the transformation
oncogene <- dplyr::filter(cancerGenes, grepl("oncogene", ROLE_IN_CANCER))

# Filter for tumour suppressor genes (without oncogenes)
# Loss of gene function drives the transformation
TSG <- dplyr::filter(cancerGenes, grepl("TSG", ROLE_IN_CANCER))
TSG_filtered <- TSG %>% anti_join(oncogene, by = "GENE_SYMBOL")

# Fusion: the gene is known to be involved in oncogenic fusions
fusion <- dplyr::filter(cancerGenes, grepl("fusion", ROLE_IN_CANCER))

oncogene <- hgnc_orthologsDF %>%
  filter(human_symbol %in% oncogene$GENE_SYMBOL) %>%
  select(all_of(columns_hgnc)) %>%
  group_by(cluster) %>%
  mutate(cluster_count = n()) %>%
  ungroup()

oncogene_distinct <- oncogene %>% 
  dplyr::distinct(cluster, .keep_all = TRUE)

TSG_filtered <- hgnc_orthologsDF %>%
  filter(human_symbol %in% TSG$GENE_SYMBOL) %>%
  select(all_of(columns_hgnc)) %>%
  group_by(cluster) %>%
  mutate(cluster_count = n()) %>%
  ungroup()

tsg_distinct <- TSG_filtered %>% 
  dplyr::distinct(cluster, .keep_all = TRUE)

fusion <- hgnc_orthologsDF %>%
  filter(human_symbol %in% fusion$GENE_SYMBOL) %>%
  select(all_of(columns_hgnc)) %>%
  group_by(cluster) %>%
  mutate(cluster_count = n()) %>%
  ungroup()

fusion_distinct <- fusion %>% 
  dplyr::distinct(cluster, .keep_all = TRUE)
  
# Add a category column to each data frame
oncogene_distinct$category <- "Oncogene"
tsg_distinct$category <- "TSG"
fusion_distinct$category <- "Fusion"

# Combine the data frames
combined_human <- bind_rows(oncogene_distinct, tsg_distinct, fusion_distinct)

# Sample a subset of points to reduce overcrowding
sampled_data <- sample_points(combined_human, proportion = 0.3)

# Create individual plots for each category
plot_oncogene <- create_raincloudPlot(combined_human, sampled_data, "Oncogene", colours[1])
plot_tsg <- create_raincloudPlot(combined_human, sampled_data, "TSG", colours[2])
plot_fusion <- create_raincloudPlot(combined_human, sampled_data, "Fusion", colours[3])
plot_full <- create_raincloudPlot(hgnc_distinct, sample_all, "All", colours[5])

# Create a common y-axis title plot
y_title <- ggplot() + 
  annotate("text", x = 1, y = 0.5, label = "Node Count per Ortholog Cluster", 
           angle = 90, size = 5) +
  theme_void() +
  coord_cartesian(clip = "off")

# Create a common x-axis title plot
x_title <- ggplot() + 
  annotate("text", x = 0.5, y = 0, label = "Human Gene Category", 
           size = 5, vjust = 0.5) +
  theme_void() +
  coord_cartesian(clip = "off")

# Combine the individual plots into a single layout with common axis titles
combined_plot <- (y_title | (plot_oncogene + plot_tsg + plot_fusion + plot_full) + 
                    plot_layout(nrow = 1)) +
  plot_layout(widths = c(0.05, 1))

final_plot <- combined_plot / x_title +
  plot_layout(heights = c(1, 0.05))

# Display the combined plot
print(final_plot)

ggsave("Human Raincloud No Border - Nodes per Cluster.png", plot = final_plot,
       width = 13.34, height = 7.5, dpi = 900, bg = "transparent")


# ALLUVIAL PLOT HERE
# https://cran.r-project.org/web/packages/ggalluvial/vignettes/ggalluvial.html
# STUDENT CURRICULA 
data(majors)



# HCOP Gene Identifiers -------------------------------------------------------

# Replace NA with empty strings
hgnc_orthologsDF <- hgnc_orthologsDF %>%
  mutate(
    human_ensembl_gene = ifelse(is.na(human_ensembl_gene), "", human_ensembl_gene),
    human_symbol = ifelse(is.na(human_symbol), "", human_symbol),
    human_entrez_gene = ifelse(is.na(human_entrez_gene), "", human_entrez_gene),
    human_gene_id = ifelse(is.na(human_gene_id), "", human_gene_id),
    zebrafish_ensembl_gene = ifelse(is.na(zebrafish_ensembl_gene), "", zebrafish_ensembl_gene),
    zebrafish_symbol = ifelse(is.na(zebrafish_symbol), "", zebrafish_symbol),
    zebrafish_entrez_gene = ifelse(is.na(zebrafish_entrez_gene), "", zebrafish_entrez_gene),
    zfin_id = ifelse(is.na(zfin_id), "", zfin_id),
    zebrafish_gene_id = ifelse(is.na(zebrafish_gene_id), "", zebrafish_gene_id),
  )

# Create combined identifier for each row
hgnc_orthologsDF <- hgnc_orthologsDF %>%
  mutate(human_combination = paste0(
    ifelse(human_ensembl_gene != "", "Ensembl ID&", ""),
    ifelse(human_symbol != "", "HGNC Symbol&", ""),
    ifelse(human_entrez_gene != "", "Entrez ID&", ""),
    ifelse(human_gene_id != "", "Human ID&", "")
  ))

hgnc_orthologsDF <- hgnc_orthologsDF %>%
  mutate(human_combination = paste0(
    ifelse(human_ensembl_gene != "", "Ensembl ID&", ""),
    ifelse(human_symbol != "", "HGNC Symbol&", ""),
    ifelse(human_entrez_gene != "", "Entrez ID&", "")
  ))

# Remove trailing '&' and space
hgnc_orthologsDF <- hgnc_orthologsDF %>%
  mutate(human_combination = sub("& $", "", human_combination))

# Count number of intersections per combination
human_counts <- hgnc_orthologsDF %>%
  count(human_combination)

# Create a named vector for the UpSet plot input
human_upset_input <- setNames(human_counts$n, 
                              human_counts$human_combination)
names(human_upset_input) <- trimws(names(human_upset_input))

# Generate the UpSet plot for human gene identifiers
upset(fromExpression(human_upset_input),
      #sets = c("Human ID", "Ensembl ID", "HGNC Symbol", "Entrez ID"),
      sets = c("Ensembl ID", "HGNC Symbol", "Entrez ID"),
      #mb.ratio = c(0.65, 0.35),
      order.by = "freq",
      decreasing = TRUE,
      text.scale = 1.1,
      point.size = 2.8,
      set_size.show = FALSE,
      line.size = 1)

# Calculate ICC ---------------------------------------------------------------
HGNChuman_zebrafish <- read_delim("HGNChuman_zebrafishOrthologs.tsv",
                                  delim = "\t", escape_double = FALSE, trim_ws = TRUE)

HGNCzebrafish_human <- read_delim("HGNCzebrafish_humanOrthologs.tsv",
                                  delim = "\t", escape_double = FALSE, trim_ws = TRUE)

# Subset the relevant columns
human_zebrafishSubset <- HGNChuman_zebrafish %>% select(median_support, weighted_average)
zebrafish_humanSubset <- HGNCzebrafish_human %>% select(median_support, weighted_average)

# Calculate ICC
irr::icc(human_zebrafishSubset, model = "twoway", type = "agreement", unit = "average")
irr::icc(zebrafish_humanSubset, model = "twoway", type = "agreement", unit = "average")

icc_details <- icc(human_zebrafishSubset, model = "twoway", type = "agreement", unit = "average")
