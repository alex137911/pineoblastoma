# Remove objects in workspace
rm(list = ls())

# Required Packages -----------------------------------------------------------
suppressMessages(library(readr))
suppressMessages(library(readxl))
suppressMessages(library(dplyr))
# devtools::install_version("dbplyr", version = "2.3.4")
# https://stackoverflow.com/questions/77370659/
suppressMessages(library(dbplyr))
suppressMessages(library(data.table))
suppressMessages(library(ggplot2))
suppressMessages(library(forcats))
suppressMessages(library(colorspace))
suppressMessages(library(patchwork))
suppressMessages(library(tidyr))
suppressMessages(library(DESeq2))
suppressMessages(library(BiocManager))
suppressMessages(library(biomaRt))
suppressMessages(library(clusterProfiler))
suppressMessages(library(org.Dr.eg.db))  # For zebrafish gene annotations
suppressMessages(library(org.Hs.eg.db))  # For human gene annotations
suppressMessages(library(AnnotationDbi)) # Used in gene ontology (https://www.youtube.com/watch?v=JPwdqdo_tRg&ab_channel=Sanbomics)
suppressMessages(library(gprofiler2))
suppressMessages(library(fgsea))
suppressMessages(library(msigdbr))
# suppressMessages(library(rrvgo))        # Reduce redundancy in GO terms
suppressMessages(library(GO.db))
suppressMessages(library(GOplot))
suppressMessages(library(topGO))        # Reduce redundant GO terms
suppressMessages(library(fuzzyjoin))    # Join GO ID to GO terms w/o considering all edge cases (e.g., differences in capitalization)
suppressMessages(library(reshape2))
suppressMessages(library(ggsignif))
suppressMessages(library(MKinfer))
suppressMessages(library(ggrepel))      # For volcano plot
suppressMessages(library(pheatmap))
suppressMessages(library(stringr))
suppressMessages(library(DOSE)) # Dotplot tutorial (https://learn.gencore.bio.nyu.edu/rna-seq-analysis/gene-set-enrichment-analysis/)
suppressMessages(library(KEGGREST))
suppressMessages(library(createKEGGdb)) #https://github.com/YuLab-SMU/clusterProfiler/issues/561#issuecomment-1467266614
suppressMessages(library(KEGG.db))
suppressMessages(library(ggvenn))

# KEGG Database Creation -------------------------------------------------------
species <-c("ath","hsa","mmu", "rno","dre","dme","cel")
createKEGGdb::create_kegg_db(species)

# Load Data -------------------------------------------------------------------
inpath <- "C:/Users/acale/OneDrive/Documents/Waterloo BME/Co-op/4C/pineoblastoma/Data/Input/Laurie/"
inDir  <- sprintf("%s/Laurie", dirname(inpath))
setwd(inDir)

# Pineoblastoma tumours bulk RNA-seq data
# Fragments Per Kilobase of transcript per Million mapped reads (FPKM)
# FPKM is calculated by normalizing the raw read counts by both the length of 
# the gene and the total number of reads in the sample. This normalization accounts 
# for different gene lengths and sequencing depths.
pineoblastoma_tumours <- read_excel("pb_myc_vs_fetalbrain.xlsx")

pineoblastoma_tumoursDF <- as.data.frame(pineoblastoma_tumours)

pb_myc_fetalbrain_tpm <- read_delim("pb_myc_fetalbrain.tpm.txt", 
                                    delim = "\t", escape_double = FALSE, 
                                    trim_ws = TRUE)

pb_myc_fetalbrainDF <- as.data.frame(pb_myc_fetalbrain_tpm)

# Old (?) Zebrafish RNA-seq data
# zebrafish_all_20240508_norm <- read_delim("zebrafish.all.20240508.norm.fpkm", 
#                                           delim = "\t", escape_double = FALSE, 
#                                           trim_ws = TRUE)
# 
# zebrafish_all_20240508_normDF <- as.data.frame(zebrafish_all_20240508_norm)

# New directory
inpath <- "C:/Users/acale/OneDrive/Documents/Waterloo BME/Co-op/4C/pineoblastoma/Data/Input/Liming/"
inDir  <- sprintf("%s/Liming", dirname(inpath))
setwd(inDir)

# Zebrafish DESeq vs wildtype analysis
mycBrain_wildtype <- read_delim("all.counts.myc_brain_vs_wildtype_brain.DESeq2_all.xls",
                                delim = "\t", escape_double = FALSE, 
                                trim_ws = TRUE)

mycBrain_wildtypeDF <- as.data.frame(mycBrain_wildtype)

# Zebrafish DESeq vs control (p53 background loss) brain analysis
mycBrain_control <- read_delim("all.counts.myc_brain_vs_ctr_brain.DESeq2_all.xls",
                               delim = "\t", escape_double = FALSE, 
                               trim_ws = TRUE)

mycBrain_controlDF <- as.data.frame(mycBrain_control)

# Drop Ensembl ID version number
mycBrain_wildtypeDF$gene_id <- sub("\\..*", "", mycBrain_wildtypeDF$gene_id)
mycBrain_controlDF$gene_id  <- sub("\\..*", "", mycBrain_controlDF$gene_id)

outpath <- "C:/Users/acale/OneDrive/Documents/Waterloo BME/Co-op/4C/pineoblastoma/Data/Output/"
outDir  <- sprintf("%s/Output", dirname(outpath))
setwd(outDir)

# Read in Ortholog Data
HGNChuman_zebrafish <- read_delim("HGNChuman_zebrafishOrthologs-final.tsv", 
                                  delim = "\t", escape_double = FALSE, 
                                  trim_ws = TRUE)

HGNCzebrafish_human <- read_delim("HGNCzebrafish_humanOrthologs-final.tsv",
                                  delim = "\t", escape_double = FALSE,
                                  trim_ws = TRUE)

# Human PB-MYC DE Genes --------------------------------------------------------
# Apply fold change threshold (1.5-fold change from Yu et al.)
fc_threshold  <- 1.5
fdr_threshold <- 0.05

# 8,720/20,535 for wildtype
de_HSAgenes <- pineoblastoma_tumoursDF %>%
  filter(abs(log2FoldChange) >= log2(fc_threshold) & padj <= fdr_threshold)

# Log Transformed Expression Level of PB Marker Genes --------------------------

# Transform TPM values to log2 scale
log2_tpm <- log2(pb_myc_fetalbrainDF[,-1] + 1)           # Adding 1 to avoid log2(0)
log2_tpm$human_gene <- pb_myc_fetalbrainDF$human_gene

# Reshape data for plotting
melted_data <- melt(log2_tpm, id.vars = "human_gene")

# Create a sample annotation dataframe
sample_annotation <- data.frame(
  sample = colnames(log2_tpm)[-ncol(log2_tpm)],  # Exclude the last column (human_gene)
  group = c(rep("PB-MYC/FOXR2", 2), rep("Fetal Brain", 11))
)

# Merge annotation with melted data
melted_data <- merge(melted_data, sample_annotation, by.x = "variable", by.y = "sample")

# Rename columns for clarity
colnames(melted_data) <- c("sample", "human_gene", "log2_norm_counts", "group")

# Filter data for genes of interest
marker_genes <- c("CRX", "RB1", "FOXR2", "SYP", "MKI67")
filtered_data <- melted_data[melted_data$human_gene %in% marker_genes,]

# Perform Wilcoxon rank-sum test for each gene
p_values <- filtered_data %>%
  group_by(human_gene) %>%
  summarise(p_value = wilcox.test(log2_norm_counts ~ group)$p.value)

# Add significance levels to the data
filtered_data <- filtered_data %>%
  left_join(p_values, by = "human_gene") %>%
  mutate(significance = ifelse(p_value < 0.05, "*", "ns"))

# Nonparametric bootstrap t-test following (Dwivedi et al., 2017) - doi: 10.1002/sim.7263
# B = number of bootstrap samples
perform_bootstrap_t_test <- function(gene, data, B) {
  gene_data <- data[data$human_gene == gene, ]
  result <- boot.t.test(log2_norm_counts ~ group, data = gene_data, B = B)
  return(result)
}

# Get unique genes
unique_genes <- unique(filtered_data$human_gene)

# Initialize a list to store results
results_list <- list()

# Perform bootstrap t-test for each gene
for (gene in unique_genes) {
  print(paste("Performing bootstrap t-test for gene:", gene))
  result <- perform_bootstrap_t_test(gene, filtered_data, B = 100)
  results_list[[gene]] <- result
}

# Print results for each gene
for (gene in names(results_list)) {
  cat("\nResults for gene:", gene, "\n")
  print(results_list[[gene]])
}

# Plot
# Reorder the factor levels of 'human_gene'
filtered_data$human_gene <- factor(filtered_data$human_gene, levels = c("CRX", "MKI67", "SYP", "FOXR2",  "RB1"))

# Calculate the sample counts for each group
sample_counts <- table(sample_annotation$group)

# Create new group labels with sample counts
group_labels <- paste0(names(sample_counts), " (n = ", as.numeric(sample_counts), ")")

# Update the group variable with new labels
filtered_data$group <- factor(filtered_data$group, levels = names(sample_counts), labels = group_labels)

# With Wilcoxon rank-sum p-value (sample size too low??)
ggplot(filtered_data, aes(x = group, y = FPKM, fill = group)) +
  geom_boxplot() +
  facet_wrap(~human_gene, scales = "fixed", ncol = 12) +  # Set scales to fixed for a standard y-axis range
  scale_fill_manual(values = c("#89D6F4", "#E3F0C6")) +  # Custom fill colors
  theme_bw() +
  labs(y = expression("log"[2]*" Normalized TPM"), x = NULL, fill = "Group") +  # Remove x-axis labels
  theme(
    axis.text.x = element_blank(),  # Remove x-axis text
    axis.ticks.x = element_blank(),  # Remove x-axis ticks
    legend.position = c(0.37, 0.10),  # x-coordinate, y-coordinate
  ) +
  coord_cartesian(ylim = c(0, 100)) +
  geom_signif(
    comparisons = list(c("Fetal Brain (n = 11)", "PB-MYC/FOXR2 (n = 2)")),
    map_signif_level = TRUE,
    test = "wilcox.test"
  )

# TPM Expression Level of PB Marker Genes --------------------------------------

# Reshape data for plotting
pineoblastoma_markerExpression <- pb_myc_fetalbrain_tpm

# Transform TPM values to log2 scale
log2_tpm <- log2(pineoblastoma_markerExpression[,-1] + 1)           # Adding 1 to avoid log2(0)
log2_tpm$human_gene <- pineoblastoma_markerExpression$human_gene

# Reshape data for plotting
melted_data <- melt(log2_tpm, id.vars = "human_gene")

# Create a sample annotation dataframe
sample_annotation <- data.frame(
  sample = colnames(log2_tpm)[-ncol(log2_tpm)],  # Exclude the last column (human_gene)
  group = c(rep("PB-MYC/FOXR2", 2), rep("Fetal Brain", 11))
)

# Merge annotation with melted data
melted_data <- merge(melted_data, sample_annotation, by.x = "variable", by.y = "sample")

# Rename columns for clarity
colnames(melted_data) <- c("sample", "human_gene", "log2_FPKM", "group")

# Filter data for genes of interest
marker_genes <- c("CRX", "FOXR2", "MYC", "RB1",
                  "LIN28A", "SYP", "CHGA", "GFAP", "OLIG2", "MKI67")

filtered_data <- melted_data[melted_data$human_gene %in% marker_genes,]

# Perform Wilcoxon rank-sum test for each gene
p_values <- filtered_data %>%
  group_by(human_gene) %>%
  summarise(p_value = wilcox.test(log2_FPKM ~ group)$p.value)

# Plot
# Reorder the factor levels of 'human_gene'
filtered_data$human_gene <- factor(filtered_data$human_gene, 
                                   levels = c("CRX", "FOXR2", "MYC", "RB1", "LIN28A", "SYP", "CHGA", "GFAP", "OLIG2", "MKI67"))

# Calculate the sample counts for each group
sample_counts <- table(sample_annotation$group)

# Create new group labels with sample counts
group_labels <- paste0(names(sample_counts), " (n = ", as.numeric(sample_counts), ")")

# Update the group variable with new labels
filtered_data$group <- factor(filtered_data$group, levels = names(sample_counts), labels = group_labels)

# With Wilcoxon rank-sum p-value (sample size too low??)
ggplot(filtered_data, aes(x = group, y = log2_FPKM, fill = group)) +
  geom_boxplot() +
  facet_wrap(~human_gene, scales = "fixed", ncol = 10) +  # Set scales to fixed for a standard y-axis range
  scale_fill_manual(values = c("#89D6F4", "#E3F0C6")) +  # Custom fill colors
  theme_bw() +
  labs(y = expression("log"[2]*" Normalized FPKM"), x = NULL, fill = "Group") +  # Remove x-axis labels
  theme(
    axis.text.x = element_blank(),  # Remove x-axis text
    axis.ticks.x = element_blank(),  # Remove x-axis ticks
    legend.position = c(0.24, 0.76),  # x-coordinate, y-coordinate
  ) +
  coord_cartesian(ylim = c(0, 9)) +
  geom_signif(
    comparisons = list(c("Fetal Brain (n = 11)", "PB-MYC/FOXR2 (n = 2)")),
    map_signif_level = TRUE,
    test = "wilcox.test"
  )

# With raw data
ggplot(filtered_data, aes(x = group, y = log2_FPKM, fill = group, color = group)) +
  geom_boxplot(alpha = 0.6) +  # Boxplot with fill color and outline color with transparency
  geom_point(position = position_dodge(width = 0.75), size = 1.5, alpha = 0.3) +  # Points without jitter
  facet_wrap(~human_gene, scales = "fixed", ncol = 10) +  # Use gene_symbol for faceting
  scale_fill_manual(values = c("#5A9FCD", "#93C6C9")) +  # Custom fill colors
  scale_color_manual(values = c("#17374D", "#203E40")) +  # Custom colors for outline and points
  theme_bw() +
  labs(y = NULL, x = NULL, fill = NULL) +  # Update y-axis label and remove fill legend title
  theme(
    axis.text.x = element_blank(),  # Remove x-axis text
    axis.ticks.x = element_blank(),  # Remove x-axis ticks
    legend.position = "top",  # Position legend at the bottom
    legend.title = element_blank(),  # Remove legend title
    legend.text = element_text(size = 10),  # Set legend text size
    legend.direction = "horizontal",  # Arrange legend items horizontally
    legend.box = "horizontal",  # Ensure legend items are in a single row
    plot.margin = margin(t = 10, r = 10, b = 10, l = 10)  # Adjust plot margins to reduce white space
  ) +
  scale_y_continuous(breaks = seq(2, 8, by = 2), limits = c(0, 9)) +  # Set y-axis marks
  geom_signif(
    comparisons = list(c("Fetal Brain (n = 11)", "PB-MYC/FOXR2 (n = 2)")),
    map_signif_level = TRUE,
    test = "wilcox.test",
    color = "black"  # Set text color for significance annotations
  ) +
  guides(fill = guide_legend(override.aes = list(alpha = 0.6)), color = "none")  # Remove color legend


ggsave("HSA-PB Marker Genes Expression.png", plot = plot, device = "png", 
       width = 23.01, height = 14.6, units = "cm", dpi = 900)

# Human Volcano Plot of DE Genes -----------------------------------------------
sum(pineoblastoma_tumoursDF$padj < 0.05) #9,409/20,535 

# Create a colour column based on the up or downregulated status
# Gray = gene does not mean DE threshold (1.5 foldchange and/or FDR < 0.05)

pineoblastoma_tumoursDF <- pineoblastoma_tumoursDF %>%
  mutate(color = case_when(
    (log2FoldChange >  1.5 & padj < 0.05) ~ "red",
    (log2FoldChange < -1.5 & padj < 0.05) ~ "blue",
    TRUE ~ "gray"
  ))

# Set thresholds for labeling
log2fc_threshold <- 6  # Threshold for log2 fold change
qvalue_threshold <- 25  # Threshold for -log10(p-value)

# Create a new column to indicate whether the gene is a marker gene
marker_genes <- c("CRX", "FOXR2", "MYC", "RB1", 
                  "LIN28A", "SYP", "CHGA", "GFAP", "OLIG2", "MKI67")

pineoblastoma_tumoursDF$marker <- ifelse(pineoblastoma_tumoursDF$gene %in% marker_genes, 
                                         "marker", "non-marker")

# Create the volcano plot
ggplot(pineoblastoma_tumoursDF, aes(x = log2FoldChange, y = -log10(padj), col = color, label = gene)) +
  geom_vline(xintercept = c(-1.5, 1.5), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') +
  geom_point(aes(alpha = marker), size = 2) +
  scale_color_manual(values = c("#00AFBB", "gray", "#bb0c00"),
                     labels = c("Downregulated", "Not significant", "Upregulated")) +
  scale_alpha_manual(values = c("marker" = 1, "non-marker" = 0.1)) +  # Adjust alpha values
  coord_cartesian(ylim = c(0, 40), xlim = c(-12, 30)) +
  labs(color = 'Expression',
       x = expression("log"[2]*"FC"),
       y = expression("-log"[10]*"q-value")) +
  scale_x_continuous(breaks = seq(-12, 30, 2)) +
  # geom_text_repel(aes(label = ifelse(-log10(padj) > qvalue_threshold & abs(log2FoldChange) > log2fc_threshold, gene, "")),
  #                 size = 3, max.overlaps = Inf) +
  geom_text_repel(aes(label = ifelse(gene %in% marker_genes, gene, "")),
                  size = 3, max.overlaps = Inf, segment.color = "grey30") +
  theme_minimal() +
  theme(axis.line = element_line(color = "gray10"),
        axis.title.x = element_text(margin = margin(t = 7)),    # Increase space between x-axis label and axis
        axis.title.y = element_text(margin = margin(r = 10)),   # Increase space between y-axis label and axis
        legend.position = c(0.90, 0.75)) +
  annotate("text", x = 1.5, y = 37, label = "FC = 1.5", size = 3, color = "black", vjust = -1) +                 # Label for right vertical line
  annotate("text", x = 23, y = -log10(0.15), label = "q-value = 0.05", size = 3, color = "black", hjust = -0.1)  # Label for horizontal line

# Human GO Enrichment Analysis -------------------------------------------------
human_gseaData <- as.data.frame(de_HSAgenes)

# Get DE genes as vector
hsa_deGenes <- human_gseaData$gene
hsa_backgroundGenes <- pineoblastoma_tumoursDF$gene

hsaGO_results <- enrichGO(gene = hsa_deGenes, universe = hsa_backgroundGenes, OrgDb = org.Hs.eg.db, 
                          keyType = "SYMBOL", ont = "BP", 
                          pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.10)

hsaGO_df <- as.data.frame(hsaGO_results)

# data(EC)
# david <- as.data.frame(EC$david)
# head(EC$genelist)

# Preparing the data for the GOplot
hsaGO_bubble <- hsaGO_df %>%
  mutate(Category = "BP") %>%
  # Renaming columns to match the desired format
  rename(Term = Description, adj_pval = p.adjust, genes = geneID) %>%
  dplyr::select(Category, ID, Term, genes, adj_pval) %>%
  mutate(genes = gsub("/", ", ", genes))

hsa_geneList <- human_gseaData %>%
  rename(ID = gene, logFC = log2FoldChange, AveExpr = baseMeanA, 
         P.value = pvalue, adj.P.Val = padj) %>%
  dplyr::select(ID, logFC, AveExpr, P.value, adj.P.Val)

hsaGO_bubble <- circle_dat(hsaGO_bubble, hsa_geneList)

# Bubble plot
GOBubble(hsaGO_bubble, labels = 1)

# Reduced plot with redundant terms (e.g. gene overlap > 0.75)
hsaGO_reduced <- reduce_overlap(hsaGO_bubble, overlap = 0.99)

GOBubble(hsaGO_reduced, labels = 1) 

# ggplot2
ggplot(hsaGO_bubble, aes(x = zscore, y = -log10(adj_pval))) +
  geom_point(aes(size = count, color = category), shape = 21, fill = "lightgreen") +
  geom_text(aes(label = ID), hjust = 1.1, vjust = 1.1, size = 3) +
  # scale_size_area(max_size = 10) +
  scale_size_continuous(range = c(2, 15)) + 
  scale_color_manual(values = c("BP" = "darkgreen")) +
  labs(x = "z-score", y = "-log(adj p-value)", size = "Count") +
  theme_minimal() +
  theme(legend.position = "bottom") +
  coord_cartesian(
    xlim = c(0, 10),
    ylim = c(0, max(-log10(hsaGO_bubble$adj_pval)))
  ) +
  guides(color = guide_legend(override.aes = list(size = 6)))

# Human Gene Set Enrichment Analysis -------------------------------------------
# Create ranks based on fold change
human_gseaRanks <- pineoblastoma_tumoursDF$log2FoldChange
names(human_gseaRanks) <- pineoblastoma_tumoursDF$gene

# Filter out non-finite values from human_gseaRanks (otherwise fgsea throws error)
human_gseaRanks <- human_gseaRanks[is.finite(human_gseaRanks)]

# Plot ranks
barplot(sort(human_gseaRanks, decreasing = TRUE))

# Load gene sets from MSigDB
# GO biological processes gene sets (following Anthony Liu and Bryan Li paper: 10.1007/s00401-021-02284-5)
goBiological_gsets <- msigdbr::msigdbr(species = "Homo sapiens", category = "C5", subcategory = "BP") %>%
  dplyr::select(gs_name, gene_symbol) %>%
  split(x = .$gene_symbol, f = .$gs_name)

# Load hallmark gene sets
hallmark_gsets <- msigdbr(species = "Homo sapiens", category = "H") %>%
  dplyr::select(gs_name, gene_symbol) %>%
  split(x = .$gene_symbol, f = .$gs_name)

# Load KEGG gene sets from MSigDB
kegg_gsets <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG") %>%
  dplyr::select(gs_name, gene_symbol) %>%
  split(x = .$gene_symbol, f = .$gs_name)

# Load Reactome gene sets from MSigDB
reactome_gsets <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME") %>%
  dplyr::select(gs_name, gene_symbol) %>%
  split(x = .$gene_symbol, f = .$gs_name)

# Combine all gene sets
combined_gsets <- c(hallmark_gsets, kegg_gsets, reactome_gsets)

# Run fgsea
human_fgsea <- fgsea(pathways = goBiological_gsets, stats = human_gseaRanks, 
                     minSize = 15, maxSize = 500, nperm = 1000)

human_fgseaCombined <- fgsea(pathways = combined_gsets, stats = human_gseaRanks, 
                             minSize = 15, maxSize = 500, nperm = 1000)

# Sort the results by adjusted p-value
human_fgsea <- human_fgsea[order(human_fgsea$padj), ]
human_fgseaCombined <- human_fgseaCombined[order(human_fgseaCombined$padj)]

head(human_fgsea[order(padj, -abs(NES)), ], n=10)
head(human_fgseaCombined[order(padj, -abs(NES)), ], n=10)

# Convert fgsea results to data frame
human_fgseaDF <- as.data.frame(human_fgsea)
human_fgseaCombinedDF <- as.data.frame(human_fgseaCombined)

# Reduce redundancy of GO terms
# Extract all GO terms
goterms <- Term(GOTERM)

# Map the GO term names to GO IDs
go_mapping <- data.frame(go_id = names(goterms), pathway = goterms, stringsAsFactors = FALSE)

go_mapping <- go_mapping %>%
  mutate(pathway = tolower(gsub("GOBP_", "", pathway)), pathway = gsub("-", " ", pathway),
         pathway = gsub(",", " ", pathway), pathway = gsub("  ", " ", pathway))


# Modify the pathway names
human_fgseaDF <- human_fgseaDF %>%
  mutate(pathway = tolower(gsub("GOBP_", "", pathway)), pathway = gsub("_", " ", pathway),
         pathway = gsub("-", " ", pathway))

# Ensure the pathways in human_fgseaDF are mapped correctly to GO IDs
human_fgseaDF <- fuzzyjoin::stringdist_left_join(human_fgseaDF, 
                                      go_mapping, by = "pathway", max_dist = 3)

human_fgseaDF <- human_fgseaDF %>%
  rename(pathway = pathway.x, go_mappingPathway = pathway.y)

# Extract the mapped GO term IDs and NES scores
go_terms <- human_fgseaDF$go_id
scores <- human_fgseaDF$NES

# Calculate semantic similarity
# simMatrix <- calculateSimMatrix(go_terms, orgdb="org.Hs.eg.db", ont="BP", method="Rel")

# Filter valid terms present in the similarity matrix
valid_terms <- intersect(go_terms, rownames(simMatrix))
valid_scores <- scores[go_terms %in% valid_terms]
names(valid_scores) <- valid_terms  # Ensure scores are named to match terms in simMatrix
simMatrix_filtered <- simMatrix[valid_terms, valid_terms]

# Reduce redundancy
reducedTerms <- reduceSimMatrix(simMatrix_filtered, valid_scores, threshold=0.7, 
                                orgdb="org.Hs.eg.db")

# Filter the results
topPathwaysUp <- human_fgsea[order(human_fgsea$NES, decreasing = TRUE), ][1:10, ]
topPathwaysDown <- human_fgsea[order(human_fgsea$NES, decreasing = FALSE), ][1:10, ]
topPathways <- rbind(topPathwaysUp, topPathwaysDown)

# Prepare the data for plotting
topPathways <- topPathways %>%
  mutate(Pathway = factor(pathway, levels = pathway[order(NES)]))

# Add color for positive and negative NES
topPathways$color <- ifelse(topPathways$NES > 0, "Upregulated", "Downregulated")

# Plot the data
ggplot(topPathways, aes(x = NES, y = Pathway, fill = color)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("Upregulated" = "red", "Downregulated" = "blue")) +
  labs(x = "Normalized Enrichment Score",
       y = "Pathway") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 7, face = "bold"),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_text(size = 12))

# Reduced GO plot
# Create unique factor levels based on both pathway and go_id
reduced_fgseaDF <- human_fgseaDF %>%
  filter(go_id %in% reducedTerms$parent) %>%
  mutate(Pathway = factor(paste(pathway, go_id, sep = "_"), levels = paste(pathway[order(NES)], go_id[order(NES)], sep = "_")),
         color = ifelse(NES > 0, "Upregulated", "Downregulated"))

# Top pathways based on NES
top_upregulated <- reduced_fgseaDF %>% filter(NES > 0) %>% arrange(desc(NES)) %>%
  slice(1:10)

top_downregulated <- reduced_fgseaDF %>% filter(NES < 0) %>% arrange(NES) %>%
  slice(1:10)

topPathways <- bind_rows(top_upregulated, top_downregulated)

# Combine top pathways
top_pathways <- bind_rows(top_upregulated, top_downregulated)

# Ensure unique factor levels for plotting
top_pathways <- top_pathways %>%
  mutate(Pathway = factor(paste(pathway, go_id, sep = "_"), levels = paste(pathway[order(NES)], go_id[order(NES)], sep = "_")),
         color = ifelse(NES > 0, "Upregulated", "Downregulated"))

ggplot(top_pathways, aes(x = NES, y = Pathway, fill = color)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("Upregulated" = "red", "Downregulated" = "blue")) +
  labs(x = "Normalized Enrichment Score",
       y = "Pathway") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 7, face = "bold"),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_text(size = 12))

# Zebrafish DE Genes -----------------------------------------------------------

# Apply fold change threshold (1.5-fold change from Yu et al.)
fc_threshold  <- 1.5
fdr_threshold <- 0.05

# 5,470/25,979 for wildtype
de_ZFgenes <- mycBrain_wildtypeDF %>%
  filter(abs(log2FC) >= log2(fc_threshold) & FDR <= fdr_threshold)

# 5,160/25,808 for control
de_ZFgenes <- mycBrain_controlDF %>%
  filter(abs(log2FC) >= log2(fc_threshold) & FDR <= fdr_threshold)

# Get gene IDs of the differentially expressed genes
gene_ids <- de_ZFgenes$gene_id

# # Save results
# write.csv(de_ZFgenes, file = "mycBrain_wildtypeDE.csv")

# TPM Expression Level of ZF PB Marker Genes -----------------------------------

# Extract only the count data columns
count_data <- mycBrain_wildtypeDF[, c("MYC3_Brain", "MYC4_Brain", "MYC5_Brain", "myc1_brain", "myc2_brain",
                                      "WT_Brain2", "WT_6_B", "WT_7_B", "WT_8_B")]

rownames(count_data) <- mycBrain_wildtypeDF$gene_id

# Convert FPKM to TPM
fpkm_sums <- colSums(count_data)
tpm_data <- sweep(count_data, 2, fpkm_sums, FUN = "/") * 1e6

# Transform to log2 scale
log2_counts <- log2(count_data + 1)  # Adding 1 to avoid log2(0)
log2_counts <- as.data.frame(log2_counts)
log2_counts$gene_id <- rownames(log2_counts)

# log2_tpm <- log2(tpm_data + 1)  # Adding 1 to avoid log2(0)
# log2_tpm <- as.data.frame(log2_tpm)
# log2_tpm$gene_id <- rownames(log2_tpm)


# Reshape data for plotting
melted_data <- melt(log2_counts, id.vars = "gene_id")

# Create a sample annotation dataframe
sample_annotation <- data.frame(
  sample = colnames(log2_counts)[-ncol(log2_counts)],  # Exclude the last column (gene_id)
  group = c(rep("PB-MYC/FOXR2", 5), rep("Wildtype", 4))
)

# Merge annotation with melted data
melted_data <- merge(melted_data, sample_annotation, by.x = "variable", by.y = "sample")

# Rename columns for clarity
colnames(melted_data) <- c("sample", "gene_id", "log2_norm_counts", "group")

# Filter data for genes of interest
gene_mapping <- data.frame(
  gene_symbol = c("crx", "foxr1", "myca", "mycb", "mych", "rb1",
                  "lin28aa", "lin28ab", "sypa", "sypb", "chga",
                  "gfap", "olig2", "mki67"),
  ensembl_id = c("ENSDARG00000011989", "ENSDARG00000004864", "ENSDARG00000045695", 
                 "ENSDARG00000007241", "ENSDARG00000077473", "ENSDARG00000006782",
                 "ENSDARG00000004328", "ENSDARG00000016999", "ENSDARG00000110528", 
                 "ENSDARG00000002230", "ENSDARG00000008829", "ENSDARG00000025301",
                 "ENSDARG00000040946", "ENSDARG00000091150")
)

zf_filtered_data <- melted_data[melted_data$gene_id %in% gene_mapping$ensembl_id,]

# Merge the gene mapping with the filtered data to add gene symbols
zf_filtered_data <- merge(zf_filtered_data, gene_mapping, by.x = "gene_id", by.y = "ensembl_id")
colnames(zf_filtered_data) <- c("ensembl_id", "sample", "log2_FPKM", "group", "gene_symbol")

# Perform Wilcoxon rank-sum test for each gene
p_values <- zf_filtered_data %>%
  group_by(gene_symbol) %>%
  summarise(p_value = wilcox.test(log2_FPKM ~ group)$p.value)

# Calculate the sample counts for each group
sample_counts <- table(sample_annotation$group)

# Update the group variable with new labels
group_labels <- paste0(names(sample_counts), " (n = ", as.numeric(sample_counts), ")")
zf_filtered_data$group <- factor(zf_filtered_data$group, levels = names(sample_counts), labels = group_labels)

# Reorder the factor levels
zf_filtered_data$gene_symbol <- factor(zf_filtered_data$gene_symbol, levels = gene_mapping$gene_symbol)
zf_filtered_data$group <- factor(zf_filtered_data$group, levels = c("Wildtype (n = 4)", "PB-MYC/FOXR2 (n = 5)"))

# Plot
ggplot(zf_filtered_data, aes(x = group, y = log2_norm_counts, fill = group)) +
  geom_boxplot() +
  facet_wrap(~gene_symbol, scales = "fixed", ncol = 15) +  # Use gene_symbol for faceting
  scale_fill_manual(values = c("#89D6F4", "#E3F0C6")) +  # Custom fill colors
  theme_bw() +
    labs(y = "log2 Normalized TPM", x = NULL, fill = "Group") +  # Update y-axis label
  theme(
    axis.text.x = element_blank(),  # Remove x-axis text
    axis.ticks.x = element_blank(),  # Remove x-axis ticks
    legend.position = c(0.30, 0.77),  # x-coordinate, y-coordinate
  ) +
  coord_cartesian(ylim = c(0, 8)) +
  geom_signif(
    comparisons = list(c("Wildtype (n = 4)", "PB-MYC/FOXR2 (n = 5)")),
    map_signif_level = TRUE,
    test = "wilcox.test"
  )

# With raw data
ggplot(zf_filtered_data, aes(x = group, y = log2_FPKM, fill = group, color = group)) +
  geom_boxplot(alpha = 0.6) +  # Boxplot with fill color and outline color with transparency
  geom_point(position = position_dodge(width = 0.75), size = 1.5, alpha = 0.3) +  # Points without jitter
  facet_wrap(~gene_symbol, scales = "fixed", ncol = 15) +  # Use gene_symbol for faceting
  scale_fill_manual(values = c("#A9D3C0", "#C7E2A9")) +  # Custom fill colors
  scale_color_manual(values = c("#2B5341", "#334B19")) +  # Custom colors for outline and points
  theme_bw() +
  labs(y = NULL, x = NULL, fill = NULL) +  # Update y-axis label and remove fill legend title
  theme(
    axis.text.x = element_blank(),  # Remove x-axis text
    axis.ticks.x = element_blank(),  # Remove x-axis ticks
    legend.position = "bottom",  # Position legend at the bottom
    legend.title = element_blank(),  # Remove legend title
    legend.text = element_text(size = 10),  # Set legend text size
    legend.direction = "horizontal",  # Arrange legend items horizontally
    legend.box = "horizontal",  # Ensure legend items are in a single row
    plot.margin = margin(t = 5, r = 10, b = 5, l = 10)  # Adjust plot margins to reduce white space
  ) +
  coord_cartesian(ylim = c(0, 8)) +
  geom_signif(
    comparisons = list(c("Wildtype (n = 4)", "PB-MYC/FOXR2 (n = 5)")),
    map_signif_level = TRUE,
    test = "wilcox.test",
    color = "black"  # Set text color for significance annotations
  ) +
  guides(fill = guide_legend(override.aes = list(alpha = 0.6)), color = "none")  # Remove color legend

# Save the plot with the specified size in cm
ggsave("ZF-PB Marker Genes Expression.png", plot = plot, device = "png", 
       width = 23.01, height = 14.6, units = "cm", dpi = 900)

# TPM ZF + HSA Marker Gene Expression Level for Poster -------------------------

# Create individual plots for human and zebrafish
create_individual_plot <- function(data, genes, fill_colors, outline_colors, gene_column, comparisons, show_y_axis, signif_y_position) {
  plots <- list()
  for (i in seq_along(genes)) {
    gene <- genes[i]
    subset_data <- subset(data, get(gene_column) == gene)
    if (nrow(subset_data) > 0) {
      p <- ggplot(subset_data, aes(x = group, y = log2_FPKM, fill = group, color = group)) +
        geom_boxplot(alpha = 0.6) +
        geom_point(position = position_dodge(width = 0.75), size = 1.5, alpha = 0.5) +
        scale_fill_manual(values = fill_colors) +
        scale_color_manual(values = outline_colors) +
        theme_bw() +
        labs(y = NULL, x = NULL, fill = NULL) +
        theme(
          axis.text.y = if (show_y_axis && i == 1) element_text() else element_blank(),  # Show y-axis text only for the first plot
          axis.ticks.y = if (show_y_axis && i == 1) element_line() else element_blank(),  # Show y-axis ticks only for the first plot
          axis.text.x = element_blank(),  # Remove x-axis text
          axis.ticks.x = element_blank(),  # Remove x-axis ticks
          legend.position = "none",  # Remove legend from individual plots
          plot.margin = margin(t = 10, r = 2, b = 10, l = 2)
        ) +
        scale_y_continuous(breaks = seq(2, 8, by = 2), limits = c(0, 10)) +
        geom_signif(
          comparisons = comparisons,
          map_signif_level = TRUE,
          test = "wilcox.test",
          color = "black",
          y_position = signif_y_position
        )
      plots[[gene]] <- p
    }
  }
  plots
}

# Define the gene groups and the corresponding data
top_genes <- c("CRX", "FOXR2", "MYC", "RB1", "LIN28A", "SYP", "CHGA", "GFAP", "OLIG2", "MKI67")
bottom_genes <- c("crx", "foxr1", "myca", "mycb", "mych", "rb1", "lin28ab", "sypa", "sypb", "chga", "gfap", "olig2", "mki67")

# Define the group names for comparison
top_comparisons <- list(c("Fetal Brain (n = 11)", "PB-MYC/FOXR2 (n = 2)"))
bottom_comparisons <- list(c("Wildtype (n = 4)", "PB-MYC/FOXR2 (n = 5)"))

# Define the fixed y-position for significance annotations
signif_y_position <- 8.8

# Create individual plots
top_plots <- create_individual_plot(filtered_data, top_genes, c("#5A9FCD", "#93C6C9"), c("#17374D", "#203E40"), "human_gene", top_comparisons, TRUE, signif_y_position)
bottom_plots <- create_individual_plot(zf_filtered_data, bottom_genes, c("#A9D3C0", "#C7E2A9"), c("#2B5341", "#334B19"), "gene_symbol", bottom_comparisons, TRUE, signif_y_position)

# Combine the plots with patchwork specifying widths
top_row <- wrap_plots(
  top_plots$CRX, top_plots$FOXR2, top_plots$MYC, top_plots$RB1, top_plots$LIN28A, top_plots$SYP, top_plots$CHGA, top_plots$GFAP, top_plots$OLIG2, top_plots$MKI67,
  widths = c(1, 1, 3, 1, 1, 2, 1, 1, 1, 1)
)

bottom_row <- wrap_plots(
  bottom_plots$crx, bottom_plots$foxr1, (bottom_plots$myca | bottom_plots$mycb | bottom_plots$mych), bottom_plots$rb1, bottom_plots$lin28ab, (bottom_plots$sypa | bottom_plots$sypb),
  bottom_plots$chga, bottom_plots$gfap, bottom_plots$olig2, bottom_plots$mki67,
  widths = c(1, 1, 3, 1, 1, 2, 1, 1, 1, 1)
)

final_plot <- (top_row) / (bottom_row) +
  plot_layout(heights = c(1, 1))

# Print the final combined plot
print(final_plot)

# Save the plot with the specified size in cm
ggsave("ZF + HSA PB Marker Genes Expression.png", plot = final_plot, device = "png", 
       width = 23.01, height = 29.2, units = "cm", dpi = 900)

# Expression Level of ZF PB Marker Genes vs Control ----------------------------

# Extract only the count data columns
count_data <- mycBrain_controlDF[, c("MYC3_Brain", "MYC4_Brain", "MYC5_Brain", "myc1_brain", "myc2_brain",
                                      "Ctr_5B", "Ctr_6_B", "Ctr_7_B")]

rownames(count_data) <- mycBrain_controlDF$gene_id

# Transform to log2 scale
log2_counts <- log2(count_data + 1)  # Adding 1 to avoid log2(0)
log2_counts <- as.data.frame(log2_counts)
log2_counts$gene_id <- rownames(log2_counts)

# Reshape data for plotting
melted_data <- melt(log2_counts, id.vars = "gene_id")

# Create a sample annotation dataframe
sample_annotation <- data.frame(
  sample = colnames(log2_counts)[-ncol(log2_counts)],  # Exclude the last column (gene_id)
  group = c(rep("PB-MYC/FOXR2", 5), rep("Control", 3))
)

# Merge annotation with melted data
melted_data <- merge(melted_data, sample_annotation, by.x = "variable", by.y = "sample")

# Rename columns for clarity
colnames(melted_data) <- c("sample", "gene_id", "log2_norm_counts", "group")

# Filter data for genes of interest
gene_mapping <- data.frame(
  gene_symbol = c("crx", "mki67", "sypa", "sypb", "chga", "gfap",
                  "foxr1", "rb1", "myca", "mycb", "mych" ),
  ensembl_id = c("ENSDARG00000011989", "ENSDARG00000091150", "ENSDARG00000110528", 
                 "ENSDARG00000002230", "ENSDARG00000008829", "ENSDARG00000025301",
                 "ENSDARG00000004864", "ENSDARG00000006782", "ENSDARG00000045695", 
                 "ENSDARG00000007241", "ENSDARG00000077473")
)

filtered_data <- melted_data[melted_data$gene_id %in% gene_mapping$ensembl_id,]

# Merge the gene mapping with the filtered data to add gene symbols
filtered_data <- merge(filtered_data, gene_mapping, by.x = "gene_id", by.y = "ensembl_id")
colnames(filtered_data) <- c("ensembl_id", "sample", "log2_norm_counts", "group", "gene_symbol")

# Perform Wilcoxon rank-sum test for each gene
p_values <- filtered_data %>%
  group_by(gene_symbol) %>%
  summarise(p_value = wilcox.test(log2_norm_counts ~ group)$p.value)

# Calculate the sample counts for each group
sample_counts <- table(sample_annotation$group)

# Update the group variable with new labels
group_labels <- paste0(names(sample_counts), " (n = ", as.numeric(sample_counts), ")")
filtered_data$group <- factor(filtered_data$group, levels = names(sample_counts), labels = group_labels)

# Reorder the factor levels
filtered_data$gene_symbol <- factor(filtered_data$gene_symbol, levels = gene_mapping$gene_symbol)
filtered_data$group <- factor(filtered_data$group, levels = c("Control (n = 3)", "PB-MYC/FOXR2 (n = 5)"))

# Plot
ggplot(filtered_data, aes(x = group, y = log2_norm_counts, fill = group)) +
  geom_boxplot() +
  facet_wrap(~gene_symbol, scales = "fixed", ncol = 11) +  # Use gene_symbol for faceting
  scale_fill_manual(values = c("#89D6F4", "#E3F0C6")) +  # Custom fill colors
  theme_bw() +
  labs(y = "log2 Normalized Counts", x = NULL, fill = "Group") +  # Update y-axis label
  theme(
    axis.text.x = element_blank(),  # Remove x-axis text
    axis.ticks.x = element_blank(),  # Remove x-axis ticks
    legend.position = c(0.85, 0.75),  # x-coordinate, y-coordinate
  ) +
  coord_cartesian(ylim = c(0, 8)) +
  geom_signif(
    comparisons = list(c("Control (n = 3)", "PB-MYC/FOXR2 (n = 5)")),
    map_signif_level = TRUE,
    test = "wilcox.test"
  )

# Zebrafish Volcano Plot of DE genes ------------------------------------------
# Wildtype
sum(mycBrain_wildtypeDF$FDR < 0.05) #6,546/25,979 

# Create a colour column based on the up or downregulated status
# Gray = gene does not mean DE threshold (1.5 foldchange and/or FDR < 0.05)
mycBrain_wildtypeDF <- mycBrain_wildtypeDF %>%
  mutate(color = case_when(
    updown == "UP" ~ "red",
    updown == "DOWN" ~ "blue",
    TRUE ~ "gray"
  ))

# Set thresholds for labeling
log2fc_threshold <- 6  # Threshold for log2 fold change
qvalue_threshold <- 12  # Threshold for -log10(p-value)

# Create a new column to indicate whether the gene is a marker gene
gene_mapping <- data.frame(
  gene_symbol = c("crx", "foxr1", "myca", "mycb", "mych", "rb1",
                  "lin28aa", "lin28ab", "sypa", "sypb", "chga",
                  "gfap", "olig2", "mki67"),
  ensembl_id = c("ENSDARG00000011989", "ENSDARG00000004864", "ENSDARG00000045695", 
                 "ENSDARG00000007241", "ENSDARG00000077473", "ENSDARG00000006782",
                 "ENSDARG00000004328", "ENSDARG00000016999", "ENSDARG00000110528", 
                 "ENSDARG00000002230", "ENSDARG00000008829", "ENSDARG00000025301",
                 "ENSDARG00000040946", "ENSDARG00000091150")
)

mycBrain_wildtypeDF <- mycBrain_wildtypeDF %>%
  left_join(gene_mapping, by = c("gene_id" = "ensembl_id"))

mycBrain_wildtypeDF$marker <- ifelse(mycBrain_wildtypeDF$gene_id %in% gene_mapping$ensembl_id, 
                                     "marker", "non-marker")

mycBrain_wildtypeDFtest <- mycBrain_wildtypeDF
mycBrain_wildtypeDFtest <- mycBrain_wildtypeDFtest[mycBrain_wildtypeDFtest$gene_id %in% gene_mapping$ensembl_id,]

# Create the volcano plot
ggplot(mycBrain_wildtypeDF, aes(x = log2FC, y = -log10(FDR), col = color, label = gene_symbol)) +
  geom_vline(xintercept = c(-1.5, 1.5), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') +
  geom_point(aes(alpha = marker), size = 2) +
  scale_color_manual(values = c("#00AFBB", "gray", "#bb0c00"),
                     labels = c("Downregulated", "Not significant", "Upregulated")) +
  scale_alpha_manual(values = c("marker" = 1, "non-marker" = 0.1)) +  # Adjust alpha values
  coord_cartesian(ylim = c(0, 30), xlim = c(-10, 10)) +
  labs(color = 'Expression',
       x = expression("log"[2]*"FC"),
       y = expression("-log"[10]*"q-value")) +
  scale_x_continuous(breaks = seq(-10, 10, 2)) +
  geom_text_repel(aes(label = ifelse(gene_symbol %in% gene_mapping$gene_symbol, 
                                     gene_symbol, "")),
                  size = 3, max.overlaps = Inf, segment.color = "grey30") +
  theme_minimal() +
  theme(axis.line = element_line(color = "gray10"),
        axis.title.x = element_text(margin = margin(t = 7)),    # Increase space between x-axis label and axis
        axis.title.y = element_text(margin = margin(r = 10)),   # Increase space between y-axis label and axis
        legend.position = c(0.90, 0.75)) +
  annotate("text", x = 1.5, y = 27, label = "FC = 1.5", size = 3, color = "black", vjust = -1) +                 # Label for right vertical line
  annotate("text", x = -11, y = -log10(0.15), label = "q-value = 0.05", size = 3, color = "black", hjust = -0.1)  # Label for horizontal line


# Control
sum(mycBrain_controlDF$FDR < 0.05) #6,265/25,808 

# Create a colour column based on the up or downregulated status
# Gray = gene does not mean DE threshold (1.5 foldchange and/or FDR < 0.05)
mycBrain_controlDF <- mycBrain_controlDF %>%
  mutate(color = case_when(
    updown == "UP" ~ "red",
    updown == "DOWN" ~ "blue",
    TRUE ~ "gray"
  ))

# Set thresholds for labeling
log2fc_threshold <- 6  # Threshold for log2 fold change
qvalue_threshold <- 12  # Threshold for -log10(p-value)

# Create the volcano plot
ggplot(mycBrain_controlDF, aes(x = log2FC, y = -log10(FDR), col = color, label = gene_id)) +
  geom_vline(xintercept = c(-1.5, 1.5), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') +
  geom_point(size = 2) +
  scale_color_manual(values = c("#00AFBB", "gray", "#bb0c00"),
                     labels = c("Downregulated", "Not significant", "Upregulated")) +
  coord_cartesian(ylim = c(0, 20), xlim = c(-10, 10)) +
  labs(color = 'Expression',
       x = expression("log"[2]*"FC"),
       y = expression("-log"[10]*"q-value")) +
  scale_x_continuous(breaks = seq(-10, 10, 2)) +
  geom_text_repel(aes(label = ifelse(-log10(FDR) > qvalue_threshold & abs(log2FC) > log2fc_threshold, gene_id, "")),
                  size = 3, max.overlaps = Inf) +
  theme_minimal() +
  theme(axis.line = element_line(color = "gray10"),
        axis.title.x = element_text(margin = margin(t = 7)),    # Increase space between x-axis label and axis
        axis.title.y = element_text(margin = margin(r = 10))) + # Increase space between y-axis label and axis
  annotate("text", x = 1.5, y = 19.8, label = "FC = 1.5", size = 3, color = "black", vjust = -1) +                 # Label for right vertical line
  annotate("text", x = -10.5, y = -log10(0.15), label = "q-value = 0.05", size = 3, color = "black", hjust = -0.1)  # Label for horizontal line


# Zebrafish GO Enrichment Analysis ---------------------------------------------
zf_gseaData <- as.data.frame(de_ZFgenes)

# Get DE genes as vector
zf_deGeneNames <- zf_gseaData$gene_id
zf_backgroundGenes <- mycBrain_wildtypeDF$gene_id

zfGO_results <- enrichGO(gene = zf_deGeneNames, universe = zf_backgroundGenes, OrgDb = org.Dr.eg.db, 
                         keyType = "ENSEMBL", ont = "BP", 
                         pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.10)

zfGO_df <- as.data.frame(zfGO_results)

# Preparing the data for the GOplot
zfGO_bubble <- zfGO_df %>%
  mutate(Category = "BP") %>%
  # Renaming columns to match the desired format
  rename(Term = Description, adj_pval = p.adjust, genes = geneID) %>%
  dplyr::select(Category, ID, Term, genes, adj_pval) %>%
  mutate(genes = gsub("/", ", ", genes))

zf_geneList <- zf_gseaData %>%
  rename(ID = gene_id, logFC = log2FC, AveExpr = wildtype_brain, 
         P.value = Pvalue, adj.P.Val = FDR) %>%
  dplyr::select(ID, logFC, AveExpr, P.value, adj.P.Val)

zfGO_bubble <- circle_dat(zfGO_bubble, zf_geneList)

# Bubble plot
GOBubble(zfGO_bubble, labels = 1)

# Reduced plot with redundant terms (i.e. gene overlap > 0.75)
zfGO_reduced <- reduce_overlap(zfGO_bubble, overlap = 0.75)

GOBubble(zfGO_reduced, labels = 1) 

# ggplot2
ggplot(zfGO_reduced, aes(x = zscore, y = -log10(adj_pval))) +
  geom_point(aes(size = count, color = category), shape = 21, fill = "lightgreen") +
  geom_text(aes(label = ID), hjust = 1.1, vjust = 1.1, size = 3) +
  # scale_size_area(max_size = 10) +
  scale_size_continuous(range = c(2, 15)) + 
  scale_color_manual(values = c("BP" = "darkgreen")) +
  labs(x = "z-score", y = "-log(adj p-value)", size = "Count") +
  theme_minimal() +
  theme(legend.position = "bottom") +
  coord_cartesian(
    xlim = c(-10, 10),
    ylim = c(0, max(-log10(zfGO_reduced$adj_pval)))
  ) +
  guides(color = guide_legend(override.aes = list(size = 6)))


# Zebrafish Hypergeometric test (Gene Ontology Enrichment) --------------------
# Annotate genes with GO terms using biomaRt
# ONLY RUN ONCE
# mart <- biomaRt::useMart("ensembl",
#                          dataset = "drerio_gene_ensembl",
#                          host = "https://useast.ensembl.org")
# 
# biomart_attributes <- listAttributes(mart)
# 
# go_annotations <- biomaRt::getBM(attributes = c('ensembl_gene_id', 'go_id'),
#                                  filters = 'ensembl_gene_id',
#                                  values = gene_ids,
#                                  mart = mart)

# Remove genes without GO annotations
no_annotationGO      <- go_annotations[go_annotations$go_id == "", ]
go_annotationsFilter <- go_annotations[go_annotations$go_id != "", ]

# 6,964 for MYC Brain vs Wildtype Brain
# 6,743 for MYC Brain vs Control Brain
cat("GO annotations after filtering:", nrow(go_annotationsFilter), "\n")

# Perform hypergeometric test for GO term enrichment using clusterProfiler
gene_list <- unique(go_annotationsFilter$ensembl_gene_id)

# Background genes (double check depending on data)
background_genes <- unique(mycBrain_wildtypeDF$gene_id) # 25,979 genes
background_genes <- unique(mycBrain_controlDF$gene_id)  # 25,808 genes


# Conduct the enrichment analysis
enrichedGO <- enrichGO(gene = gene_list,
                       OrgDb = org.Dr.eg.db,
                       keyType = "ENSEMBL",
                       ont = "ALL",
                       pAdjustMethod = "BH",
                       universe = background_genes,
                       pvalueCutoff = 0.05, # threshold for significant GO term
                       qvalueCutoff = 0.10, # Benjamini-Hochberg p-value adjustment 
                                            # (to control FDR)
                       readable = TRUE)

enrichedGO_df <- as.data.frame(enrichedGO)

# 14 for MYC Brain vs Wildtype Brain
# 23 for MYC Brain vs Control Brain
cat("Number of enriched GO terms:", nrow(enrichedGO_df), "\n")

# Save results
write.csv(enrichedGO_df, file = "mycBrain_wildtypeGO_enrichment.csv")
write.csv(enrichedGO_df, file = "mycBrain_controlGO_enrichment.csv")

# Read GO Enrichment Results
mycBrain_wildtypeGO <- read_csv("mycBrain_wildtypeGO_enrichment.csv")
mycBrain_controlGO  <- read_csv("mycBrain_controlGO_enrichment.csv")

wildtype_goID <- mycBrain_wildtypeGO$ID
control_goID  <- mycBrain_controlGO$ID

diffGO_id <- setdiff(wildtype_goID, control_goID) # Returns empty: no GO IDs not included in control
diffGO_id <- setdiff(control_goID, wildtype_goID) # Returns 9 GO IDs

# GO IDs only in MYC Brain vs Control Brain
controlOnly_GO <- mycBrain_controlGO %>% 
  filter(ID %in% diffGO_id)

# Save results
write.csv(controlOnly_GO, file = "mycBrain_controlGO_only.csv")

# Functional Enrichment Analysis ----------------------------------------------
# 5,470/25,979 for wildtype
de_ZFgenes <- mycBrain_wildtypeDF %>%
  filter(abs(log2FC) >= log2(fc_threshold) & FDR <= fdr_threshold)

# 5,160/25,808 for control
de_ZFgenes <- mycBrain_controlDF %>%
  filter(abs(log2FC) >= log2(fc_threshold) & FDR <= fdr_threshold)

# Background genes (double check depending on data)
background_genes <- unique(mycBrain_wildtypeDF$gene_id) # 25,979 genes
background_genes <- unique(mycBrain_controlDF$gene_id)  # 25,808 genes

# Get gene IDs of the differentially expressed genes
gene_ids <- de_ZFgenes$gene_id

# Functional enrichment analysis using g:Profiler
diffExFunctional <- gost(query = gene_list,
                         organism = "drerio",
                         custom_bg = background_genes,
                         significant = TRUE,            # Only return significant results
                         user_threshold = 0.05,         # Adjusted p-value threshold
                         correction_method = "fdr")     # Benjamini-Hochberg correction

diffExFunctional_df <- as.data.frame(diffExFunctional$result)

# Convert "parents" column to string (instead of list) to save table
diffExFunctional_df$parents <- sapply(diffExFunctional_df$parents, 
                                      function(y) paste(y, collapse = ", "))
 
# Save results
write.csv(diffExFunctional_df, file = "mycBrain_wtFunctionalEnrichment.csv")
write.csv(diffExFunctional_df, file = "mycBrain_ctrlFunctionalEnrichment.csv")

# Gene Expression Response (MYC)  ----------------------------------------------

# Define gene symbols for zebrafish
zebrafish_symbols <- c("ENSDARG00000045695" = "myca", 
                       "ENSDARG00000007241" = "mycb", 
                       "ENSDARG00000077473" = "mych")

# Extract and average the expression values (HEALTHY)
# Human MYC in healthy tissue
human_MYChealthy <- pineoblastoma_tumoursDF %>%
  filter(human_gene == "MYC") %>%
  dplyr::select(starts_with("NFB")) %>%
  unlist() %>%
  as.numeric()

human_MYChealthy_mean <- mean(human_MYChealthy)
human_MYChealthy_sd <- sd(human_MYChealthy)
human_MYChealthy_CI <- qt(0.975, df=length(human_MYChealthy)-1) * human_MYChealthy_sd / sqrt(length(human_MYChealthy))

# Extract the expression values for zebrafish myca, mycb, mych in wildtype + control tissue
wildtype_data <- mycBrain_wildtypeDF %>%
  filter(gene_id %in% c("ENSDARG00000045695", "ENSDARG00000007241", "ENSDARG00000077473"))

control_data <- mycBrain_controlDF %>%
  filter(gene_id %in% c("ENSDARG00000045695", "ENSDARG00000007241", "ENSDARG00000077473"))

zebrafish_wildtype <- mycBrain_wildtypeDF %>%
  filter(gene_id %in% names(zebrafish_symbols)) %>%
  mutate(gene_id = zebrafish_symbols[gene_id]) %>%
  pivot_longer(cols = starts_with("WT"),
               names_to = "sample", values_to = "expression") %>%
  group_by(gene_id) %>%
  summarise(mean_expression = mean(expression),
            sd_expression = sd(expression),
            ci_expression = qt(0.975, df=n()-1) * sd_expression / sqrt(n()))

zebrafish_control <- mycBrain_controlDF %>%
  filter(gene_id %in% names(zebrafish_symbols)) %>%
  mutate(gene_id = zebrafish_symbols[gene_id]) %>%
  pivot_longer(cols = starts_with("ctr"),
               names_to = "sample", values_to = "expression") %>%
  group_by(gene_id) %>%
  summarise(mean_expression = mean(expression),
            sd_expression = sd(expression),
            ci_expression = qt(0.975, df=n()-1) * sd_expression / sqrt(n()))

# Extract and average the expression values (TUMOURS)
# Human MYC in tumor tissue
human_MYCtumour <- pineoblastoma_tumoursDF %>%
  filter(human_gene == "MYC") %>%
  dplyr::select(starts_with("RBTC")) %>%
  unlist() %>%
  as.numeric()

human_MYCtumour_mean <- mean(human_MYCtumour)
human_MYCtumour_sd <- sd(human_MYCtumour)
human_MYCtumour_CI <- qt(0.975, df=length(human_MYCtumour)-1) * human_MYCtumour_sd / sqrt(length(human_MYCtumour))

# Zebrafish myca, mycb, mych in tumor tissue
zebrafish_tumour <- mycBrain_wildtypeDF %>%
  filter(gene_id %in% names(zebrafish_symbols)) %>%
  mutate(gene_id = zebrafish_symbols[gene_id]) %>%
  dplyr::select(gene_id, starts_with("MYC"), starts_with("myc")) %>%
  pivot_longer(cols = starts_with("MYC") | starts_with("myc"),
               names_to = "sample", values_to = "expression") %>%
  group_by(gene_id) %>%
  summarise(mean_expression = mean(expression),
            sd_expression = sd(expression),
            ci_expression = qt(0.975, df=n()-1) * sd_expression / sqrt(n()))

# Combine data for plotting
human_healthy <- data.frame(gene_id = "MYC", sample = "Healthy",
  mean_expression = human_MYChealthy_mean,
  sd_expression = human_MYChealthy_sd,
  ci_expression = human_MYChealthy_CI)

human_tumour <- data.frame(gene_id = "MYC", sample = "Tumour",
  mean_expression = human_MYCtumour_mean,
  sd_expression = human_MYCtumour_sd,
  ci_expression = human_MYCtumour_CI)

zebrafish_wildtype <- zebrafish_wildtype %>%
  mutate(sample = "Wildtype")

zebrafish_control <- zebrafish_control %>%
  mutate(sample = "Control")

zebrafish_tumour <- zebrafish_tumour %>%
  mutate(sample = "Tumour")

# PCA of Human + Zebrafish Genes (Marker Genes Only) ---------------------------

# Human TPM Data
humanTPM_data <- as.data.frame(pb_myc_fetalbrain_tpm)

# Zebrafish TPM Data
zebrafishTPM_data <- mycBrain_wildtypeDF[, c("MYC3_Brain", "MYC4_Brain", "MYC5_Brain", "myc1_brain", "myc2_brain",
                                             "WT_Brain2", "WT_6_B", "WT_7_B", "WT_8_B")]
zebrafishTPM_data$gene_id <- mycBrain_wildtypeDF$gene_id
zebrafishTPM_data <- zebrafishTPM_data %>%
  relocate(gene_id, .before = MYC3_Brain)

# Define marker genes
marker_genes_human <- c("CRX", "FOXR2", "MYC", "RB1",
                        "LIN28A", "SYP", "CHGA", "GFAP", "OLIG2", "MKI67")

# Gene mapping for zebrafish orthologs to human genes
gene_mapping <- data.frame(
  gene_symbol = c("crx", "mki67", "sypa", "sypb", "chga", "gfap",
                  "foxr1", "rb1", "myca", "mycb", "mych"),
  ensembl_id = c("ENSDARG00000011989", "ENSDARG00000091150", "ENSDARG00000110528", 
                 "ENSDARG00000002230", "ENSDARG00000008829", "ENSDARG00000025301",
                 "ENSDARG00000004864", "ENSDARG00000006782", "ENSDARG00000045695", 
                 "ENSDARG00000007241", "ENSDARG00000077473")
)

# Add corresponding human gene names to the gene mapping
gene_mapping$human_gene <- recode(gene_mapping$gene_symbol,
                                  "crx" = "CRX",
                                  "mki67" = "MKI67",
                                  "sypa" = "SYP",
                                  "sypb" = "SYP",
                                  "chga" = "CHGA",
                                  "gfap" = "GFAP",
                                  "foxr1" = "FOXR2",
                                  "rb1" = "RB1",
                                  "myca" = "MYC",
                                  "mycb" = "MYC",
                                  "mych" = "MYC")

# Filter human TPM data to include only marker genes
human_tpm_filtered <- humanTPM_data %>%
  filter(human_gene %in% marker_genes_human) %>%
  dplyr::select(human_gene, starts_with("RBTC"))

# Filter zebrafish TPM data to include only marker genes
zebrafish_tpm_filtered <- zebrafishTPM_data %>%
  filter(gene_id %in% gene_mapping$ensembl_id) %>%
  dplyr::select(gene_id, starts_with("MYC", ignore.case = TRUE))

# Add corresponding human gene names to zebrafish TPM data
zebrafish_tpm_filtered <- zebrafish_tpm_filtered %>%
  left_join(gene_mapping, by = c("gene_id" = "ensembl_id"))

# Average zebrafish orthologs for each human gene
zebrafish_tpm_avg <- zebrafish_tpm_filtered %>%
  group_by(human_gene) %>%
  summarize(across(starts_with("MYC"), mean, na.rm = TRUE))

# Convert back to matrix form and set row names
zebrafish_tpm_avg <- as.data.frame(zebrafish_tpm_avg)
rownames(zebrafish_tpm_avg) <- zebrafish_tpm_avg$human_gene
zebrafish_tpm_avg <- zebrafish_tpm_avg[, -1]

# Align TPM Data Based on Marker Genes
# Subset human TPM data to match zebrafish data genes
human_tpm_aligned <- human_tpm_filtered[match(rownames(zebrafish_tpm_avg), human_tpm_filtered$human_gene), , drop = FALSE]
human_tpm_aligned <- human_tpm_aligned %>% dplyr::select(-human_gene)

# Transpose and Combine Data
human_tpm_transposed <- t(human_tpm_aligned)
zebrafish_tpm_transposed <- t(zebrafish_tpm_avg)

# Combine data
combined_tpm_de <- rbind(human_tpm_transposed, zebrafish_tpm_transposed)

# Remove columns with zero variance
combined_tpm_de <- combined_tpm_de[, apply(combined_tpm_de, 2, var) != 0]

# Perform PCA
pca_result_de <- prcomp(combined_tpm_de, scale. = TRUE)

# Create a data frame with PCA results
num_human_samples <- nrow(human_tpm_transposed)
num_zebrafish_samples <- nrow(zebrafish_tpm_transposed)

pca_data_de <- data.frame(
  Sample = rownames(pca_result_de$x),
  PC1 = pca_result_de$x[, 1],
  PC2 = pca_result_de$x[, 2],
  Species = rep(c("Human", "Zebrafish"), c(num_human_samples, num_zebrafish_samples))
)

# Plot PCA
ggplot(pca_data_de, aes(x = PC1, y = PC2, color = Species)) +
  geom_point(size = 3) +
  scale_color_manual(values = c("Human" = "#E3F0C6", "Zebrafish" = "#89D6F4")) +
  theme_minimal() +
  theme(
    axis.title.x = element_blank(),  # Remove x axis title
    axis.title.y = element_blank(),  # Remove y axis title
    panel.border = element_rect(color = "gray20", fill = NA, size = 1)  # Dark gray border around the plot
  ) +
  labs(color = "Species")  # Label for the color legend

# With ellipses
border_colors <- c("Human" = "#C0D09C", "Zebrafish" = "#5CB9E1")

# Plot
ggplot(pca_data_de, aes(x = PC1, y = PC2, color = Species)) +
  geom_point(size = 4, shape = 21, fill = pca_data_de$Species %>% recode("Human" = "#E3F0C6", "Zebrafish" = "#89D6F4"), 
             color = pca_data_de$Species %>% recode("Human" = "#C0D09C", "Zebrafish" = "#5CB9E1"), alpha = 0.6) + 
  scale_color_manual(values = c("Human" = "#E3F0C6", "Zebrafish" = "#89D6F4")) +  # Custom colors for points
  stat_ellipse(level = 0.95, geom = "polygon", aes(fill = Species), alpha = 0.1, size = 1, color = NA) +  # Add ellipses with transparency
  theme_minimal() +
  theme(
    axis.title.x = element_blank(),  # Remove x axis title
    axis.title.y = element_blank(),  # Remove y axis title
    panel.border = element_rect(color = "gray20", fill = NA, size = 1),  # Dark gray border around the plot
    axis.ticks = element_line(color = "black"),  # Add tick marks
    axis.text = element_text(color = "black")  # Set tick mark text color
  ) +
  labs(color = "Species", fill = "Species") + # Label for the color legend 
  guides(color = guide_legend(override.aes = list(shape = 21, fill = fill_colors, alpha = 0.6)))


# Define the colors
fill_colors <- c("Human" = "#E3F0C6", "Zebrafish" = "#89D6F4")
border_colors <- c("Human" = "#C0D09C", "Zebrafish" = "#5CB9E1")

# Plot PCA with custom colors, borders, and styling
plot<-ggplot(pca_data_de, aes(x = PC1, y = PC2, color = Species, fill = Species)) +
  geom_point(size = 3, shape = 16, stroke = 1, alpha = 0.6, aes(fill = Species)) +  # Adjust circle size, border color, and fill alpha
  scale_color_manual(values = border_colors) +  # Custom colors for borders
  scale_fill_manual(values = fill_colors) +  # Fill colors for points and ellipses
  stat_ellipse(level = 0.95, geom = "polygon", aes(fill = Species, color = Species), alpha = 0.1, size = 0) +  # Add ellipses with fill and border color
  theme_minimal() +
  theme(
    axis.title.x = element_blank(),  # Remove x axis title
    axis.title.y = element_blank(),  # Remove y axis title
    panel.border = element_rect(color = "gray20", fill = NA, size = 1),  # Dark gray border around the plot
    panel.background = element_rect(fill = "white", color = NA),
    axis.ticks = element_line(color = "black"),  # Add tick marks
    axis.text = element_text(color = "black"),  # Set tick mark text color
    legend.position = "none"  # Remove the legend
  ) +
  labs(color = "Species", fill = "Species") +  # Label for the color legend
  guides(color = guide_legend(override.aes = list(shape = 21, fill = fill_colors, alpha = 0.6)))

# Save the plot with the specified size in cm
ggsave("PCA.png", plot = plot, device = "png", 
       width = 23.01, height = 14, units = "cm", dpi = 900)

# PCA of Human + Zebrafish Genes (Orthologs Only) ------------------------------

# Human TPM Data
humanTPM_data <- as.data.frame(pb_myc_fetalbrain_tpm)

# Zebrafish TPM Data
# Extract only the count data columns
zebrafishTPM_data <- mycBrain_wildtypeDF[, c("MYC3_Brain", "MYC4_Brain", "MYC5_Brain", "myc1_brain", "myc2_brain",
                                             "WT_Brain2", "WT_6_B", "WT_7_B", "WT_8_B")]

zebrafishTPM_data$gene_id <- mycBrain_wildtypeDF$gene_id

zebrafishTPM_data <- zebrafishTPM_data %>%
  relocate(gene_id, .before = MYC3_Brain)

# Step 1: Filter DE Genes based on the specified thresholds
de_genes_human <- pineoblastoma_tumoursDF %>%
  filter(abs(log2FoldChange) >= 1.5 & padj < 0.1) %>%
  pull(gene)

de_genes_zebrafish <- mycBrain_wildtypeDF %>%
  filter(abs(log2FC) >= 1.5 & FDR < 0.1) %>%
  pull(gene_id)

# Step 2: Identify DE Orthologous Genes
# Filter the orthologs based on the DE genes
human_orthologs_de <- HGNChuman_zebrafish %>%
  separate_rows(consensus_ensembl_gene, sep = ",") %>%
  mutate(consensus_ensembl_gene = str_trim(consensus_ensembl_gene)) %>%  # Trim whitespace
  filter(human_symbol %in% de_genes_human & !is.na(consensus_ensembl_gene) & consensus_ensembl_gene != "NA") %>%
  dplyr::select(human_symbol, consensus_ensembl_gene)

zebrafish_orthologs_de <- HGNCzebrafish_human %>%
  separate_rows(consensus_ensembl_gene, sep = ",") %>%
  mutate(consensus_ensembl_gene = str_trim(consensus_ensembl_gene)) %>%  # Trim whitespace
  filter(consensus_ensembl_gene %in% de_genes_zebrafish & !is.na(consensus_ensembl_gene) & consensus_ensembl_gene != "NA") %>%
  dplyr::select(human_symbol, consensus_ensembl_gene)

# Merge to get a unified list of DE orthologs
common_orthologs_de <- merge(human_orthologs_de, zebrafish_orthologs_de, by = "consensus_ensembl_gene")

# Step 3: Filter and Align TPM Data Based on DE Orthologous Genes

# Filter human TPM data
human_tpm_filtered <- humanTPM_data %>%
  filter(human_gene %in% common_orthologs_de$human_symbol.y) %>%
  dplyr::select(human_gene, starts_with("RBTC"))

# Filter zebrafish TPM data
zebrafish_tpm_filtered <- zebrafishTPM_data %>%
  filter(gene_id %in% common_orthologs_de$consensus_ensembl_gene) %>%
  dplyr::select(gene_id, starts_with("MYC", ignore.case = TRUE))

# Ensure the rownames are the gene IDs
rownames(human_tpm_filtered) <- human_tpm_filtered$human_gene
rownames(zebrafish_tpm_filtered) <- zebrafish_tpm_filtered$gene_id

# Remove the gene columns
human_tpm_filtered <- human_tpm_filtered[, -1]
zebrafish_tpm_filtered <- zebrafish_tpm_filtered[, -1]

# Step 4: Align TPM Data Based on Sorted DE Orthologs
# Sort the rows in the same order for both human and zebrafish data frames
human_tpm_aligned <- human_tpm_filtered[match(common_orthologs_de$human_symbol.y, rownames(human_tpm_filtered)), ]
zebrafish_tpm_aligned <- zebrafish_tpm_filtered[match(common_orthologs_de$consensus_ensembl_gene, rownames(zebrafish_tpm_filtered)), ]

# Verify alignment and drop misaligned genes
aligned_indices <- which(rownames(human_tpm_aligned) == common_orthologs_de$human_symbol.x & rownames(zebrafish_tpm_aligned) == common_orthologs_de$consensus_ensembl_gene)
human_tpm_aligned <- human_tpm_aligned[aligned_indices, ]
zebrafish_tpm_aligned <- zebrafish_tpm_aligned[aligned_indices, ]
common_orthologs_de <- common_orthologs_de[aligned_indices, ]

# Step 5: Transpose and Combine Data
# Transpose data to have genes as columns and samples as rows
human_tpm_transposed <- t(human_tpm_aligned)
zebrafish_tpm_transposed <- t(zebrafish_tpm_aligned)

# Combine data
combined_tpm_de <- rbind(human_tpm_transposed, zebrafish_tpm_transposed)

# Remove columns with zero variance
combined_tpm_de <- combined_tpm_de[, apply(combined_tpm_de, 2, var) != 0]

# Step 6: Perform PCA
pca_result_de <- prcomp(combined_tpm_de, scale. = TRUE)

pca_genes <- colnames(combined_tpm_de)

# Create a data frame with PCA results
num_human_samples <- nrow(human_tpm_de_transposed)
num_zebrafish_samples <- nrow(zebrafish_tpm_de_transposed)

pca_data_de <- data.frame(
  Sample = rownames(pca_result_de$x),
  PC1 = pca_result_de$x[, 1],
  PC2 = pca_result_de$x[, 2],
  Species = rep(c("Human", "Zebrafish"), c(num_human_samples, num_zebrafish_samples))
)

# Step 7: Plot PCA
ggplot(pca_data_de, aes(x = PC1, y = PC2, color = Species)) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(x = "Principal Component 1",
       y = "Principal Component 2")

# Heatmap Clustering Orthologs -------------------------------------------------

# Step 1: Filter DE Genes
de_genes_human <- pineoblastoma_tumoursDF %>%
  filter(abs(log2FoldChange) >= 2 & padj < 0.05) %>%
  pull(gene)

de_genes_zebrafish <- mycBrain_wildtypeDF %>%
  filter(abs(log2FC) >= 2 & FDR < 0.05) %>%
  pull(gene_id)

# Step 2: Identify DE Orthologous Genes
orthologs_de <- HGNChuman_zebrafish %>%
  separate_rows(consensus_ensembl_gene, sep = ",") %>%
  mutate(consensus_ensembl_gene = str_trim(consensus_ensembl_gene)) %>%
  filter(human_symbol %in% de_genes_human & 
           !is.na(consensus_ensembl_gene) & 
           consensus_ensembl_gene != "NA") %>%
  dplyr::select(human_symbol, consensus_ensembl_gene)

# Step 3: Filter and Align TPM Data Based on DE Orthologous Genes
# Filter human TPM data
human_tpm_de_filtered <- humanTPM_data %>%
  filter(human_gene %in% orthologs_de$human_symbol) %>%
  dplyr::select(human_gene, starts_with("RBTC"))

# Filter zebrafish TPM data
zebrafish_tpm_de_filtered <- zebrafishTPM_data %>%
  filter(gene_id %in% orthologs_de$consensus_ensembl_gene) %>%
  dplyr::select(gene_id, starts_with("MYC", ignore.case = TRUE))

# Ensure row names are gene identifiers
rownames(human_tpm_de_filtered) <- human_tpm_de_filtered$human_gene
human_tpm_de_filtered <- human_tpm_de_filtered %>%
  dplyr::select(-human_gene)

rownames(zebrafish_tpm_de_filtered) <- zebrafish_tpm_de_filtered$gene_id
zebrafish_tpm_de_filtered <- zebrafish_tpm_de_filtered %>%
  dplyr::select(-gene_id)

# Convert dataframes to matrices
human_tpm_matrix <- as.matrix(human_tpm_de_filtered)
zebrafish_tpm_matrix <- as.matrix(zebrafish_tpm_de_filtered)

# Ensure all data in the matrices are numeric
human_tpm_matrix <- apply(human_tpm_matrix, 2, as.numeric)
zebrafish_tpm_matrix <- apply(zebrafish_tpm_matrix, 2, as.numeric)

# Set row names for the matrices
rownames(human_tpm_matrix) <- rownames(human_tpm_de_filtered)
rownames(zebrafish_tpm_matrix) <- rownames(zebrafish_tpm_de_filtered)

# Function to remove rows with NA/NaN/Inf values
remove_na_inf <- function(matrix) {
  matrix <- matrix[complete.cases(matrix) & apply(matrix, 1, function(row) all(is.finite(row))), ]
  return(matrix)
}

# Remove rows with NA/NaN/Inf values
human_tpm_matrix <- remove_na_inf(human_tpm_matrix)
zebrafish_tpm_matrix <- remove_na_inf(zebrafish_tpm_matrix)

# Align matrices by ensuring rownames are mapped correctly using orthologs_de
# Filter orthologs_de to keep only those orthologs present in both matrices
orthologs_de <- orthologs_de %>%
  filter(human_symbol %in% rownames(human_tpm_matrix) &
           consensus_ensembl_gene %in% rownames(zebrafish_tpm_matrix))

# Align matrices based on orthologs
aligned_human_matrix <- human_tpm_matrix[orthologs_de$human_symbol, ]
aligned_zebrafish_matrix <- zebrafish_tpm_matrix[orthologs_de$consensus_ensembl_gene, ]

# Combine the matrices side-by-side
combined_matrix <- cbind(aligned_human_matrix, aligned_zebrafish_matrix)

# Verify if there are any NA/NaN/Inf values in the combined matrix
print(any(!is.finite(combined_matrix)))  # Diagnostic check

if (any(!is.finite(combined_matrix))) {
  stop("There are still NA/NaN/Inf values in the combined matrix.")
}

# Check
print(summary(combined_matrix))
is.infinite(combined_matrix) %>% table()
is.na(combined_matrix) %>% table()

# Average duplicate rows by gene name
combined_matrix <- aggregate(combined_matrix, by = list(rownames(combined_matrix)), FUN = mean)
rownames(combined_matrix) <- combined_matrix$Group.1
combined_matrix <- combined_matrix[, -1]

# Remove rows with all 0 values or same values across all samples
# Error in hclust(d, method = method) : 
# NA/NaN/Inf in foreign function call (arg 10)
# https://www.biostars.org/p/446761/
remove_zero_and_same_rows <- function(matrix) {
  matrix <- matrix[apply(matrix, 1, function(row) !all(row == 0)), ]
  matrix <- matrix[apply(matrix, 1, function(row) length(unique(row)) > 1), ]
  return(matrix)
}

combined_matrix <- remove_zero_and_same_rows(combined_matrix)

# Generate heatmap with hierarchical clustering
pheatmap(combined_matrix, cluster_rows = TRUE, cluster_cols = TRUE, 
         scale = "row", show_rownames = TRUE, show_colnames = TRUE)

# Pathway Enrichment -----------------------------------------------------------

# FDR < 0.05
# humanDE_genes <- pineoblastoma_tumoursDF %>% filter(padj < 0.05)
humanDE_genes <- pineoblastoma_tumoursDF %>% filter(padj < 0.05 & abs(log2FoldChange) > 1.5)
gene_list <- humanDE_genes$log2FoldChange
names(gene_list) <- humanDE_genes$gene

# Convert gene IDs from HGNC symbols to ENTREZID
# 13.73% of genes fail to map
gene_conversion <- bitr(names(gene_list), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# Filter the gene list to include only those that were successfully converted
gene_list <- gene_list[gene_conversion$SYMBOL]
names(gene_list) <- gene_conversion$ENTREZID

# Sort gene list in decreasing order
gene_list <- sort(gene_list, decreasing = TRUE)

# Define the organism code for KEGG
kegg_organism <- "hsa"

# Perform KEGG pathway enrichment analysis
human_pathways <- gseKEGG(geneList = gene_list,
                          organism = "hsa",
                          nPerm = 10000,
                          minGSSize = 3,
                          maxGSSize = 800,
                          pvalueCutoff = 0.05,
                          pAdjustMethod = "none",
                          keyType = "ncbi-geneid")

# # Access the results dataframe
# human_pathways_results <- human_pathways@result
# 
# # Filter pathways to include only those with 40+ genes in the count
# filtered_humanPathways <- human_pathways_results %>% filter(setSize >= 40)

# Create the dot plot with facets for activated and suppressed pathways
plot<-dotplot(human_pathways, showCategory = 10, split = ".sign") + 
  facet_grid(. ~ .sign, labeller = as_labeller(c("activated" = "Activated", "suppressed" = "Suppressed"))) +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    strip.text = element_text(size = 14),
    axis.text.y = element_text(size = 12)
  ) +
  scale_color_gradient(low = "#B0DE7E", high = "#0070C0", name = "FDR") 

# Save the plot with the specified size in cm
ggsave("Human Pathways.png", plot = plot, device = "png", 
       width = 23.01, height = 23, units = "cm", dpi = 900)



# FDR < 0.05
zebrafishDE_genes <- mycBrain_wildtypeDF %>% filter(FDR < 0.05 & abs(log2FC) > 1.5)

gene_list <- zebrafishDE_genes$log2FC
names(gene_list) <- zebrafishDE_genes$gene_id

# Remove infinite values from gene_list
gene_list <- gene_list[is.finite(gene_list)]

# Convert gene IDs from ENSEMBL to SYMBOL
# 14.48% fail to map
# gene_conversion <- bitr(names(gene_list), fromType = "ENSEMBL", toType = "SYMBOL", OrgDb = org.Dr.eg.db)
gene_conversion <- bitr(names(gene_list), fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Dr.eg.db)

# Filter the gene list to include only those that were successfully converted
mapped_genes <- gene_conversion %>% filter(!is.na(ENTREZID))
gene_list <- gene_list[mapped_genes$ENSEMBL]
names(gene_list) <- mapped_genes$ENTREZID

# Create a mapping from zebrafish gene IDs to human orthologs
gene_mapping <- HGNCzebrafish_human %>%
  select(zebrafish_gene_id, human_entrez_gene) %>%
  distinct()

gene_mapping <- gene_mapping %>%
  group_by(zebrafish_gene_id) %>%
  summarise(human_entrez_gene = first(human_entrez_gene))

# Merge the gene list with the mapping to replace names
mapped_genes <- merge(data.frame(zebrafish_gene_id = names(gene_list), value = gene_list, stringsAsFactors = FALSE),
                      gene_mapping, by.x = "zebrafish_gene_id", by.y = "zebrafish_gene_id", all.x = TRUE)

# Replace the names of gene_list with the mapped human orthologs
mapped_gene_list <- mapped_genes$value
names(mapped_gene_list) <- mapped_genes$human_entrez_gene

# Ensure the gene list is sorted in decreasing order
mapped_gene_list <- sort(mapped_gene_list, decreasing = TRUE)
gene_list <- sort(gene_list, decreasing = TRUE)

# # Add the "dre:" prefix to each gene ID in the gene_list
# names(gene_list) <- paste0("dre:", names(gene_list))


# Define the organism code for KEGG
kegg_organism <- "hsa"

# Perform KEGG pathway enrichment analysis
zebrafish_pathways <- gseKEGG(geneList = gene_list,
                              organism = "dre",
                              nPerm = 10000,
                              minGSSize = 3,
                              maxGSSize = 800,
                              pvalueCutoff = 0.05,
                              pAdjustMethod = "none",
                              keyType = "ncbi-geneid")

# Create the dot plot with facets for activated and suppressed pathways
plot<-dotplot(zebrafish_pathways, showCategory = 10, split = ".sign") + 
  facet_grid(. ~ .sign, labeller = as_labeller(c("activated" = "Activated", "suppressed" = "Suppressed"))) +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    strip.text = element_text(size = 14),
    axis.text.y = element_text(size = 12)
  ) +
  scale_color_gradient(low = "#B0DE7E", high = "#0070C0", name = "FDR")  # Customize the gradient

# Save the plot with the specified size in cm
ggsave("Zebrafish Pathways.png", plot = plot, device = "png", 
       width = 23.01, height = 23, units = "cm", dpi = 900)

# GO enrichment analysis ------------------------------------------------------
# human_biologicalProcesses <- enrichGO(gene = names(gene_list),
#                                       OrgDb = org.Hs.eg.db,
#                                       keyType = "ENTREZID",
#                                       ont = "BP",  # Biological Process
#                                       pAdjustMethod = "BH",
#                                       pvalueCutoff = 0.05,
#                                       qvalueCutoff = 0.2)

human_biologicalProcesses <- gseGO(geneList=gene_list, 
                                   ont ="BP", 
                                   keyType = "ENTREZID", 
                                   nPerm = 10000, 
                                   minGSSize = 3, 
                                   maxGSSize = 800, 
                                   pvalueCutoff = 0.05, 
                                   verbose = TRUE,
                                   OrgDb = org.Hs.eg.db, 
                                   pAdjustMethod = "none")


# Split results by sign of log2FoldChange
ego@result$sign <- ifelse(gene_list[ego@result$ID] > 0, "activated", "suppressed")

plot<-dotplot(human_biologicalProcesses, showCategory = 10, split = ".sign") + 
  facet_grid(. ~ .sign, labeller = as_labeller(c("activated" = "Activated", "suppressed" = "Suppressed"))) +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    strip.text = element_text(size = 12),
    axis.text.y = element_text(size = 8)
  ) +
  scale_color_gradient(low = "#B0DE7E", high = "#0070C0", name = "FDR")  # Customize the gradient

ggsave("Human Enriched GO-BP.png", plot = plot, device = "png", 
       width = 17.54, height = 21.38, units = "cm", dpi = 900)



zebrafish_biologicalProcesses <- gseGO(geneList=gene_list, 
                                   ont ="BP", 
                                   keyType = "ENTREZID", 
                                   nPerm = 10000, 
                                   minGSSize = 3, 
                                   maxGSSize = 800, 
                                   pvalueCutoff = 0.05, 
                                   verbose = TRUE,
                                   OrgDb = org.Dr.eg.db, 
                                   pAdjustMethod = "none")

# Extract descriptions and wrap them
descriptions <- slot(zebrafish_biologicalProcesses, "Description")
wrapped_descriptions <- str_wrap(descriptions, width = 40)
slot(zebrafish_biologicalProcesses, "Description") <- wrapped_descriptions

plot<-dotplot(zebrafish_biologicalProcesses, showCategory = 10, split = ".sign") + 
  facet_grid(. ~ .sign, labeller = as_labeller(c("activated" = "Activated", "suppressed" = "Suppressed"))) +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    strip.text = element_text(size = 12),
    axis.text.y = element_text(size = 8),
    #plot.margin = margin(t = 5, r = 5, b = 5, l = 100)  # Increase the left margin of the plot
  ) +
  scale_color_gradient(low = "#B0DE7E", high = "#0070C0", name = "FDR")  # Customize the gradient

ggsave("Zebrafish Enriched GO-BP.png", plot = plot, device = "png", 
       width = 17.54, height = 21.38, units = "cm", dpi = 900)


# Venn Diagram -----------------------------------------------------------------
# Human processing
humanDE_genes$gene <- trimws(humanDE_genes$gene)
HGNChuman_zebrafish$human_symbol <- trimws(HGNChuman_zebrafish$human_symbol)

# Zebrafish processing
zebrafishDE_genes$gene_id <- trimws(zebrafishDE_genes$gene_id)

split_zebrafishGenes <- strsplit(HGNChuman_zebrafish$consensus_ensembl_gene, split = ",")

# Unlist the resulting list to get a vector of individual gene IDs
all_zebrafishGenes <- unlist(split_zebrafishGenes)
all_zebrafishGenes <- trimws(all_zebrafishGenes)

sum(humanDE_genes$gene %in% HGNChuman_zebrafish$human_symbol)
sum(zebrafishDE_genes$gene_id %in% all_zebrafishGenes)

# Define the sets of genes
human_orthologous_DE_genes <- humanDE_genes$gene[humanDE_genes$gene %in% HGNChuman_zebrafish$human_symbol]
human_exclusive_DE_genes <- humanDE_genes$gene[!humanDE_genes$gene %in% HGNChuman_zebrafish$human_symbol]

# Convert zebrafish genes to human genes based on orthologs
zebrafish_DE_orthologous_human_genes <- unique(HGNChuman_zebrafish$human_symbol[
  HGNChuman_zebrafish$consensus_ensembl_gene %in% all_zebrafishGenes
])

# Flatten the mapping data frame
flattened_mapping <- HGNChuman_zebrafish %>%
  mutate(consensus_ensembl_gene = str_split(consensus_ensembl_gene, ",")) %>%
  unnest(consensus_ensembl_gene)

zebrafish_DE_orthologous_human_genes <- unique(flattened_mapping$human_symbol[unique(zebrafishDE_genes$gene_id) %in% flattened_mapping$consensus_ensembl_gene])

n_distinct(flattened_mapping$human_symbol[zebrafishDE_genes$gene_id %in% flattened_mapping$consensus_ensembl_gene])


# Create list of sets
vennDiagram_data <- list(
  "Human DE Genes" = human_orthologous_DE_genes,
  "Zebrafish DE Genes" = zebrafish_DE_orthologous_human_genes
)

ggvenn(vennDiagram_data, 
  fill_color = c("#5A9FCD", "#A9D3C0", "#93C6C9", "#C7E2A9"),
  stroke_size = 0.5, set_name_size = 2
)

ggsave("Venn Diagram.png", plot = plot, device = "png", 
       width = 32.18, height = 15.25, units = "cm", dpi = 900)

# Find orthologous genes in human and zebrafish DE sets
human_orthologs_in_DE <- humanDE_genes$gene[humanDE_genes$gene %in% HGNChuman_zebrafish$human_symbol]
zebrafish_orthologs_in_DE <- zebrafishDE_genes$gene_id[zebrafishDE_genes$gene_id %in% all_zebrafishGenes]