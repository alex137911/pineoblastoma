suppressMessages(library(dbplyr))
suppressMessages(library(biomaRt))


# Fetch all attributes available in the Ensembl Biomart
# Want "drerio_homolog_ensembl_gene" and "drerio_homolog_associated_gene_name"
mart <- useEnsembl(biomart = "genes", 
                   dataset = "hsapiens_gene_ensembl",
                   version = 109)

biomart_attributes <- listAttributes(mart)

zebrafish_geneLookup <- getBM(mart = mart,
                              attributes = c('drerio_homolog_associated_gene_name', 
                                             'drerio_homolog_ensembl_gene'),
                              uniqueRows = TRUE)
