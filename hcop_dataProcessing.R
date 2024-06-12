suppressMessages(library(readr))
suppressMessages(library(dplyr))
suppressMessages(library(irr))


outpath <- "C:/Users/acale/OneDrive/Documents/Waterloo BME/Co-op/4C/pineoblastoma/Data/Output/"
outDir  <- sprintf("%s/Output", dirname(outpath))
setwd(outDir)

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
