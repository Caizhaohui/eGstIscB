# Merge fragment tables, keeping the maximum absolute fold change per mutation.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) stop("Usage: Rscript 07_merge_fragment_foldchanges.R OUTPUT.csv FRAGMENT1.csv FRAGMENT2.csv [FRAGMENT...]")
suppressPackageStartupMessages({ library(dplyr); library(readr); library(stringr) })
files <- args[-1]
merged <- lapply(seq_along(files), function(i) {
  data <- read_csv(files[[i]], show_col_types = FALSE)
  data$Fragment <- paste0("fragment", i)
  data
}) %>% bind_rows() %>% mutate(global_row = row_number())
result <- merged %>% group_by(Mutation) %>%
  slice_max(order_by = abs(log2_foldchange), n = 1, with_ties = FALSE) %>%
  ungroup() %>% arrange(global_row) %>% select(-global_row)
write_csv(result, args[[1]])
