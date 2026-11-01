# Compare two fragment-level Mutation/Count tables using the local analysis rule.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) stop("Usage: Rscript 06_compare_fragment_counts.R EARLIER.csv LATER.csv OUTPUT.csv")
suppressPackageStartupMessages(library(dplyr))
process_data <- function(data) {
  data %>% filter(Count >= 5) %>% mutate(Frequency = Count / sum(Count))
}
first <- process_data(read.csv(args[[1]]))
second <- process_data(read.csv(args[[2]]))
merged <- merge(first, second, by = "Mutation", suffixes = c("_file1", "_file2"))
result <- merged %>%
  mutate(log2_foldchange = log2((Frequency_file2 + 1e-6) / (Frequency_file1 + 1e-6))) %>%
  group_by(Mutation) %>% slice_max(order_by = abs(log2_foldchange), n = 1, with_ties = FALSE) %>% ungroup()
write.csv(result, args[[3]], row.names = FALSE)
