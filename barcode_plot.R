# setwd("~/git.repos/auto_barcode/")
# file.name <- "log_barcodes_observed.fq_multi_fq.bar_RK11.barcode.tsv";
# expt.name <- "RK11"

barcode_plot <- function(file.name, expt.name) {
  library("ggplot2")
  log <- read.table(file = file.name, header = TRUE, sep = "\t")
  log.plot <-
    ggplot(
      data = log,
      aes(x = factor(matched), y = count / 1000000)) + 
    geom_boxplot(outlier.size = 0) +
    geom_jitter() +
    ggtitle(label = expt.name) + 
    xlab(label = "") +
    scale_y_continuous(name = "Count\n(million reads)")
  ggsave(
    filename = paste(sep = "", expt.name, ".barcodes.png"),
    plot     = log.plot,
    width    = 4,
    height   = 5
  )
}