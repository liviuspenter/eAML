read_tcr_bcr <- function(file, file2, prefix, bcrtcr) {
  TCR <- read.csv(file)
  TCR$barcode <- paste0(prefix, "_", TCR$barcode)
  TCR <- TCR[!duplicated(TCR$barcode), ]
  TCR <- TCR[, c("barcode", "raw_clonotype_id")]
  TCR$raw_clonotype_id <- paste0(prefix, "_", TCR$raw_clonotype_id)
  names(TCR)[names(TCR) == "raw_clonotype_id"] <- "clonotype_id"
  Clonotypes <- read.csv(file2)
  Clonotypes$clonotype_id <- paste0(prefix, "_", Clonotypes$clonotype_id)
  TCR <- merge(TCR, Clonotypes[, c("clonotype_id", "cdr3s_aa")])
  TCR <- TCR[, c(2, 1, 3)]
  rownames(TCR) <- TCR[, 1]
  TCR[, 1] <- NULL
  if (bcrtcr == "TCR") {
    names(TCR) <- c("TCR_clonotype_id", "TCR_cdr3s_aa")
  } else {
    names(TCR) <- c("BCR_clonotype_id", "BCR_cdr3s_aa")
  }
  return(TCR)
}
