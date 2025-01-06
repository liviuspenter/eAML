# calculate mitochondrial variants from ATAC-seq bulk data
# filtering steps:
# (i) coverage of position > 200x
# (ii) strand coordination > 0.6
# (iii) coverage of mutation > 100x

variant_calling_bulk <- function(sample.bulk, total.coverage = 200, strand.coordination = 0.6,
                                 coverage.position = 100, sample.name = "samplename") {
  # Determinie key coverage statistics every which way
  cov.sample <- assays(sample.bulk)[["coverage"]]
  ref_allele <- toupper(as.character(rowRanges(sample.bulk)$refAllele))

  # Process mutation for one alternate letter
  process_letter <- function(letter) {
    message(letter)
    boo <- ref_allele != letter & ref_allele != "N"
    pos <- start(rowRanges(sample.bulk))
    variant_name <- paste0(as.character(pos), ref_allele, ">", letter)[boo]
    nucleotide <- paste0(ref_allele, ">", letter)[boo]
    position_filt <- pos[boo]

    comparison.matrix <- cbind(
      Matrix::as.matrix((assays(sample.bulk)[[paste0(letter, "_counts_fw")]] + assays(sample.bulk)[[paste0(letter, "_counts_rev")]]) /
        Matrix::as.matrix(cov.sample)),
      Matrix::as.matrix(assays(sample.bulk)[[paste0(letter, "_counts_fw")]] + assays(sample.bulk)[[paste0(letter, "_counts_rev")]]),
      Matrix::as.matrix(assays(sample.bulk)[[paste0(letter, "_counts_fw")]] / assays(sample.bulk)[[paste0(letter, "_counts_rev")]]),
      Matrix::as.matrix(cov.sample)
    )
    comparison.matrix <- comparison.matrix[boo, ]
    colnames(comparison.matrix) <- c("heteroplasmy", "coverage.variant", "ratio.reads", "coverage")

    update_missing_w_zero <- function(vec) {
      ifelse(is.na(vec) | is.nan(vec), 0, vec)
    }
    comparison.matrix <- update_missing_w_zero(comparison.matrix)

    # apply filters
    insufficient.cov <- which(comparison.matrix[, 2] <= coverage.position | comparison.matrix[, 4] <= total.coverage)
    insufficient.statistics <- which(comparison.matrix[, 3] < strand.coordination | comparison.matrix[, 3] > (1 / strand.coordination))
    comparison.matrix[insufficient.cov, 1] <- NA
    comparison.matrix[insufficient.statistics, 1] <- NA

    results <- data.frame(
      variant = variant_name,
      heteroplasmy = comparison.matrix[, 1],
      coverage = comparison.matrix[, 2],
      strand.ratio = comparison.matrix[, 3]
    )

    return(results)

    # Compute per-mutation summary statistics
    var_summary_df <- data.frame(
      position = position_filt,
      nucleotide = nucleotide,
      variant = variant_name,
      strand_correlation = comparison.matrix[, 3],
      mean_coverage = comparison.matrix[, 4],
      stringsAsFactors = FALSE, row.names = variant_name
    )
    col_data <- data.frame(sample = sample.name, depth = sample.bulk$depth)
    rownames(col_data) <- sample.name

    se_new <- SummarizedExperiment(
      rowData = var_summary_df,
      colData = col_data,
      assays = list(
        heteroplasmy = as.matrix(comparison.matrix[, 1]),
        coverage = as.matrix(comparison.matrix[, 2])
      )
    )

    return(se_new)
  }

  return(rbind(process_letter("A"), process_letter("C"), process_letter("G"), process_letter("T")))
}
