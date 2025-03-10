read_clonotypes <- function(list_libraries, AML.combined, TCR.clones.file = NA, exclude.samples = NA) {
  sample.list <- readxl::read_excel(list_libraries)

  clonotypes.combined <- data.frame()
  for (sample in unique(sample.list$Sample)) {
    if (!sample %in% exclude.samples) {
      pool <- sample.list$Pool[which(sample.list$Sample == sample & sample.list$Assay == "TCR")]
      clonotypes <- data.table::fread(file = paste0("./data/", sample.list$Pool[which(sample.list$Assay == "TCR" & sample.list$Sample == sample)], "/clonotypes.csv"))
      clonotypes$sample <- sample
      clonotypes.combined <- rbind(clonotypes.combined, clonotypes)
    }
  }

  # identify clones across the samples
  TCR.clones <- data.frame(
    CloneID = as.numeric(),
    CDR3.beta1 = as.character(), CDR3.beta2 = as.character()
  )
  TCR.clones.line <- data.frame(as.numeric(0))
  for (s in unique(sample.list$Sample)) {
    if (!s %in% exclude.samples) {
      print(s)
      TCR.clones <- cbind(TCR.clones, data.frame(s = as.numeric()))
      colnames(TCR.clones)[ncol(TCR.clones)] <- s
      TCR.clones.line <- cbind(TCR.clones.line, data.frame(s = as.numeric(0)))
      colnames(TCR.clones.line)[ncol(TCR.clones.line)] <- s
    }
  }
  TCR.clones.line <- TCR.clones.line[, -1]
  if (length(unique(sample.list$Sample)) == 1) {
    TCR.clones.line <- data.frame(s = 0)
    colnames(TCR.clones.line) <- unique(sample.list$Sample)
  }

  # extract CDR3 alpha and CDR3 beta
  clonotypes.combined$CDR3.alpha1 <- NA
  clonotypes.combined$CDR3.alpha2 <- NA
  clonotypes.combined$CDR3.beta1 <- NA
  clonotypes.combined$CDR3.beta2 <- NA
  for (r in seq(1, nrow(clonotypes.combined))) {
    CDR3.splitted <- strsplit(clonotypes.combined$cdr3s_aa[r], ";")
    for (CDR3 in CDR3.splitted[[1]]) {
      CDR3.individual.type <- strsplit(CDR3[1], ":")[[1]][1]
      CDR3.individual.sequence <- strsplit(CDR3[1], ":")[[1]][2]
      if (CDR3.individual.type == "TRA" & is.na(clonotypes.combined$CDR3.alpha1[r])) {
        clonotypes.combined$CDR3.alpha1[r] <- CDR3.individual.sequence
      } else if (CDR3.individual.type == "TRA" & is.na(clonotypes.combined$CDR3.alpha2[r])) {
        clonotypes.combined$CDR3.alpha2[r] <- CDR3.individual.sequence
      } else if (CDR3.individual.type == "TRB" & is.na(clonotypes.combined$CDR3.beta1[r])) {
        clonotypes.combined$CDR3.beta1[r] <- CDR3.individual.sequence
      } else if (CDR3.individual.type == "TRB" & is.na(clonotypes.combined$CDR3.beta2[r])) {
        clonotypes.combined$CDR3.beta2[r] <- CDR3.individual.sequence
      }
    }
  }

  # assess all clones with 2x CDR3 beta
  cloneid <- 1
  clonotypes.combined$CloneID <- 0

  for (r in which(!is.na(clonotypes.combined$CDR3.beta2))) {
    CDR3.beta1 <- sort(c(clonotypes.combined$CDR3.beta1[r], clonotypes.combined$CDR3.beta2[r]))[1]
    CDR3.beta2 <- sort(c(clonotypes.combined$CDR3.beta1[r], clonotypes.combined$CDR3.beta2[r]))[2]
    if (length(which(TCR.clones$CDR3.beta1 == CDR3.beta1 & TCR.clones$CDR3.beta2 == CDR3.beta2)) > 0) {
      clone <- as.numeric(which(TCR.clones$CDR3.beta1 == CDR3.beta1 & TCR.clones$CDR3.beta2 == CDR3.beta2))
      TCR.clones[clone, clonotypes.combined$sample[r]] <- TCR.clones[clone, clonotypes.combined$sample[r]] +
        clonotypes.combined$frequency[r]
      clonotypes.combined$CloneID[r] <- TCR.clones$CloneID[clone]
    } else {
      TCR.clones <- rbind(TCR.clones, data.frame(CloneID = cloneid, CDR3.beta1 = CDR3.beta1, CDR3.beta2 = CDR3.beta2, TCR.clones.line))
      TCR.clones[nrow(TCR.clones), clonotypes.combined$sample[r]] <- clonotypes.combined$frequency[r]
      clonotypes.combined$CloneID[r] <- cloneid
      cloneid <- cloneid + 1
    }
  }

  # assess all clones with 1x CDR3 beta
  for (r in which(clonotypes.combined$CloneID == 0)) {
    CDR3.beta1 <- sort(c(clonotypes.combined$CDR3.beta1[r], clonotypes.combined$CDR3.beta2[r]))[1]
    CDR3.beta2 <- sort(c(clonotypes.combined$CDR3.beta1[r], clonotypes.combined$CDR3.beta2[r]))[2]
    clone <- which(TCR.clones$CDR3.beta1 == clonotypes.combined$CDR3.beta1[r] | TCR.clones$CDR3.beta2 == clonotypes.combined$CDR3.beta1[r])
    if (length(clone) == 1) {
      clonotypes.combined$CloneID[r] <- TCR.clones$CloneID[clone]
      TCR.clones[clone, clonotypes.combined$sample[r]] <- TCR.clones[clone, clonotypes.combined$sample[r]] +
        clonotypes.combined$frequency[r]
    } else {
      TCR.clones <- rbind(TCR.clones, data.frame(CloneID = cloneid, CDR3.beta1 = CDR3.beta1, CDR3.beta2 = CDR3.beta2, TCR.clones.line))
      TCR.clones[nrow(TCR.clones), clonotypes.combined$sample[r]] <- clonotypes.combined$frequency[r]
      clonotypes.combined$CloneID[r] <- cloneid
      cloneid <- cloneid + 1
    }
  }

  for (s in unique(sample.list$Sample)) {
    if (!s %in% exclude.samples) {
      TCR.clones <- cbind(TCR.clones, TCR.clones[, s] / sum(TCR.clones[, s]))
      colnames(TCR.clones)[ncol(TCR.clones)] <- paste0(s, ".freq")
    }
  }

  AML.combined$CloneID <- NA
  AML.combined$clonotype.freq <- NA
  AML.combined$clonotype.cells <- NA
  for (c in seq(1, nrow(clonotypes.combined))) {
    index <- paste0(clonotypes.combined$sample[c], "_", clonotypes.combined$clonotype_id[c])
    AML.combined$CloneID[which(AML.combined$TCR_clonotype_id == index)] <- clonotypes.combined$CloneID[c]
    AML.combined$clonotype.freq[which(AML.combined$TCR_clonotype_id == index)] <- clonotypes.combined$proportion[c]
    AML.combined$clonotype.cells[which(AML.combined$TCR_clonotype_id == index)] <- clonotypes.combined$frequency[c]
  }

  write.table(TCR.clones, TCR.clones.file, sep = "\t", quote = F, row.names = F)

  AML.combined
}
