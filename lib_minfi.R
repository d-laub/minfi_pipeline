suppressPackageStartupMessages(library(minfi))
suppressPackageStartupMessages(library(sva))
suppressPackageStartupMessages(library(data.table))

process_data <- function(idat_dir, pheno_file, sex_chr=TRUE, sex_col='sex') {
    pheno = fread(pheno_file)

    message("Processing ", nrow(pheno), " samples.")
    rg_set <- read.metharray.exp(idat_dir, target=pheno, recursive=TRUE, force=TRUE)
    message("Found ", nrow(rg_set), " CpGs.")

    # Filter bad probes
    message("Removing probes that have a poor signal in any sample...")
    rg_pvals <- detectionP(rg_set)
    valid <- rowSums(rg_pvals > 0.01) == 0L
    rg_set <- rg_set[valid, ]
    message("Removed ", sum(!valid), " probes.")

    b4 = nrow(rg_set)
    message('Removing probes that have no variance across samples...')
    rg_set <- rg_set[rowVars(getGreen(rg_set)) != 0, ]
    rg_set <- rg_set[rowVars(getRed(rg_set)) != 0, ]
    message('Removed ', b4 - nrow(rg_set), ' probes.')

    # Preprocess
    message("Removing SNP-associated, cross-reactive, and non-CpG probes...")
    mset <- preprocessNoob(rg_set)
    gmset <- mapToGenome(mset)
    b4 = nrow(gmset)
    gmset <- dropLociWithSnps(gmset, snps=c("SBE", "CpG"))
    gmset <- dropMethylationLoci(gmset, dropRS=TRUE, dropCH=TRUE)
    gmset <- maxprobes::dropXreactiveLoci(gmset)
    message("Removed ", b4 - nrow(gmset), ' probes.')
    
    message('Removing samples with poor overall signal (see minfi::getQC)...')
    mset_raw <- preprocessRaw(rg_set)
    gmset_raw <- mapToGenome(mset_raw)
    qc_raw = getQC(gmset_raw)
    good_samples <- qc_raw$mMed + qc_raw$uMed > 21
    gmset_raw <- gmset_raw[, good_samples]
    gmset <- gmset[, good_samples]
    message('Removed ', sum(!good_samples), ' samples.')

    if (sex_col %chin% names(colData(gmset))) {
        message('Found sex column (', sex_col, ') in sample phenotype file.\n',
        'Removing samples with incorrectly predicted sex (see minfi::getSex)...')
        pred_sex <- getSex(gmset_raw)
        b4 = ncol(gmset)
        gmset = gmset[, colData(gmset)$sex == pred_sex$predictedSex]
        message('Removed ', b4 - ncol(gmset), ' samples.')
    }

    if (sex_chr) {
        message('Removing probes from sex chromosomes...')
        b4 = nrow(gmset)
        gmset <- gmset[seqnames(gmset) != "chrX", ]
        gmset <- gmset[seqnames(gmset) != "chrY", ]
        message("Removed ", b4 - nrow(gmset), ' probes.')
    }

    gmset
}

process_both <- function(idat_dir) {
    gr_450k <- preprocess_data(file.path(idat_dir, "450k"))
    gr_epic <- preprocess_data(file.path(idat_dir, "EPIC"))
    combineArrays(gr_450k, gr_epic)
}

cg_annotations <- function(gmset) {
    as.data.table(getAnnotation(gmset))
}

m_values <- function(gmset) {
    dt <- as.data.table(getM(gmset), keep.rownames="cg")
    setnames(dt, names(dt)[-1], colData(gmset)$Basename)
    dt
}