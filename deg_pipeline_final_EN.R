# =============================================================================
# COMPLETE DEG PIPELINE — mRNA & DNA Methylation
# TCGA Gliomas (GBM + LGG) — MOFA Cohort (318 patients)
# =============================================================================
#
# PART A — DATA PREPARATION (from raw TCGA data)
#   STEP 1  Download TCGA data (TCGAbiolinks)
#   STEP 2  Extract matrices + barcode → submitter_id deduplication
#   STEP 3  Load metadata (318 patients, MOFA cohort)
#   STEP 4  ★ Strict filtering to metadata patients ★
#   STEP 5  Missing values (methylation: CpGs with ≥ 90% NA removed)
#   STEP 6  Outlier removal (mRNA: low-count filter; methylation: IQR)
#   STEP 7  Gene type split (coding / lncRNA / miRNA)
#   STEP 8  Transformation (mRNA: TMM+voom→log2-CPM; methylation: beta→M)
#   STEP 9  Variance-based feature selection
#   STEP 10 Final z-score scaling
#   STEP 11 Save processed assays
#
# PART B — DIFFERENTIAL ANALYSIS
#   STEP 12 Load data for DEG (raw counts + beta values)
#   STEP 13 Filter to 307 classified patients (remove Unknown)
#   STEP 14 DMP — differentially methylated positions (ChAMP)
#   STEP 15 DEG — differentially expressed genes (edgeR glmQLF)
#   STEP 16 Final summary
#
# NOTE: PART A produces files for MOFA and other downstream analyses.
#       PART B restarts from RAW data (counts + beta values) because edgeR
#       and champ.DMP each require their own internal normalisation.
#       The aligned raw files saved at STEP 11 ensure cohort consistency.
# =============================================================================

setwd("D:/DEG analysis")

# =============================================================================
# LIBRARIES
# =============================================================================
library(TCGAbiolinks)
library(SummarizedExperiment)
library(edgeR)
library(limma)
library(ChAMP)
library(dplyr)
library(ggplot2)
library(reshape2)

dir.create("DataSets/processed_assays", recursive = TRUE, showWarnings = FALSE)
dir.create("dge-results",               recursive = TRUE, showWarnings = FALSE)
dir.create("info_data",                 recursive = TRUE, showWarnings = FALSE)


# =============================================================================
# ░░░░░░░░░░░░░░░░░░░░░░  PART A — DATA PREPARATION  ░░░░░░░░░░░░░░░░░░░░░░░
# =============================================================================

# =============================================================================
# STEP 1 — DOWNLOAD TCGA DATA
# =============================================================================
cat("\n========== STEP 1: Download TCGA data ==========\n")

# ---- mRNA (raw counts, STAR aligner) ---------------------------------------
query_mrna <- GDCquery(
  project       = c("TCGA-GBM", "TCGA-LGG"),
  data.category = "Transcriptome Profiling",
  data.type     = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  access        = "open"
)
GDCdownload(query_mrna)
se_mrna <- GDCprepare(query_mrna, summarizedExperiment = TRUE)

# ---- DNA methylation (beta values, 450k array) -----------------------------
query_methy <- GDCquery(
  project       = c("TCGA-GBM", "TCGA-LGG"),
  data.category = "DNA Methylation",
  data.type     = "Methylation Beta Value",
  platform      = "Illumina Human Methylation 450",
  access        = "open"
)
GDCdownload(query_methy)
se_methy <- GDCprepare(query_methy, summarizedExperiment = TRUE)


# =============================================================================
# STEP 2 — EXTRACT MATRICES + BARCODE → SUBMITTER_ID DEDUPLICATION
# =============================================================================
# TCGA barcode format: TCGA-XX-XXXX-XXX-...
# Submitter ID = first 12 characters = unique patient-level identifier
# Multiple barcodes can map to the same patient → keep the first occurrence only

cat("\n========== STEP 2: Barcode deduplication ==========\n")

barcode_to_submitter <- function(barcodes) substr(barcodes, 1, 12)

deduplicate_by_submitter <- function(mat) {
  # mat: samples (rows) × features (columns)
  submitter_ids <- barcode_to_submitter(rownames(mat))
  is_dup        <- duplicated(submitter_ids)
  cat("  Duplicates removed:", sum(is_dup), "\n")
  mat_clean <- mat[!is_dup, , drop = FALSE]
  rownames(mat_clean) <- submitter_ids[!is_dup]
  return(mat_clean)
}

# mRNA: genes × samples in the SE → transpose to samples × genes
cat("--- mRNA ---\n")
raw_counts <- assay(se_mrna, "unstranded")    # genes × samples
assay_mrna <- t(raw_counts)                    # samples × genes
cat("  Before dedup:", nrow(assay_mrna), "samples\n")
assay_mrna <- deduplicate_by_submitter(assay_mrna)
cat("  After  dedup:", nrow(assay_mrna), "samples,", ncol(assay_mrna), "genes\n")
info_mrna  <- as.data.frame(rowData(se_mrna))  # gene-level metadata

# Methylation: CpGs × samples in the SE → transpose to samples × CpGs
cat("--- Methylation ---\n")
beta_mat    <- assay(se_methy)                 # CpGs × samples
assay_methy <- t(beta_mat)                     # samples × CpGs
cat("  Before dedup:", nrow(assay_methy), "samples\n")
assay_methy <- deduplicate_by_submitter(assay_methy)
cat("  After  dedup:", nrow(assay_methy), "samples,", ncol(assay_methy), "CpGs\n")
info_methy  <- as.data.frame(rowData(se_methy)) # CpG-level metadata


# =============================================================================
# STEP 3 — LOAD METADATA (MOFA cohort, 318 patients)
# =============================================================================
# The metadata is the single source of truth for which patients are retained.
# It encodes all upstream QC decisions: deduplication, age filter (≥ 18),
# and intersection across all omics layers used in MOFA.

cat("\n========== STEP 3: Load metadata ==========\n")

metadata <- readRDS("DataSets/Metadata.rds")

cat("  Patients in metadata:", nrow(metadata), "\n")
cat("  Subtype distribution:\n"); print(table(metadata$Subtype, useNA = "ifany"))
cat("  Type distribution:\n");    print(table(metadata$Type,    useNA = "ifany"))


# =============================================================================
# STEP 4 — ★ STRICT FILTERING TO METADATA PATIENTS ★
# =============================================================================
# Only patients present in the metadata (MOFA cohort) are kept.
# This is mandatory for consistency between MOFA and the DEG analysis:
# both analyses must operate on the exact same patient population.

cat("\n========== STEP 4: Filter to metadata patients ==========\n")

filter_to_metadata <- function(mat, patients, omic_name) {
  # mat     : samples × features (rownames = patient IDs)
  # patients: vector of retained IDs (rownames of metadata)
  common  <- intersect(rownames(mat), patients)
  removed <- setdiff(rownames(mat), patients)   # in assay, not in metadata
  missing <- setdiff(patients, rownames(mat))   # in metadata, not in assay

  cat("---", omic_name, "---\n")
  cat("  Patients in assay              :", nrow(mat),        "\n")
  cat("  Patients in metadata           :", length(patients), "\n")
  cat("  Common (kept)                  :", length(common),   "\n")
  cat("  In assay but NOT in metadata   :", length(removed),  "(removed)\n")
  cat("  In metadata but NOT in assay   :", length(missing),  "(missing)\n")

  mat_filtered <- mat[common, , drop = FALSE]
  # Reorder rows to match metadata order for downstream consistency
  mat_filtered <- mat_filtered[intersect(patients, common), , drop = FALSE]
  return(mat_filtered)
}

assay_mrna  <- filter_to_metadata(assay_mrna,  rownames(metadata), "mRNA")
assay_methy <- filter_to_metadata(assay_methy, rownames(metadata), "Methylation")

# Final intersection across both assays and metadata
common_patients <- Reduce(intersect, list(
  rownames(assay_mrna),
  rownames(assay_methy),
  rownames(metadata)
))
cat("\n  Patients common to both assays + metadata:", length(common_patients), "\n")

assay_mrna  <- assay_mrna[common_patients,  , drop = FALSE]
assay_methy <- assay_methy[common_patients, , drop = FALSE]
metadata    <- metadata[common_patients,    , drop = FALSE]

# Alignment verification
stopifnot(identical(rownames(assay_mrna),  rownames(metadata)))
stopifnot(identical(rownames(assay_methy), rownames(metadata)))
cat("  ✓ Alignment verified: all patient IDs match\n")
cat("  ✓ mRNA        :", dim(assay_mrna),  "\n")
cat("  ✓ Methylation :", dim(assay_methy), "\n")


# =============================================================================
# STEP 5 — MISSING VALUES
# =============================================================================
# mRNA  : count data → no true NAs expected
# Methy : remove CpG sites with ≥ 90% NAs across samples

cat("\n========== STEP 5: Missing values ==========\n")

na_filter_features <- function(mat, threshold = 0.90, label = "") {
  na_prop <- colMeans(is.na(mat))
  keep    <- na_prop < threshold
  cat(" ", label, "- removed:", sum(!keep),
      "features (≥", threshold * 100, "% NA) | retained:", sum(keep), "\n")
  return(mat[, keep, drop = FALSE])
}

assay_methy <- na_filter_features(assay_methy, 0.90, "Methylation")


# =============================================================================
# STEP 6 — OUTLIER REMOVAL
# =============================================================================
cat("\n========== STEP 6: Outlier removal ==========\n")

# ---- mRNA: low-count gene filter -------------------------------------------
# Retain genes with counts > 1 in at least 5 samples
cat("  mRNA: low-count filter\n")
count_mat  <- t(assay_mrna)                      # genes × samples
keep_genes <- rowSums(count_mat > 1) >= 5
cat("  Genes removed:", sum(!keep_genes), "| Retained:", sum(keep_genes), "\n")
assay_mrna <- t(count_mat[keep_genes, ])          # back to samples × genes

# ---- Methylation: IQR filter (±3×IQR per CpG site) ------------------------
# Remove CpG sites whose distribution contains extreme values
cat("  Methylation: IQR filter\n")
iqr_filter_features <- function(mat) {
  Q1      <- apply(mat, 2, quantile, 0.25, na.rm = TRUE)
  Q3      <- apply(mat, 2, quantile, 0.75, na.rm = TRUE)
  IQR_val <- Q3 - Q1
  lower   <- Q1 - 3 * IQR_val
  upper   <- Q3 + 3 * IQR_val
  has_outlier <- sapply(seq_len(ncol(mat)), function(j) {
    any(mat[, j] < lower[j] | mat[, j] > upper[j], na.rm = TRUE)
  })
  cat("  CpG sites with outliers:", sum(has_outlier),
      "| Retained:", sum(!has_outlier), "\n")
  return(mat[, !has_outlier, drop = FALSE])
}
assay_methy <- iqr_filter_features(assay_methy)


# =============================================================================
# STEP 7 — GENE TYPE SPLIT
# =============================================================================
# Split mRNA assay into protein-coding, lncRNA, and miRNA
# using gene_type from the SummarizedExperiment rowData

cat("\n========== STEP 7: Gene type split ==========\n")

gene_type_map <- setNames(info_mrna$gene_type, rownames(info_mrna))
gene_types    <- gene_type_map[colnames(assay_mrna)]

idx_coding <- which(gene_types == "protein_coding")
idx_lncrna <- which(gene_types == "lncRNA")
idx_mirna  <- which(gene_types == "miRNA")

cat("  Protein-coding:", length(idx_coding), "\n")
cat("  lncRNA        :", length(idx_lncrna), "\n")
cat("  miRNA         :", length(idx_mirna),  "\n")

assay_rna_coding <- assay_mrna[, idx_coding, drop = FALSE]
assay_rna_lnc    <- assay_mrna[, idx_lncrna, drop = FALSE]
assay_rna_mirna  <- assay_mrna[, idx_mirna,  drop = FALSE]


# =============================================================================
# STEP 8 — TRANSFORMATION
# =============================================================================
cat("\n========== STEP 8: Transformation ==========\n")

# ---- 8A. mRNA: TMM normalisation + voom → log2-CPM -------------------------
# TMM corrects for differences in sequencing depth and composition bias.
# voom models the mean-variance relationship and produces precision weights.
normalize_rnaseq <- function(count_mat_s_x_g, label = "") {
  # Input : samples × genes (integer counts)
  # Output: samples × genes (log2-CPM, voom-stabilised)
  mat <- t(count_mat_s_x_g)                        # genes × samples for edgeR
  dge <- DGEList(counts = mat)
  dge <- calcNormFactors(dge, method = "TMM")       # TMM: sequencing depth correction
  v   <- voom(dge, plot = FALSE)                    # voom: variance stabilisation
  cat(" ", label, "→ log2-CPM:", nrow(v$E), "genes ×", ncol(v$E), "samples\n")
  return(t(v$E))                                    # back to samples × genes
}

assay_rna_coding_norm <- normalize_rnaseq(assay_rna_coding, "mRNA coding")
assay_rna_lnc_norm    <- normalize_rnaseq(assay_rna_lnc,    "lncRNA")
assay_rna_mirna_norm  <- normalize_rnaseq(assay_rna_mirna,  "miRNA")

# ---- 8B. Methylation: beta values → M-values --------------------------------
# M = log2(β / (1−β))
# Small epsilon added to avoid log(0) when β = 0 or β = 1
beta_to_mvalue <- function(mat, epsilon = 1e-6) {
  mat_adj <- pmin(pmax(mat, epsilon), 1 - epsilon)  # clamp to (ε, 1−ε)
  log2(mat_adj / (1 - mat_adj))
}
assay_methy_mval <- beta_to_mvalue(assay_methy)
cat("  Methylation M-values:", dim(assay_methy_mval), "\n")


# =============================================================================
# STEP 9 — VARIANCE-BASED FEATURE SELECTION
# =============================================================================
# Thresholds:
#   Protein-coding mRNA / lncRNA : top 50% most variable genes
#   miRNA                        : top 80% most variable
#   Methylation                  : top  2% most variable CpG sites

cat("\n========== STEP 9: Variance-based feature selection ==========\n")

variance_filter <- function(mat, top_fraction, label = "") {
  vars      <- apply(mat, 2, var, na.rm = TRUE)
  threshold <- quantile(vars, 1 - top_fraction, na.rm = TRUE)
  keep      <- vars >= threshold
  cat(" ", label, "→ top", top_fraction * 100, "%:",
      sum(keep), "/", length(keep), "features retained\n")
  return(mat[, keep, drop = FALSE])
}

assay_rna_coding_filt <- variance_filter(assay_rna_coding_norm, 0.50, "mRNA coding")
assay_rna_lnc_filt    <- variance_filter(assay_rna_lnc_norm,    0.50, "lncRNA")
assay_rna_mirna_filt  <- variance_filter(assay_rna_mirna_norm,  0.80, "miRNA")
assay_methy_filt      <- variance_filter(assay_methy_mval,      0.02, "Methylation")


# =============================================================================
# STEP 10 — FINAL Z-SCORE SCALING
# =============================================================================
# Centre (mean = 0) and scale (sd = 1) each feature independently.
# Applied to mRNA and methylation only (not mutations).

cat("\n========== STEP 10: Z-score scaling ==========\n")

scale_assay <- function(mat, label = "") {
  mat_scaled <- scale(mat, center = TRUE, scale = TRUE)
  cat(" ", label, "→ mean ≈",
      round(mean(mat_scaled, na.rm = TRUE), 4),
      "| sd ≈", round(sd(mat_scaled, na.rm = TRUE), 4), "\n")
  return(mat_scaled)
}

assay_rna_coding_final <- scale_assay(assay_rna_coding_filt, "mRNA coding")
assay_rna_lnc_final    <- scale_assay(assay_rna_lnc_filt,    "lncRNA")
assay_rna_mirna_final  <- scale_assay(assay_rna_mirna_filt,  "miRNA")
assay_methy_final      <- scale_assay(assay_methy_filt,      "Methylation")


# =============================================================================
# STEP 11 — SAVE PROCESSED ASSAYS
# =============================================================================
# Processed (z-scored) assays → for MOFA and other downstream analyses
# Aligned raw data               → reused in Part B for DEG/DMP

cat("\n========== STEP 11: Save processed assays ==========\n")

write.csv(assay_rna_coding_final, "DataSets/processed_assays/assay_rna_coding.csv",    row.names = TRUE)
write.csv(assay_rna_lnc_final,    "DataSets/processed_assays/assay_rna_lnc.csv",       row.names = TRUE)
write.csv(assay_rna_mirna_final,  "DataSets/processed_assays/assay_rna_mirna.csv",     row.names = TRUE)
write.csv(assay_methy_final,      "DataSets/processed_assays/assay_methylation.csv",   row.names = TRUE)
write.csv(info_mrna,              "DataSets/processed_assays/info_rna_coding.csv",     row.names = TRUE)
write.csv(info_methy,             "DataSets/processed_assays/info_methylation.csv",    row.names = TRUE)
write.csv(metadata,               "DataSets/processed_assays/metadata_final.csv",      row.names = TRUE)

# Raw counts and beta values aligned to the MOFA cohort → used in Part B
# (genes × patients and CpGs × patients format for edgeR / champ.DMP)
write.csv(t(assay_mrna),  "DataSets/processed_assays/raw_counts_coding_aligned.csv",  row.names = TRUE)
write.csv(t(assay_methy), "DataSets/processed_assays/beta_values_aligned.csv",        row.names = TRUE)

cat("  ✓ Files saved to DataSets/processed_assays/\n")
cat(sprintf("  %-22s %d patients × %d features\n", "mRNA coding:",  nrow(assay_rna_coding_final), ncol(assay_rna_coding_final)))
cat(sprintf("  %-22s %d patients × %d features\n", "lncRNA:",       nrow(assay_rna_lnc_final),    ncol(assay_rna_lnc_final)))
cat(sprintf("  %-22s %d patients × %d features\n", "miRNA:",        nrow(assay_rna_mirna_final),  ncol(assay_rna_mirna_final)))
cat(sprintf("  %-22s %d patients × %d features\n", "Methylation:",  nrow(assay_methy_final),      ncol(assay_methy_final)))

stopifnot(identical(rownames(assay_rna_coding_final), rownames(metadata)))
stopifnot(identical(rownames(assay_methy_final),      rownames(metadata)))
cat("  ✓ Final alignment confirmed\n")


# =============================================================================
# ░░░░░░░░░░░░░░░░░░░░  PART B — DIFFERENTIAL ANALYSIS  ░░░░░░░░░░░░░░░░░░░░
# =============================================================================
# IMPORTANT: Part B starts from RAW aligned data, not from the processed assays.
#   - edgeR (mRNA DEG) : requires integer raw counts, not log-CPM
#   - champ.DMP (DMP)  : requires beta values in (0,1), not M-values
# The files raw_counts_coding_aligned.csv and beta_values_aligned.csv saved
# at STEP 11 provide exactly these, already restricted to the MOFA cohort.

cat("\n\n========== PART B: DIFFERENTIAL ANALYSIS ==========\n")


# =============================================================================
# STEP 12 — LOAD DATA FOR DEG
# =============================================================================
cat("\n========== STEP 12: Load DEG data ==========\n")

# ---- Option A: running directly after Part A (objects already in memory) ---
# assay_mrna  : samples × genes  (raw counts, filtered to MOFA cohort)
# assay_methy : samples × CpGs   (beta values, filtered to MOFA cohort)
# metadata    : 318-patient metadata data frame

# ---- Option B: restarting the script from Part B only ----------------------
# Uncomment the lines below:
# assay_mrna  <- as.matrix(read.csv("DataSets/processed_assays/raw_counts_coding_aligned.csv", row.names = 1))
# assay_methy <- as.matrix(read.csv("DataSets/processed_assays/beta_values_aligned.csv",       row.names = 1))
# metadata    <- read.csv("DataSets/processed_assays/metadata_final.csv", row.names = 1)
# assay_mrna  <- t(assay_mrna)    # → samples × genes
# assay_methy <- t(assay_methy)   # → samples × CpGs

cat("  mRNA (raw counts) :", nrow(assay_mrna),  "samples ×", ncol(assay_mrna),  "genes\n")
cat("  Methylation (beta):", nrow(assay_methy), "samples ×", ncol(assay_methy), "CpGs\n")
cat("  Metadata          :", nrow(metadata),    "patients\n")


# =============================================================================
# STEP 13 — FILTER TO 307 CLASSIFIED PATIENTS (remove Unknown)
# =============================================================================
# 11 patients in the metadata have no WHO 2021 subtype assigned (Unknown).
# Including them in the DEG would introduce an unclassified group that would
# corrupt the GBM/ASTRO/OLIGO design matrix and group comparisons.
# We work on 307 patients: 143 ASTRO + 84 OLIGO + 80 GBM.

cat("\n========== STEP 13: Remove Unknown → 307 classified patients ==========\n")

keep_classified <- metadata$Subtype %in% c("GBM", "ASTRO", "OLIGO")
metadata_deg    <- metadata[keep_classified, ]
cat("  Patients after removing Unknown:", nrow(metadata_deg), "\n")
cat("  Final distribution:\n"); print(table(metadata_deg$Subtype))

# Align both assays to these 307 patients
# Note: champ.DMP expects CpGs × samples → transpose from samples × CpGs
patients_deg <- rownames(metadata_deg)

mrna_deg  <- t(assay_mrna[patients_deg,  ])   # genes × patients (for edgeR)
methy_deg <- t(assay_methy[patients_deg, ])   # CpGs  × patients (for champ.DMP)

cat("  mRNA  for DEG:", nrow(mrna_deg),  "genes  ×", ncol(mrna_deg),  "patients\n")
cat("  Methy for DMP:", nrow(methy_deg), "CpGs   ×", ncol(methy_deg), "patients\n")

# Labels aligned to columns (= patients)
labels <- metadata_deg$Subtype
stopifnot(identical(colnames(mrna_deg),  patients_deg))
stopifnot(identical(colnames(methy_deg), patients_deg))
cat("  ✓ Patient / label alignment verified\n")


# =============================================================================
# STEP 14 — DMP: DIFFERENTIALLY METHYLATED POSITIONS (ChAMP)
# =============================================================================
# champ.DMP arguments:
#   beta         : CpGs × samples matrix (values in 0-1)
#   pheno        : label vector aligned to columns of beta
#   compare.group = c("A", "B") → result stored in $<B>_to_<A>
#     logFC > 0 means hypermethylated in A relative to B
#   adjPVal = 1  → no internal pre-filtering; thresholds applied afterwards
#
# Thresholds applied: adj.P.Val (BH) < 0.05 AND |logFC| > 0.3
#   logFC here is the difference in M-values ≈ delta-beta ~ 0.10

cat("\n========== STEP 14: DMP Methylation (ChAMP) ==========\n")

dmp_gbm_vs_astro <- champ.DMP(
  beta          = methy_deg,
  pheno         = labels,
  compare.group = c("GBM", "ASTRO"),   # logFC > 0 = hypermethylated in GBM
  adjPVal       = 1,
  adjust.method = "BH",
  arraytype     = "450K"
)

dmp_gbm_vs_oligo <- champ.DMP(
  beta          = methy_deg,
  pheno         = labels,
  compare.group = c("GBM", "OLIGO"),   # logFC > 0 = hypermethylated in GBM
  adjPVal       = 1,
  adjust.method = "BH",
  arraytype     = "450K"
)

dmp_astro_vs_oligo <- champ.DMP(
  beta          = methy_deg,
  pheno         = labels,
  compare.group = c("ASTRO", "OLIGO"), # logFC > 0 = hypermethylated in ASTRO
  adjPVal       = 1,
  adjust.method = "BH",
  arraytype     = "450K"
)

# ---- Result extraction and annotation --------------------------------------
# Result slot: $<group2>_to_<group1>
#   compare.group = c("GBM", "ASTRO") → $ASTRO_to_GBM
#   logFC > 0 = hypermethylated in GBM (group1) vs ASTRO (group2)

annotate_dmp <- function(df, group1, group2) {
  df$comparison     <- paste0(group1, "_vs_", group2)
  df$diffmethylated <- "NO"
  df$diffmethylated[!is.na(df$adj.P.Val) & df$adj.P.Val < 0.05 & df$logFC >  0.3] <- "HYPERme"
  df$diffmethylated[!is.na(df$adj.P.Val) & df$adj.P.Val < 0.05 & df$logFC < -0.3] <- "HYPOme"
  cat(sprintf("  %-22s : %5d HYPERme | %5d HYPOme | %6d non-sig\n",
      paste(group1, "vs", group2),
      sum(df$diffmethylated == "HYPERme"),
      sum(df$diffmethylated == "HYPOme"),
      sum(df$diffmethylated == "NO")))
  return(df)
}

gbm_astro_methy   <- annotate_dmp(dmp_gbm_vs_astro$ASTRO_to_GBM,    "GBM",   "ASTRO")
gbm_oligo_methy   <- annotate_dmp(dmp_gbm_vs_oligo$OLIGO_to_GBM,    "GBM",   "OLIGO")
astro_oligo_methy <- annotate_dmp(dmp_astro_vs_oligo$OLIGO_to_ASTRO, "ASTRO", "OLIGO")

write.csv(gbm_astro_methy,   "dge-results/gbm_vs_astro_methy.csv",   row.names = TRUE)
write.csv(gbm_oligo_methy,   "dge-results/gbm_vs_oligo_methy.csv",   row.names = TRUE)
write.csv(astro_oligo_methy, "dge-results/astro_vs_oligo_methy.csv", row.names = TRUE)
cat("  ✓ DMP results saved\n")


# =============================================================================
# STEP 15 — DEG: DIFFERENTIALLY EXPRESSED GENES (edgeR glmQLF)
# =============================================================================
# Protocol: edgeR User Guide, GLM quasi-likelihood approach (section 4.1)
# No voom here: counts are modelled directly with a negative binomial (NB) GLM
# Thresholds: FDR (BH) < 0.05 AND |logFC| > 1 (≥ 2-fold on log2 scale)

cat("\n========== STEP 15: DEG mRNA (edgeR glmQLF) ==========\n")

# ---- 15a. Low-count gene filter --------------------------------------------
keep_genes <- rowSums(cpm(mrna_deg) > 1) >= 5
cat("  Genes after low-count filter:", sum(keep_genes), "/", nrow(mrna_deg), "\n")
mrna_filt <- mrna_deg[keep_genes, ]

# ---- 15b. DGEList + TMM normalisation --------------------------------------
dge <- DGEList(counts = mrna_filt)
dge <- normLibSizes(dge, method = "TMM")

# ---- 15c. Assign groups from metadata_deg (single source of truth) ---------
# Labels come from metadata$Subtype, not from the classification CSV
dge$samples$group <- factor(labels, levels = c("GBM", "ASTRO", "OLIGO"))
cat("  Groups in DGEList:\n"); print(table(dge$samples$group))

# ---- 15d. Design matrix without intercept ----------------------------------
design <- model.matrix(~ 0 + group, data = dge$samples)
colnames(design) <- levels(dge$samples$group)   # "GBM", "ASTRO", "OLIGO"

# ---- 15e. Dispersion estimation --------------------------------------------
dge <- estimateDisp(dge, design, robust = TRUE)
cat("  Common dispersion:", round(dge$common.dispersion, 4), "\n")

# ---- 15f. Quasi-likelihood negative binomial fit ---------------------------
fit <- glmQLFit(dge, design, robust = TRUE)

# ---- 15g. Contrasts (logFC > 0 = upregulated in group1) -------------------
con_gbm_astro   <- makeContrasts(GBM - ASTRO,  levels = design)
con_gbm_oligo   <- makeContrasts(GBM - OLIGO,  levels = design)
con_astro_oligo <- makeContrasts(ASTRO - OLIGO, levels = design)

qlf_gbm_astro   <- glmQLFTest(fit, contrast = con_gbm_astro)
qlf_gbm_oligo   <- glmQLFTest(fit, contrast = con_gbm_oligo)
qlf_astro_oligo <- glmQLFTest(fit, contrast = con_astro_oligo)

# ---- 15h. Annotation -------------------------------------------------------
annotate_deg <- function(qlf_result, group1, group2) {
  df <- topTags(qlf_result, n = Inf, sort.by = "none")$table
  df$comparison    <- paste0(group1, "_vs_", group2)
  df$diffexpressed <- "NO"
  df$diffexpressed[!is.na(df$FDR) & df$FDR < 0.05 & df$logFC >  1] <- "UP"
  df$diffexpressed[!is.na(df$FDR) & df$FDR < 0.05 & df$logFC < -1] <- "DOWN"
  cat(sprintf("  %-22s : %5d UP | %5d DOWN | %6d non-sig\n",
      paste(group1, "vs", group2),
      sum(df$diffexpressed == "UP"),
      sum(df$diffexpressed == "DOWN"),
      sum(df$diffexpressed == "NO")))
  return(df)
}

gbm_astro_mrna   <- annotate_deg(qlf_gbm_astro,   "GBM",   "ASTRO")
gbm_oligo_mrna   <- annotate_deg(qlf_gbm_oligo,   "GBM",   "OLIGO")
astro_oligo_mrna <- annotate_deg(qlf_astro_oligo, "ASTRO", "OLIGO")

write.csv(gbm_astro_mrna,   "dge-results/gbm_vs_astro_mrna.csv",   row.names = TRUE)
write.csv(gbm_oligo_mrna,   "dge-results/gbm_vs_oligo_mrna.csv",   row.names = TRUE)
write.csv(astro_oligo_mrna, "dge-results/astro_vs_oligo_mrna.csv", row.names = TRUE)
cat("  ✓ DEG results saved\n")


# =============================================================================
# STEP 16 — FINAL SUMMARY
# =============================================================================
cat("\n========== STEP 16: FINAL SUMMARY ==========\n")

cat("\n--- Methylation (DMP): adj.P.Val < 0.05, |logFC| > 0.3 ---\n")
for (res in list(
  list(df = gbm_astro_methy,   g1 = "GBM",   g2 = "ASTRO"),
  list(df = gbm_oligo_methy,   g1 = "GBM",   g2 = "OLIGO"),
  list(df = astro_oligo_methy, g1 = "ASTRO", g2 = "OLIGO")
)) {
  cat(sprintf("  %-22s : %5d HYPERme | %5d HYPOme\n",
      paste(res$g1, "vs", res$g2),
      sum(res$df$diffmethylated == "HYPERme"),
      sum(res$df$diffmethylated == "HYPOme")))
}

cat("\n--- mRNA (DEG): FDR < 0.05, |logFC| > 1 ---\n")
for (res in list(
  list(df = gbm_astro_mrna,   g1 = "GBM",   g2 = "ASTRO"),
  list(df = gbm_oligo_mrna,   g1 = "GBM",   g2 = "OLIGO"),
  list(df = astro_oligo_mrna, g1 = "ASTRO", g2 = "OLIGO")
)) {
  cat(sprintf("  %-22s : %5d UP | %5d DOWN\n",
      paste(res$g1, "vs", res$g2),
      sum(res$df$diffexpressed == "UP"),
      sum(res$df$diffexpressed == "DOWN")))
}

cat("\n=== COMPLETE DEG PIPELINE — DONE ===\n")

# =============================================================================
# OUTPUT FILES
# =============================================================================
# DataSets/processed_assays/
#   assay_rna_coding.csv            mRNA coding   log2-CPM z-scored  (for MOFA)
#   assay_rna_lnc.csv               lncRNA        log2-CPM z-scored  (for MOFA)
#   assay_rna_mirna.csv             miRNA         log2-CPM z-scored  (for MOFA)
#   assay_methylation.csv           Methylation   M-values z-scored  (for MOFA)
#   raw_counts_coding_aligned.csv   Raw counts    aligned to cohort  (for DEG)
#   beta_values_aligned.csv         Beta values   aligned to cohort  (for DEG)
#   metadata_final.csv              Metadata 318 patients
#
# dge-results/
#   gbm_vs_astro_methy.csv          DMP GBM vs ASTRO
#   gbm_vs_oligo_methy.csv          DMP GBM vs OLIGO
#   astro_vs_oligo_methy.csv        DMP ASTRO vs OLIGO
#   gbm_vs_astro_mrna.csv           DEG GBM vs ASTRO
#   gbm_vs_oligo_mrna.csv           DEG GBM vs OLIGO
#   astro_vs_oligo_mrna.csv         DEG ASTRO vs OLIGO
# =============================================================================
