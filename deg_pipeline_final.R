# =============================================================================
# PIPELINE DEG COMPLET — mRNA & DNA Methylation
# Gliomes TCGA (GBM + LGG) — Cohorte MOFA (318 patients)
# =============================================================================
#
# PARTIE A — PRÉPARATION DES DONNÉES (depuis TCGA brut)
#   STEP 1  Téléchargement TCGA (TCGAbiolinks)
#   STEP 2  Extraction matrices + déduplication barcode → submitter_id
#   STEP 3  Chargement metadata (318 patients, cohorte MOFA)
#   STEP 4  ★ Filtrage strict aux patients metadata ★
#   STEP 5  Valeurs manquantes (méthylation : CpGs ≥ 90% NA supprimés)
#   STEP 6  Outliers (mRNA : low-count ; méthylation : IQR)
#   STEP 7  Séparation par type de gène (coding / lncRNA / miRNA)
#   STEP 8  Transformation (mRNA : TMM+voom→log2-CPM ; méthylation : beta→M)
#   STEP 9  Sélection de features par variance
#   STEP 10 Normalisation z-score finale
#   STEP 11 Sauvegarde des assays processés
#
# PARTIE B — ANALYSE DIFFÉRENTIELLE
#   STEP 12 Chargement des données pour le DEG (raw counts + beta values)
#   STEP 13 Filtrage aux 307 patients classifiés (sans Unknown)
#   STEP 14 DMP — positions différentiellement méthylées (ChAMP)
#   STEP 15 DEG — gènes différentiellement exprimés (edgeR glmQLF)
#   STEP 16 Résumé final
#
# NOTE : PARTIE A produit les fichiers pour MOFA et d'autres analyses.
#        PARTIE B repart des données BRUTES (counts + beta) car edgeR
#        et champ.DMP ont leur propre normalisation interne.
# =============================================================================

setwd("D:/DEG analysis")

# =============================================================================
# LIBRAIRIES
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
# ░░░░░░░░░░░░░░░  PARTIE A — PRÉPARATION DES DONNÉES  ░░░░░░░░░░░░░░░░░░░░░
# =============================================================================

# =============================================================================
# STEP 1 — TÉLÉCHARGEMENT TCGA
# =============================================================================
cat("\n========== STEP 1: Téléchargement TCGA ==========\n")

query_mrna <- GDCquery(
  project       = c("TCGA-GBM", "TCGA-LGG"),
  data.category = "Transcriptome Profiling",
  data.type     = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  access        = "open"
)
GDCdownload(query_mrna)
se_mrna <- GDCprepare(query_mrna, summarizedExperiment = TRUE)

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
# STEP 2 — EXTRACTION + DÉDUPLICATION BARCODE → SUBMITTER_ID
# =============================================================================
# Barcode TCGA : TCGA-XX-XXXX-XXX-...
# Submitter ID = 12 premiers caractères = identifiant patient unique
# Plusieurs barcodes peuvent correspondre au même patient → garder le premier

cat("\n========== STEP 2: Déduplication barcodes ==========\n")

barcode_to_submitter <- function(barcodes) substr(barcodes, 1, 12)

deduplicate_by_submitter <- function(mat) {
  # mat : samples (lignes) × features (colonnes)
  submitter_ids <- barcode_to_submitter(rownames(mat))
  is_dup        <- duplicated(submitter_ids)
  cat("  Doublons supprimés :", sum(is_dup), "\n")
  mat_clean <- mat[!is_dup, , drop = FALSE]
  rownames(mat_clean) <- submitter_ids[!is_dup]
  return(mat_clean)
}

# mRNA : genes × samples dans le SE → transposer en samples × genes
cat("--- mRNA ---\n")
raw_counts <- assay(se_mrna, "unstranded")    # genes × samples
assay_mrna <- t(raw_counts)                    # samples × genes
cat("  Avant dédup :", nrow(assay_mrna), "samples\n")
assay_mrna <- deduplicate_by_submitter(assay_mrna)
cat("  Après dédup :", nrow(assay_mrna), "samples,", ncol(assay_mrna), "gènes\n")
info_mrna  <- as.data.frame(rowData(se_mrna))  # métadonnées gènes

# Méthylation : CpGs × samples dans le SE → transposer en samples × CpGs
cat("--- Méthylation ---\n")
beta_mat    <- assay(se_methy)                 # CpGs × samples
assay_methy <- t(beta_mat)                     # samples × CpGs
cat("  Avant dédup :", nrow(assay_methy), "samples\n")
assay_methy <- deduplicate_by_submitter(assay_methy)
cat("  Après dédup :", nrow(assay_methy), "samples,", ncol(assay_methy), "CpGs\n")
info_methy  <- as.data.frame(rowData(se_methy)) # métadonnées CpGs


# =============================================================================
# STEP 3 — CHARGEMENT DU METADATA (cohorte MOFA, 318 patients)
# =============================================================================
cat("\n========== STEP 3: Chargement metadata ==========\n")

metadata <- readRDS("DataSets/Metadata.rds")

cat("  Patients dans le metadata :", nrow(metadata), "\n")
cat("  Distribution Subtype :\n"); print(table(metadata$Subtype, useNA = "ifany"))
cat("  Distribution Type    :\n"); print(table(metadata$Type,    useNA = "ifany"))


# =============================================================================
# STEP 4 — ★ FILTRAGE STRICT AUX PATIENTS METADATA ★
# =============================================================================
# Principe : seuls les patients présents DANS le metadata (cohorte MOFA)
# sont conservés. Garantit la cohérence entre MOFA et le DEG.

cat("\n========== STEP 4: Filtrage aux patients metadata ==========\n")

filter_to_metadata <- function(mat, patients, omic_name) {
  # mat     : samples × features (rownames = patient IDs)
  # patients: vecteur des IDs retenus (rownames du metadata)
  common    <- intersect(rownames(mat), patients)
  removed   <- setdiff(rownames(mat), patients)
  missing   <- setdiff(patients, rownames(mat))

  cat("---", omic_name, "---\n")
  cat("  Dans l'assay                :", nrow(mat),         "\n")
  cat("  Dans le metadata            :", length(patients),  "\n")
  cat("  Communs (conservés)         :", length(common),    "\n")
  cat("  Dans assay, absent metadata :", length(removed),   "(supprimés)\n")
  cat("  Dans metadata, absent assay :", length(missing),   "(manquants)\n")

  mat_filtered <- mat[common, , drop = FALSE]
  # Réordonner pour correspondre à l'ordre du metadata
  mat_filtered <- mat_filtered[intersect(patients, common), , drop = FALSE]
  return(mat_filtered)
}

assay_mrna  <- filter_to_metadata(assay_mrna,  rownames(metadata), "mRNA")
assay_methy <- filter_to_metadata(assay_methy, rownames(metadata), "Méthylation")

# Intersection finale sur les deux assays + metadata
common_patients <- Reduce(intersect, list(
  rownames(assay_mrna),
  rownames(assay_methy),
  rownames(metadata)
))
cat("\n  Patients communs aux deux assays + metadata :", length(common_patients), "\n")

assay_mrna  <- assay_mrna[common_patients,  , drop = FALSE]
assay_methy <- assay_methy[common_patients, , drop = FALSE]
metadata    <- metadata[common_patients,    , drop = FALSE]

# Vérification alignement
stopifnot(identical(rownames(assay_mrna),  rownames(metadata)))
stopifnot(identical(rownames(assay_methy), rownames(metadata)))
cat("  ✓ Alignement vérifié : tous les IDs concordent\n")
cat("  ✓ mRNA         :", dim(assay_mrna),  "\n")
cat("  ✓ Méthylation  :", dim(assay_methy), "\n")


# =============================================================================
# STEP 5 — VALEURS MANQUANTES
# =============================================================================
# mRNA  : données de comptage → pas de vrais NAs attendus
# Méthy : supprimer les CpGs avec ≥ 90% de NAs sur les samples

cat("\n========== STEP 5: Valeurs manquantes ==========\n")

na_filter_features <- function(mat, threshold = 0.90, label = "") {
  na_prop <- colMeans(is.na(mat))
  keep    <- na_prop < threshold
  cat(" ", label, "- supprimés :", sum(!keep),
      "features (≥", threshold * 100, "% NA) | conservés :", sum(keep), "\n")
  return(mat[, keep, drop = FALSE])
}

assay_methy <- na_filter_features(assay_methy, 0.90, "Méthylation")


# =============================================================================
# STEP 6 — OUTLIERS
# =============================================================================
cat("\n========== STEP 6: Outliers ==========\n")

# ---- mRNA : filtre low-count -----------------------------------------------
# Conserver les gènes avec counts > 1 dans au moins 5 samples
cat("  mRNA : filtre low-count\n")
count_mat  <- t(assay_mrna)                      # genes × samples
keep_genes <- rowSums(count_mat > 1) >= 5
cat("  Gènes supprimés :", sum(!keep_genes), "| Conservés :", sum(keep_genes), "\n")
assay_mrna <- t(count_mat[keep_genes, ])          # retour samples × genes

# ---- Méthylation : filtre IQR (±3×IQR par CpG) ----------------------------
cat("  Méthylation : filtre IQR\n")
iqr_filter_features <- function(mat) {
  Q1      <- apply(mat, 2, quantile, 0.25, na.rm = TRUE)
  Q3      <- apply(mat, 2, quantile, 0.75, na.rm = TRUE)
  IQR_val <- Q3 - Q1
  lower   <- Q1 - 3 * IQR_val
  upper   <- Q3 + 3 * IQR_val
  has_outlier <- sapply(seq_len(ncol(mat)), function(j) {
    any(mat[, j] < lower[j] | mat[, j] > upper[j], na.rm = TRUE)
  })
  cat("  CpGs avec outliers :", sum(has_outlier),
      "| Conservés :", sum(!has_outlier), "\n")
  return(mat[, !has_outlier, drop = FALSE])
}
assay_methy <- iqr_filter_features(assay_methy)


# =============================================================================
# STEP 7 — SÉPARATION PAR TYPE DE GÈNE
# =============================================================================
cat("\n========== STEP 7: Séparation par type de gène ==========\n")

gene_type_map <- setNames(info_mrna$gene_type, rownames(info_mrna))
gene_types    <- gene_type_map[colnames(assay_mrna)]

idx_coding <- which(gene_types == "protein_coding")
idx_lncrna <- which(gene_types == "lncRNA")
idx_mirna  <- which(gene_types == "miRNA")

cat("  Protein-coding :", length(idx_coding), "\n")
cat("  lncRNA         :", length(idx_lncrna), "\n")
cat("  miRNA          :", length(idx_mirna),  "\n")

assay_rna_coding <- assay_mrna[, idx_coding, drop = FALSE]
assay_rna_lnc    <- assay_mrna[, idx_lncrna, drop = FALSE]
assay_rna_mirna  <- assay_mrna[, idx_mirna,  drop = FALSE]


# =============================================================================
# STEP 8 — TRANSFORMATION
# =============================================================================
cat("\n========== STEP 8: Transformation ==========\n")

# ---- 8A. mRNA : TMM + voom → log2-CPM -------------------------------------
normalize_rnaseq <- function(count_mat_s_x_g, label = "") {
  # Entrée  : samples × genes (counts entiers)
  # Sortie  : samples × genes (log2-CPM, variance stabilisée par voom)
  mat <- t(count_mat_s_x_g)                        # genes × samples pour edgeR
  dge <- DGEList(counts = mat)
  dge <- calcNormFactors(dge, method = "TMM")       # TMM : correction profondeur
  v   <- voom(dge, plot = FALSE)                    # voom : modélisation variance
  cat(" ", label, "→ log2-CPM :", nrow(v$E), "gènes ×", ncol(v$E), "samples\n")
  return(t(v$E))                                    # retour samples × genes
}

assay_rna_coding_norm <- normalize_rnaseq(assay_rna_coding, "mRNA coding")
assay_rna_lnc_norm    <- normalize_rnaseq(assay_rna_lnc,    "lncRNA")
assay_rna_mirna_norm  <- normalize_rnaseq(assay_rna_mirna,  "miRNA")

# ---- 8B. Méthylation : beta → M-values -------------------------------------
# M = log2(β / (1−β))  avec epsilon pour éviter log(0)
beta_to_mvalue <- function(mat, epsilon = 1e-6) {
  mat_adj <- pmin(pmax(mat, epsilon), 1 - epsilon)
  log2(mat_adj / (1 - mat_adj))
}
assay_methy_mval <- beta_to_mvalue(assay_methy)
cat("  Méthylation M-values :", dim(assay_methy_mval), "\n")


# =============================================================================
# STEP 9 — SÉLECTION DE FEATURES PAR VARIANCE
# =============================================================================
# Seuils du pipeline original :
#   mRNA coding / lncRNA : top 50% les plus variables
#   miRNA                : top 80%
#   Méthylation          : top  2% (CpG sites)

cat("\n========== STEP 9: Sélection par variance ==========\n")

variance_filter <- function(mat, top_fraction, label = "") {
  vars      <- apply(mat, 2, var, na.rm = TRUE)
  threshold <- quantile(vars, 1 - top_fraction, na.rm = TRUE)
  keep      <- vars >= threshold
  cat(" ", label, "→ top", top_fraction * 100, "% :",
      sum(keep), "/", length(keep), "features conservées\n")
  return(mat[, keep, drop = FALSE])
}

assay_rna_coding_filt <- variance_filter(assay_rna_coding_norm, 0.50, "mRNA coding")
assay_rna_lnc_filt    <- variance_filter(assay_rna_lnc_norm,    0.50, "lncRNA")
assay_rna_mirna_filt  <- variance_filter(assay_rna_mirna_norm,  0.80, "miRNA")
assay_methy_filt      <- variance_filter(assay_methy_mval,      0.02, "Méthylation")


# =============================================================================
# STEP 10 — NORMALISATION Z-SCORE FINALE
# =============================================================================
# Centrage (moyenne = 0) et réduction (écart-type = 1) par feature
# Appliqué à mRNA et méthylation (pas aux mutations)

cat("\n========== STEP 10: Z-score ==========\n")

scale_assay <- function(mat, label = "") {
  mat_scaled <- scale(mat, center = TRUE, scale = TRUE)
  cat(" ", label, "→ moyenne ≈",
      round(mean(mat_scaled, na.rm = TRUE), 4),
      "| sd ≈", round(sd(mat_scaled, na.rm = TRUE), 4), "\n")
  return(mat_scaled)
}

assay_rna_coding_final <- scale_assay(assay_rna_coding_filt, "mRNA coding")
assay_rna_lnc_final    <- scale_assay(assay_rna_lnc_filt,    "lncRNA")
assay_rna_mirna_final  <- scale_assay(assay_rna_mirna_filt,  "miRNA")
assay_methy_final      <- scale_assay(assay_methy_filt,      "Méthylation")


# =============================================================================
# STEP 11 — SAUVEGARDE DES ASSAYS PROCESSÉS (pour MOFA et autres analyses)
# =============================================================================
cat("\n========== STEP 11: Sauvegarde assays processés ==========\n")

write.csv(assay_rna_coding_final, "DataSets/processed_assays/assay_rna_coding.csv",    row.names = TRUE)
write.csv(assay_rna_lnc_final,    "DataSets/processed_assays/assay_rna_lnc.csv",       row.names = TRUE)
write.csv(assay_rna_mirna_final,  "DataSets/processed_assays/assay_rna_mirna.csv",     row.names = TRUE)
write.csv(assay_methy_final,      "DataSets/processed_assays/assay_methylation.csv",   row.names = TRUE)
write.csv(info_mrna,              "DataSets/processed_assays/info_rna_coding.csv",     row.names = TRUE)
write.csv(info_methy,             "DataSets/processed_assays/info_methylation.csv",    row.names = TRUE)
write.csv(metadata,               "DataSets/processed_assays/metadata_final.csv",      row.names = TRUE)

# Sauvegarde également des raw counts + beta values alignés sur le metadata
# → seront réutilisés en Partie B pour le DEG
write.csv(t(assay_mrna),  "DataSets/processed_assays/raw_counts_coding_aligned.csv",  row.names = TRUE)
write.csv(t(assay_methy), "DataSets/processed_assays/beta_values_aligned.csv",        row.names = TRUE)

cat("  ✓ Fichiers sauvegardés dans DataSets/processed_assays/\n")
cat(sprintf("  %-22s %d patients × %d features\n", "mRNA coding :",  nrow(assay_rna_coding_final), ncol(assay_rna_coding_final)))
cat(sprintf("  %-22s %d patients × %d features\n", "lncRNA :",       nrow(assay_rna_lnc_final),    ncol(assay_rna_lnc_final)))
cat(sprintf("  %-22s %d patients × %d features\n", "miRNA :",        nrow(assay_rna_mirna_final),  ncol(assay_rna_mirna_final)))
cat(sprintf("  %-22s %d patients × %d features\n", "Méthylation :",  nrow(assay_methy_final),      ncol(assay_methy_final)))

stopifnot(identical(rownames(assay_rna_coding_final), rownames(metadata)))
stopifnot(identical(rownames(assay_methy_final),      rownames(metadata)))
cat("  ✓ Alignement final confirmé\n")


# =============================================================================
# ░░░░░░░░░░░░░░░  PARTIE B — ANALYSE DIFFÉRENTIELLE  ░░░░░░░░░░░░░░░░░░░░░░
# =============================================================================
# IMPORTANT : on repart ici des données BRUTES alignées sur le metadata.
#   - edgeR (DEG mRNA) : nécessite des raw counts entiers, pas de log-CPM
#   - champ.DMP (DMP)  : nécessite des beta values (0-1), pas des M-values
# Les fichiers raw_counts_coding_aligned.csv et beta_values_aligned.csv
# produits en STEP 11 contiennent exactement ces données, filtrées sur la
# cohorte MOFA et dédupliquées.

cat("\n\n========== PARTIE B : ANALYSE DIFFÉRENTIELLE ==========\n")


# =============================================================================
# STEP 12 — CHARGEMENT DES DONNÉES POUR LE DEG
# =============================================================================
cat("\n========== STEP 12: Chargement données DEG ==========\n")

# ---- Option A : si on enchaîne directement après la Partie A ---------------
# Les objets assay_mrna, assay_methy et metadata sont déjà en mémoire.
# assay_mrna  : samples × genes  (raw counts, filtrés sur metadata)
# assay_methy : samples × CpGs   (beta values, filtrées sur metadata)

# ---- Option B : si on relance le script séparément -------------------------
# Décommenter les lignes suivantes :
# assay_mrna  <- as.matrix(read.csv("DataSets/processed_assays/raw_counts_coding_aligned.csv",  row.names = 1))
# assay_methy <- as.matrix(read.csv("DataSets/processed_assays/beta_values_aligned.csv",        row.names = 1))
# metadata    <- read.csv("DataSets/processed_assays/metadata_final.csv", row.names = 1)
# assay_mrna  <- t(assay_mrna)   # → samples × genes
# assay_methy <- t(assay_methy)  # → samples × CpGs

cat("  mRNA (raw counts)  :", nrow(assay_mrna),  "samples ×", ncol(assay_mrna),  "gènes\n")
cat("  Méthylation (beta) :", nrow(assay_methy), "samples ×", ncol(assay_methy), "CpGs\n")
cat("  Metadata           :", nrow(metadata),    "patients\n")


# =============================================================================
# STEP 13 — FILTRAGE AUX 307 PATIENTS CLASSIFIÉS (suppression Unknown)
# =============================================================================
# Les 11 patients Unknown n'ont pas de subtype WHO assigné.
# Les inclure dans le DEG polluerait les groupes GBM/ASTRO/OLIGO.
# On travaille donc sur 307 patients (143 ASTRO + 84 OLIGO + 80 GBM).

cat("\n========== STEP 13: Filtrage Unknown → 307 patients classifiés ==========\n")

keep_classified <- metadata$Subtype %in% c("GBM", "ASTRO", "OLIGO")
metadata_deg    <- metadata[keep_classified, ]
cat("  Patients après suppression Unknown :", nrow(metadata_deg), "\n")
cat("  Distribution finale :\n"); print(table(metadata_deg$Subtype))

# Aligner les deux assays sur ces 307 patients
# méthylation : champ.DMP attend CpGs × samples → transposer
patients_deg  <- rownames(metadata_deg)

mrna_deg  <- t(assay_mrna[patients_deg,  ])   # genes   × patients (pour edgeR)
methy_deg <- t(assay_methy[patients_deg, ])   # CpGs    × patients (pour champ.DMP)

cat("  mRNA  pour DEG :", nrow(mrna_deg),  "gènes  ×", ncol(mrna_deg),  "patients\n")
cat("  Méthy pour DMP :", nrow(methy_deg), "CpGs   ×", ncol(methy_deg), "patients\n")

# Labels alignés sur les colonnes (= patients)
labels <- metadata_deg$Subtype
stopifnot(identical(colnames(mrna_deg),  patients_deg))
stopifnot(identical(colnames(methy_deg), patients_deg))
cat("  ✓ Alignement patients / labels vérifié\n")


# =============================================================================
# STEP 14 — DMP : POSITIONS DIFFÉRENTIELLEMENT MÉTHYLÉES (ChAMP)
# =============================================================================
# champ.DMP :
#   beta         : CpGs × samples (valeurs 0-1)
#   pheno        : vecteur de labels aligné sur les colonnes de beta
#   compare.group = c("A", "B") → résultat dans $<B>_to_<A>
#     logFC > 0 = hyperméthylé dans A par rapport à B
#   adjPVal = 1  → pas de pré-filtre interne, on applique nos propres seuils
#   Seuils retenus : adj.P.Val BH < 0.05 ET |logFC| > 0.3
#     (logFC ici = différence de M-values ≈ delta-beta ~ 0.10)

cat("\n========== STEP 14: DMP Méthylation (ChAMP) ==========\n")

dmp_gbm_vs_astro <- champ.DMP(
  beta          = methy_deg,
  pheno         = labels,
  compare.group = c("GBM", "ASTRO"),
  adjPVal       = 1,
  adjust.method = "BH",
  arraytype     = "450K"
)

dmp_gbm_vs_oligo <- champ.DMP(
  beta          = methy_deg,
  pheno         = labels,
  compare.group = c("GBM", "OLIGO"),
  adjPVal       = 1,
  adjust.method = "BH",
  arraytype     = "450K"
)

dmp_astro_vs_oligo <- champ.DMP(
  beta          = methy_deg,
  pheno         = labels,
  compare.group = c("ASTRO", "OLIGO"),
  adjPVal       = 1,
  adjust.method = "BH",
  arraytype     = "450K"
)

# ---- Annotation des résultats ----------------------------------------------
# Slot retourné : $<groupe2>_to_<groupe1>
#   compare.group = c("GBM", "ASTRO") → $ASTRO_to_GBM
#   logFC > 0 = hyperméthylé dans GBM (groupe1) vs ASTRO (groupe2)

annotate_dmp <- function(df, group1, group2) {
  df$comparison    <- paste0(group1, "_vs_", group2)
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
cat("  ✓ Résultats DMP sauvegardés\n")


# =============================================================================
# STEP 15 — DEG : GÈNES DIFFÉRENTIELLEMENT EXPRIMÉS (edgeR glmQLF)
# =============================================================================
# Protocole edgeR User Guide (GLM quasi-likelihood, section 4.1)
# Pas de voom ici : on modélise les counts directement (NB)
# Seuils retenus : FDR BH < 0.05 ET |logFC| > 1 (≥ 2x sur échelle log2)

cat("\n========== STEP 15: DEG mRNA (edgeR glmQLF) ==========\n")

# ---- 15a. Filtre low-count -------------------------------------------------
keep_genes <- rowSums(cpm(mrna_deg) > 1) >= 5
cat("  Gènes après filtre low-count :", sum(keep_genes), "/", nrow(mrna_deg), "\n")
mrna_filt <- mrna_deg[keep_genes, ]

# ---- 15b. DGEList + normalisation TMM --------------------------------------
dge <- DGEList(counts = mrna_filt)
dge <- normLibSizes(dge, method = "TMM")

# ---- 15c. Groupes depuis metadata_deg (source unique de vérité) ------------
dge$samples$group <- factor(labels, levels = c("GBM", "ASTRO", "OLIGO"))
cat("  Groupes dans DGEList :\n"); print(table(dge$samples$group))

# ---- 15d. Design matrix sans intercept -------------------------------------
design <- model.matrix(~ 0 + group, data = dge$samples)
colnames(design) <- levels(dge$samples$group)   # "GBM", "ASTRO", "OLIGO"

# ---- 15e. Dispersion -------------------------------------------------------
dge <- estimateDisp(dge, design, robust = TRUE)
cat("  Dispersion commune :", round(dge$common.dispersion, 4), "\n")

# ---- 15f. Fit quasi-likelihood négative binomiale --------------------------
fit <- glmQLFit(dge, design, robust = TRUE)

# ---- 15g. Contrastes (logFC > 0 = surexprimé dans groupe1) -----------------
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
cat("  ✓ Résultats DEG sauvegardés\n")


# =============================================================================
# STEP 16 — RÉSUMÉ FINAL
# =============================================================================
cat("\n========== STEP 16: RÉSUMÉ FINAL ==========\n")

cat("\n--- Méthylation (DMP) : adj.P.Val < 0.05, |logFC| > 0.3 ---\n")
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

cat("\n--- mRNA (DEG) : FDR < 0.05, |logFC| > 1 ---\n")
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

cat("\n=== PIPELINE DEG COMPLET — TERMINÉ ===\n")

# =============================================================================
# FICHIERS PRODUITS
# =============================================================================
# DataSets/processed_assays/
#   assay_rna_coding.csv          mRNA coding  log2-CPM z-scoré  (pour MOFA)
#   assay_rna_lnc.csv             lncRNA       log2-CPM z-scoré  (pour MOFA)
#   assay_rna_mirna.csv           miRNA        log2-CPM z-scoré  (pour MOFA)
#   assay_methylation.csv         Méthylation  M-values z-scorées (pour MOFA)
#   raw_counts_coding_aligned.csv Raw counts   alignés metadata  (pour DEG)
#   beta_values_aligned.csv       Beta values  alignées metadata  (pour DEG)
#   metadata_final.csv            Metadata 318 patients
#
# dge-results/
#   gbm_vs_astro_methy.csv        DMP GBM vs ASTRO
#   gbm_vs_oligo_methy.csv        DMP GBM vs OLIGO
#   astro_vs_oligo_methy.csv      DMP ASTRO vs OLIGO
#   gbm_vs_astro_mrna.csv         DEG GBM vs ASTRO
#   gbm_vs_oligo_mrna.csv         DEG GBM vs OLIGO
#   astro_vs_oligo_mrna.csv       DEG ASTRO vs OLIGO
# =============================================================================
