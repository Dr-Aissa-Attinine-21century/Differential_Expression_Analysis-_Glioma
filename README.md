# DEG Pipeline — TCGA Glioma (GBM + LGG)

## English

This repository contains a complete end-to-end pipeline for differential expression and methylation analysis of TCGA glioma data (GBM and LGG), restricted to the 318-patient MOFA cohort. The two R scripts are identical in logic and differ only in language of comments and messages.

| File | Description |
|------|-------------|
| `deg_pipeline_final_EN.R` | Full pipeline in **English** |
| `deg_pipeline_final.R` | Full pipeline in **French** |

**Part A — Data preparation** downloads raw data from TCGA via `TCGAbiolinks`, deduplicates barcodes to patient-level submitter IDs, filters strictly to the metadata cohort, removes missing values and outliers, splits genes by type (coding / lncRNA / miRNA), applies TMM + voom normalization for mRNA and beta-to-M-value transformation for methylation, selects features by variance, and z-scores all assays before saving. **Part B — Differential analysis** restarts from the aligned raw data (counts + beta values) and runs `edgeR glmQLFTest` for DEG (FDR < 0.05, |logFC| > 1) and `champ.DMP` for DMP (adj.P.Val < 0.05, |logFC| > 0.3) across three pairwise comparisons: GBM vs ASTRO, GBM vs OLIGO, and ASTRO vs OLIGO.

---

## Français

Ce dépôt contient un pipeline complet pour l'analyse différentielle d'expression et de méthylation sur les données TCGA gliome (GBM et LGG), restreinte à la cohorte MOFA de 318 patients. Les deux scripts R sont identiques dans leur logique et ne diffèrent que par la langue des commentaires et des messages.

| Fichier | Description |
|---------|-------------|
| `deg_pipeline_final_EN.R` | Pipeline complet en **anglais** |
| `deg_pipeline_final.R` | Pipeline complet en **français** |

**Partie A — Préparation des données** télécharge les données brutes TCGA via `TCGAbiolinks`, déduplique les barcodes en identifiants patient, filtre strictement sur la cohorte metadata, supprime les valeurs manquantes et les outliers, sépare les gènes par type (coding / lncRNA / miRNA), applique la normalisation TMM + voom pour le mRNA et la transformation beta → M-values pour la méthylation, sélectionne les features par variance et applique un z-score avant sauvegarde. **Partie B — Analyse différentielle** repart des données brutes alignées (counts + beta values) et exécute `edgeR glmQLFTest` pour le DEG (FDR < 0.05, |logFC| > 1) et `champ.DMP` pour les DMP (adj.P.Val < 0.05, |logFC| > 0.3) sur trois comparaisons : GBM vs ASTRO, GBM vs OLIGO, et ASTRO vs OLIGO.
