# working directory
working_dir <- file.path("/home/ianchen/Documents/Projects/MVA_Hackathon/MVA")
# source directory
source_dir  <- file.path(working_dir, "data/raw")
source_dir2 <- file.path(working_dir, "data/01_interim")
source_dir3 <- file.path(working_dir, "data/02_processed")

# set output directory
output_dir  <- file.path(working_dir, "out")

cat("working_dir: ", working_dir)
cat("source_dir: ",  source_dir)
cat("output_dir: ",  output_dir)

setwd(source_dir3)
library(data.table)
library(dplyr)

rc   <- fread("proband.rare_coding.analysis.txt")  # call QC (03_02) + gnomAD v4 (03_03)
pan  <- fread("gene_panels.txt")
feat <- fread("gene_feature_map.txt")

#* [1] gene-level panel annotation (keep the strongest panel per gene)
pan_best <- pan[order(-panel_weight)][, .SD[1], by = SYMBOL]
cand <- merge(rc, pan_best, by = "SYMBOL", all.x = TRUE)

#* [2] flag inheritance model membership (from stage 03 outputs)
read_model <- function(f) if (file.exists(f)) unique(fread(f)$MarkerID) else character(0)
m_homo   <- read_model("proband.model_homozygous.txt")
m_ch     <- read_model("proband.model_comphet.txt")
m_x      <- read_model("proband.model_xlinked.txt")
m_dn     <- read_model("proband.model_denovo_candidate.txt")
cand[, `:=`(
  is_homozygous = MarkerID %in% m_homo,
  is_comphet    = MarkerID %in% m_ch,
  is_xlinked    = MarkerID %in% m_x,
  is_denovo_cand= MarkerID %in% m_dn)]

#* [3] phenotype-overlap score: how many proband feature-classes the gene explains
feat[, n_feat := lengths(strsplit(features, ";"))]
cand <- merge(cand, feat[, .(SYMBOL, gene_features = features, n_feat)],
              by = "SYMBOL", all.x = TRUE)
cand[is.na(n_feat), n_feat := 0L]

#* [4] scoring
# consequence
csq_score <- c(
  transcript_ablation = 6, splice_acceptor_variant = 6, splice_donor_variant = 6,
  stop_gained = 6, frameshift_variant = 6, start_lost = 5, stop_lost = 5,
  inframe_deletion = 3, inframe_insertion = 3, protein_altering_variant = 3,
  missense_variant = 2, splice_region_variant = 2)
cand[, s_csq := csq_score[csq1]]
cand[is.na(s_csq), s_csq := 1]
# missense in-silico bump
cand[csq1 == "missense_variant" & !is.na(missense_support), s_csq := s_csq + missense_support]

# rarity (absent from gnomAD is the strongest)
cand[, s_rare := fifelse(is.na(af_use), 3,
               fifelse(af_use <= 1e-5, 2,
               fifelse(af_use <= 1e-4, 1, 0)))]

# ClinVar
cand[, s_clin := fifelse(clinvar_patho == TRUE, 4, 0)]

# panel
cand[, s_panel := fifelse(is.na(panel_weight), 0, panel_weight)]

# inheritance model plausibility (recessive strongly favoured by parental miscarriage history)
cand[, s_model := 0]
cand[is_homozygous == TRUE, s_model := s_model + 4]
cand[is_comphet    == TRUE, s_model := s_model + 3]
cand[is_xlinked    == TRUE, s_model := s_model + 2]
cand[is_denovo_cand== TRUE, s_model := s_model + 2]

# phenotype overlap
cand[, s_pheno := pmin(n_feat, 4)]

cand[, score := s_csq + s_rare + s_clin + s_panel + s_model + s_pheno]

#* [5] rank; genes on a panel OR carrying a strong model float to the top
setorder(cand, -score, csq_rank, af_use, na.last = TRUE)
cols <- c("score","SYMBOL","Gene","MarkerID","CHROM","POS","REF_ALLELE","Allele",
          "Consequence","IMPACT","HGVSc","HGVSp","zyg","GT","DP","GQ",
          "csq_rank","af_use","gnomAD_AF","gnomad4_af","gnomad4_faf95_popmax","gnomad4_hom",
          "CLIN_SIG","SIFT","PolyPhen",
          "panel","panel_weight","gene_features",
          "is_homozygous","is_comphet","is_xlinked","is_denovo_cand",
          "s_csq","s_rare","s_clin","s_panel","s_model","s_pheno")
cols <- intersect(cols, names(cand))
ranked <- cand[, ..cols]

fwrite(ranked, "proband.ranked_candidates.txt", sep = "\t")
cat("ranked candidates: ", nrow(ranked), "\n")
cat("\n--- top 25 ---\n")
print(ranked[1:min(25, .N)])
cat("\n--- on a curated panel ---\n")
print(ranked[!is.na(panel)][1:min(30, .N),
      .(score, SYMBOL, Consequence, zyg, HGVSp, af_use, panel, gene_features)])
