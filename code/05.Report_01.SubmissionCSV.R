# working directory
working_dir <- file.path("/home/ianchen/Documents/Projects/MAV_Hackathon/MVA")
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

#! Track 1 submission settings
team       <- "Nancyandher3littleelves"
approach   <- "vep105-gnomad4-panelrank"
proband_id <- "WGS_EX2312012"
N_ROWS     <- 10                          # hackathon accepts up to 10 rows

ranked <- fread("proband.ranked_candidates.txt", na.strings = c("NA", ""),
                colClasses = list(character = c("MarkerID", "CHROM")))

#* [1] VCF alleles from the stable MarkerID (CHROM_POS_REF_ALT, bcftools-normalised GRCh38)
ranked[, c("v_chrom", "v_pos", "v_ref", "v_alt") := tstrsplit(MarkerID, "_", fixed = TRUE)]

#* [2] drop highly polymorphic / paralogous gene families
#! operational choice, NOT an established standard: these families fill rare-variant
#! lists through mapping artefacts in segmental duplications. Only affects the
#! secondary rows here (the top unit outscores them either way, see stage 04).
noisy <- "^(HLA-|KIR[23]D|LILR[AB]|MUC[0-9]|PRSS[0-9]|TAS2R|OR[0-9]+[A-Z]|ZNF[0-9]|FRG2|PRAMEF|NBPF|GOLGA6L|GOLGA8|KRTAP|TPTE|PSG[0-9]|USP17L|FAM90A|DEFB)"
ranked[, noisy_family := !is.na(SYMBOL) & grepl(noisy, SYMBOL)]
cat("dropped noisy-family variants: ", ranked[noisy_family == TRUE, .N],
    " in ", ranked[noisy_family == TRUE, uniqueN(SYMBOL)], " genes\n")
cand <- ranked[noisy_family == FALSE & !is.na(SYMBOL)]
setorder(cand, -score, POS)

#* [3] candidate units
# recessive gene with >=2 potential comp-het variants -> ONE pair row (top two by
# score, ordered by position); unit score = the weaker of the two, both alleles are needed
ch_genes <- cand[is_comphet == TRUE, .N, by = SYMBOL][N >= 2, SYMBOL]
pairs <- cand[is_comphet == TRUE & SYMBOL %in% ch_genes, head(.SD, 2), by = SYMBOL]
setorder(pairs, SYMBOL, POS)
pair_units <- pairs[, .(score = min(score), id1 = MarkerID[1], id2 = MarkerID[2],
                        model = "compound_het",
                        csq_rank = min(csq_rank), af = max(af_use, na.rm = TRUE)), by = SYMBOL]
pair_units[is.infinite(af), af := NA_real_]

# every other variant is its own unit
single_units <- cand[!MarkerID %in% c(pair_units$id1, pair_units$id2),
  .(SYMBOL, score, id1 = MarkerID, id2 = NA_character_,
    model = fifelse(is_homozygous == TRUE, "homozygous",
            fifelse(is_xlinked == TRUE, "X-linked",
            fifelse(is_denovo_cand == TRUE, "private_het", "heterozygous"))),
    csq_rank, af = af_use)]

units <- rbind(pair_units, single_units)
# ties: more severe consequence, then rarer (absent = rarest), then gene name
units[, af_sort := fifelse(is.na(af), 0, af)]
setorder(units, -score, csq_rank, af_sort, SYMBOL)

# one row per gene (its best unit); genes on a curated HPO panel rank ahead of
# off-panel genes, then by score and the tie-breakers above
#! operational choice, not an established standard
units <- units[, .SD[1], by = SYMBOL]
units[, on_panel := SYMBOL %in% cand[!is.na(panel), SYMBOL]]
setorder(units, -on_panel, -score, csq_rank, af_sort, SYMBOL)

#* [4] EPCR = softmax of the unit score over ALL candidate units
#! a transform of the stage-04 score, NOT a calibrated probability; only the
#! ranking is scored by the hackathon. Made non-increasing down the rows so the
#! EPCR order always matches the row order (panel-first ordering above).
units[, epcr := exp(score - max(score))]
units[, epcr := signif(cummin(epcr / sum(epcr)), 3)]

top <- units[1:min(N_ROWS, .N)]
top[, finding_type := fifelse(.I == 1L, "primary", "secondary")]

#* [5] human-readable notes
v <- cand[, .(MarkerID, v_chrom, v_pos, v_ref, v_alt, Consequence, zyg, CLIN_SIG, panel,
              prot = sub("^[^:]+:", "", fifelse(is.na(HGVSp) | HGVSp == "", HGVSc, HGVSp)),
              af = fifelse(is.na(af_use), "absent", formatC(af_use, format = "e", digits = 2)))]
setkey(v, MarkerID)
vdesc <- function(id) {
  r <- v[J(id)]
  sprintf("%s %s %s AF=%s ClinVar=%s", r$prot, r$Consequence, r$zyg, r$af,
          fifelse(is.na(r$CLIN_SIG) | r$CLIN_SIG == "", "none", r$CLIN_SIG))
}
top[, panel := v[J(id1), panel]]
top[, notes := sprintf("%s [%s]: %s%s; panel=%s; score=%s",
                       SYMBOL, model, vdesc(id1),
                       fifelse(is.na(id2), "", paste0(" + ", vdesc(id2))),
                       fifelse(is.na(panel), "none", panel), score)]

#* [6] assemble the submission table (GRCh38, "chr" prefix as in the template)
a1 <- v[J(top$id1)]
a2 <- v[J(top$id2)]
sub_dt <- data.table(
  proband_id = proband_id,
  chrom_1 = paste0("chr", a1$v_chrom), pos_1 = as.integer(a1$v_pos),
  ref_1 = a1$v_ref, alt_1 = a1$v_alt,
  chrom_2 = fifelse(is.na(a2$v_chrom), NA_character_, paste0("chr", a2$v_chrom)),
  pos_2 = as.integer(a2$v_pos), ref_2 = a2$v_ref, alt_2 = a2$v_alt,
  epcr = top$epcr, finding_type = top$finding_type, notes = top$notes)

# naming rule: <team>_<approach>.csv
out_csv <- file.path(output_dir, sprintf("%s_%s.csv", team, approach))
fwrite(sub_dt, out_csv, na = "", scipen = 50)
fwrite(units, file.path(source_dir3, "proband.submission_units.txt"), sep = "\t")
cat("wrote ", out_csv, "\n")
print(sub_dt[, .(chrom_1, pos_1, ref_1, alt_1, chrom_2, pos_2, epcr, finding_type)])
