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
library(dplyr)

rc <- fread("proband.rare_coding.analysis.txt")    # call QC (03_02) + gnomAD v4 (03_03)

#* [1] autosomal recessive — homozygous alt
homo <- rc[zyg == "hom_alt" & CHROM %in% as.character(1:22)]
# a real recessive hit is essentially absent as homozygote in gnomAD; AF filter
# from stage 01 already enforces rarity. Prefer GQ>=30, DP>=10.
homo <- homo[GQ >= 30 & DP >= 10]
setorder(homo, csq_rank)
fwrite(homo, "proband.model_homozygous.txt", sep = "\t")
cat("homozygous-alt candidates: ", nrow(homo), " in ", uniqueN(homo$SYMBOL), " genes\n")

#* [2] autosomal recessive — potential compound heterozygous
# >=2 rare het impactful variants in the SAME gene (phase unknown, trio-less)
het <- rc[zyg == "het" & CHROM %in% as.character(1:22) & GQ >= 30 & DP >= 10]
gene_n <- het[!is.na(SYMBOL), .N, by = SYMBOL][N >= 2]
comphet <- het[SYMBOL %in% gene_n$SYMBOL]
setorder(comphet, SYMBOL, csq_rank)
fwrite(comphet, "proband.model_comphet.txt", sep = "\t")
cat("potential compound-het: ", nrow(comphet), " variants in ", nrow(gene_n), " genes\n")
print(gene_n[order(-N)])

#* [3] X-linked — hemizygous / homozygous on X (relevant if proband is male)
xsex <- readLines(file.path(source_dir2, "proband.zyg_by_chrom.txt"))  # informational
xl <- rc[CHROM == "X" & GQ >= 30 & DP >= 10]
# on a male X almost everything is called hom_alt("1/1"); keep both, flag zyg
setorder(xl, csq_rank)
fwrite(xl, "proband.model_xlinked.txt", sep = "\t")
cat("X-linked candidates: ", nrow(xl), " in ", uniqueN(xl$SYMBOL), " genes\n")

#* [4] candidate de novo (dominant) — trio-less proxy
# absent from gnomAD entirely (v2.1.1 exome in the VEP cache AND v4.1 from 03_03),
# heterozygous, very high genotype quality, balanced allele fraction (not an artefact)
#! without parents this is a PRIVATE heterozygous variant, not a confirmed de novo
denovo <- rc[zyg == "het" &
             (is.na(gnomAD_AF) | gnomAD_AF == 0) &
             (is.na(MAX_AF)    | MAX_AF == 0) &
             gnomad4_found == FALSE &
             GQ >= 50 & DP >= 15 &
             !is.na(ab) & ab >= 0.3 & ab <= 0.7]
setorder(denovo, csq_rank)
fwrite(denovo, "proband.model_denovo_candidate.txt", sep = "\t")
cat("candidate de novo (needs parental confirmation): ", nrow(denovo),
    " in ", uniqueN(denovo$SYMBOL), " genes\n")

#! Summary counts to paste back here after the run:
#!   homozygous  = ____   comphet genes = ____   X-linked = ____   de novo cand = ____
