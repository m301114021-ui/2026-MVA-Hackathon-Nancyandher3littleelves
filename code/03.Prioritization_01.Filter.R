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

ann <- fread("proband.annot.txt.gz", na.strings = c("NA", ""))
cat("input annotated variants: ", nrow(ann), "\n")

#* [1] frequency filter
# MAX_AF is the max over 1000G + ESP + gnomAD exomes; treat missing as ultra-rare
AF_CUT <- 1e-3
ann[, af_use := pmax(MAX_AF, gnomAD_popmax, gnomAD_AF, na.rm = TRUE)]
ann[is.infinite(af_use), af_use := NA_real_]
rare <- ann[is.na(af_use) | af_use <= AF_CUT]
cat("rare (MAX_AF <= ", AF_CUT, " or absent): ", nrow(rare), "\n")

#* [2] consequence filter
# keep protein-truncating / missense / inframe / splice-region, OR ClinVar patho
impactful_csq <- c("transcript_ablation","splice_acceptor_variant","splice_donor_variant",
                   "stop_gained","frameshift_variant","stop_lost","start_lost",
                   "inframe_insertion","inframe_deletion","missense_variant",
                   "protein_altering_variant","splice_region_variant",
                   "splice_donor_5th_base_variant","splice_donor_region_variant",
                   "splice_polypyrimidine_tract_variant")
rare[, impactful := csq1 %in% impactful_csq | IMPACT %in% c("HIGH","MODERATE")]
rare_coding <- rare[impactful == TRUE | clinvar_patho == TRUE]
cat("rare + impactful (or ClinVar pathogenic): ", nrow(rare_coding), "\n")

#* [3] light protein-coding biotype restriction, keep the genotype
rare_coding <- rare_coding[is.na(BIOTYPE) | BIOTYPE == "protein_coding"]
rare_coding <- rare_coding[!is.na(GT) & zyg %in% c("het","hom_alt")]
cat("with a real proband genotype: ", nrow(rare_coding), "\n")

#* [4] a simple in-silico damage tier for missense
rare_coding[, missense_support := fifelse(csq1 == "missense_variant",
              (as.integer(sift_del) + as.integer(polyphen_dm)), NA_integer_)]

setorder(rare_coding, csq_rank, af_use, na.last = TRUE)
fwrite(rare_coding, "proband.rare_coding.txt", sep = "\t")
cat("wrote proband.rare_coding.txt\n")
print(rare_coding[, .N, by = csq1][order(-N)])
