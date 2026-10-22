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

setwd(source_dir2)
library(data.table)
library(dplyr)
library(stringr)

#* [1] read the VEP --tab output
# real header is the last "##"-prefixed line, starting "#Uploaded_variation"
hdr <- readLines("proband.vep.txt", n = 300)
hdr_line <- grep("^#Uploaded_variation", hdr, value = TRUE)
vep_cols <- strsplit(sub("^#", "", hdr_line), "\t")[[1]]

vep <- fread("proband.vep.txt", skip = "#Uploaded_variation", header = FALSE,
             sep = "\t", na.strings = c("-", "", "NA"),
             colClasses = "character")
setnames(vep, vep_cols)
cat("VEP rows (1 per variant after --pick): ", nrow(vep), "\n")

#* [2] MarkerID = Uploaded_variation (we set VCF ID = CHROM_POS_REF_ALT in stage 01)
vep[, MarkerID := Uploaded_variation]
vep <- vep[str_detect(MarkerID, "^[0-9XY]+_[0-9]+_")]        # drop any odd rows

num <- function(x) suppressWarnings(as.numeric(x))
vep[, `:=`(gnomAD_AF = num(gnomAD_AF), MAX_AF = num(MAX_AF))]
pop_cols <- grep("^gnomAD_.*_AF$", names(vep), value = TRUE)
if (length(pop_cols)) {
  vep[, (pop_cols) := lapply(.SD, num), .SDcols = pop_cols]
  vep[, gnomAD_popmax := do.call(pmax, c(.SD, list(na.rm = TRUE))), .SDcols = pop_cols]
  vep[is.infinite(gnomAD_popmax), gnomAD_popmax := NA_real_]
} else {
  vep[, gnomAD_popmax := gnomAD_AF]
}

#* [3] consequence severity rank (Sequence Ontology, high -> low)
severe <- c("transcript_ablation","splice_acceptor_variant","splice_donor_variant",
            "stop_gained","frameshift_variant","stop_lost","start_lost",
            "transcript_amplification","inframe_insertion","inframe_deletion",
            "missense_variant","protein_altering_variant","splice_region_variant",
            "splice_donor_5th_base_variant","splice_donor_region_variant",
            "splice_polypyrimidine_tract_variant","incomplete_terminal_codon_variant",
            "start_retained_variant","stop_retained_variant","synonymous_variant")
vep[, csq1 := tstrsplit(Consequence, ",", fixed = TRUE)[[1]]]
vep[, csq_rank := match(csq1, severe)]
vep[is.na(csq_rank), csq_rank := 999L]

#* [4] ClinVar / predictor flags
# mixed assertions (e.g. "benign,pathogenic") are NOT counted as pathogenic
vep[, clinvar_patho := !is.na(CLIN_SIG) &
      str_detect(tolower(CLIN_SIG), "pathogenic") &
      !str_detect(tolower(CLIN_SIG), "conflicting|benign")]
vep[, sift_del    := !is.na(SIFT) & str_detect(SIFT, "deleterious")]
vep[, polyphen_dm := !is.na(PolyPhen) & str_detect(PolyPhen, "damaging")]

#* [5] attach the proband genotype (AB table from stage 01, keyed by MarkerID)
gt <- fread(file.path(source_dir2, "proband.ab_table.txt"),
            colClasses = list(character = c("MarkerID", "CHROM", "REF", "ALT")))
gt[, `:=`(AD_REF = suppressWarnings(as.integer(AD_REF)),
          AD_ALT = suppressWarnings(as.integer(AD_ALT)))]
gt[, zyg := fifelse(GT %in% c("1/1","1|1"), "hom_alt",
            fifelse(GT %in% c("0/1","1/0","0|1","1|0"), "het", "other"))]
gt[, ab := fifelse(!is.na(AD_REF) & (AD_REF + AD_ALT) > 0,
                   AD_ALT / (AD_REF + AD_ALT), NA_real_)]

ann <- merge(vep, gt[, .(MarkerID, CHROM, POS, REF, ALT, GT, zyg, AD_REF, AD_ALT, DP, GQ, ab)],
             by = "MarkerID", all.x = TRUE)
cat("annotated variants with a genotype: ", ann[!is.na(GT), .N], "\n")

keep <- c("MarkerID","CHROM","POS","REF","ALT","REF_ALLELE","Allele","SYMBOL","Gene",
          "Feature","Consequence","csq1","csq_rank","IMPACT","BIOTYPE","CANONICAL",
          "EXON","INTRON","HGVSc","HGVSp","Existing_variation","CLIN_SIG","clinvar_patho",
          "SIFT","sift_del","PolyPhen","polyphen_dm",
          "gnomAD_AF","gnomAD_popmax","MAX_AF","MAX_AF_POPS",
          "GT","zyg","AD_REF","AD_ALT","DP","GQ","ab")
ann <- ann[, intersect(keep, names(ann)), with = FALSE]
setorder(ann, csq_rank, gnomAD_popmax, na.last = TRUE)

# gzip on write — the full annotation is ~1 GB uncompressed and only stage 03_01
# reads it; data.table::fread reads .gz transparently
fwrite(ann, file.path(source_dir3, "proband.annot.txt.gz"), sep = "\t", compress = "gzip")
cat("wrote data/02_processed/proband.annot.txt.gz : ", nrow(ann), " rows\n")
print(ann[, .N, by = IMPACT])
