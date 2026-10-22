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

#* [1] FILTER breakdown of the RAW vcf (counted by 01.QC_01.Normalise.py)
print(fread("raw_filter_breakdown.txt", sep = "\t"))

#* [2] parse bcftools stats of the analysable set
st <- readLines("proband.norm.pass.stats.txt")
sn <- grep("^SN", st, value = TRUE)
writeLines(sn)                                   # records / SNPs / indels / MNPs
titv <- grep("^TSTV", st, value = TRUE); writeLines(titv)   # expect ~2.0-2.1 genome-wide

#* [3] genotype composition from the allele-balance table
ab <- fread("proband.ab_table.txt", colClasses = list(character = "CHROM"))
ab <- ab[CHROM %in% c(as.character(1:22), "X", "Y")]
ab[, zyg := fifelse(GT %in% c("1/1", "1|1"), "hom_alt",
             fifelse(GT %in% c("0/1", "1/0", "0|1", "1|0"), "het", "other"))]
print(ab[, .N, by = zyg])
cat("het/hom-alt ratio: ", ab[zyg == "het", .N] / ab[zyg == "hom_alt", .N], "\n")  # expect ~1.5-2.0


#* [4] sex check from X / Y
xhet <- ab[CHROM == "X" & DP >= 10, mean(zyg == "het")]
ycnt <- ab[CHROM == "Y" & DP >= 10, .N]
cat(sprintf("X het fraction = %.3f ; Y calls = %d  -> inferred sex: %s\n",
            xhet, ycnt, ifelse(xhet < 0.15 & ycnt > 50, "male (XY)", "female (XX) / uncertain")))

fwrite(ab[, .N, by = .(CHROM, zyg)], file.path(source_dir2, "proband.zyg_by_chrom.txt"), sep = "\t")
