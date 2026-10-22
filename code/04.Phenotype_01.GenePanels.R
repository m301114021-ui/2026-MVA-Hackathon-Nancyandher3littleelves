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

#! Curated offline gene panels. Sources: OMIM / GeneReviews / Genomics England
#! PanelApp ("Mosaic variegated aneuploidy", "Childhood solid tumours",
#! "Primordial dwarfism", "Nephrocalcinosis or nephrolithiasis"). Hand-maintained
#! here so the project is self-contained; refresh against PanelApp before writing up.

panels <- list(

  # tier-1: the direct disease hypothesis
  MVA_syndrome = c("BUB1B","CEP57","TRIP13"),

  # mitotic spindle-assembly checkpoint / chromosome segregation / centrosome
  chromosomal_instability = c(
    "BUB1","BUB3","MAD1L1","MAD2L1","TTK","CENPE","CENPF","KNL1","CASC5",
    "AURKB","AURKC","PLK4","CEP192","SPICE1","ZW10","ZWILCH","SPINDLY","KNTC1",
    "NUF2","NDC80","SPC24","SPC25","CDC20","ESPL1","PTTG1"),

  # childhood cancer predisposition (rhabdomyosarcoma / soft-tissue / embryonal)
  cancer_predisposition = c(
    "TP53","NF1","DICER1","SMARCB1","SMARCA4","PTCH1","SUFU","HRAS","FBXW7","CBL",
    "PTPN11","SOS1","RASA2","NBN","BLM","RECQL4","WRN","RECQL",
    "FANCA","FANCB","FANCC","FANCD2","FANCE","FANCF","FANCG","FANCI","FANCL","FANCM",
    "BRCA1","BRCA2","PALB2","BRIP1","RAD51C","RAD51D","SLX4","XRCC2","ERCC4",
    "MSH2","MLH1","MSH6","PMS2","EPCAM","POLE","POLD1"),

  # microcephalic primordial dwarfism / Meier-Gorlin / severe IUGR + DNA replication/repair
  primordial_dwarfism = c(
    "PCNT","ATR","ATRIP","RBBP8","CENPJ","CEP152","CEP63","CDK5RAP2","ASPM","WDR62",
    "DNA2","MCM5","ORC1","ORC4","ORC6","CDT1","CDC6","GMNN","CDC45","GINS1",
    "LIG4","XRCC4","NHEJ1","NSMCE2","BRCA1","BLM","POLE","POLA1","DONSON","MCM3AP",
    "IGF1","IGF1R","IGF2","GH1","OBSL1","CUL7","CCDC8"),

  # nephrocalcinosis / nephrolithiasis / renal tubular
  nephrocalcinosis = c(
    "SLC34A1","SLC34A3","SLC9A3R1","VDR","CYP24A1","CYP2R1","CASR","CLDN16","CLDN19",
    "CLCN5","OCRL","SLC4A1","ATP6V0A4","ATP6V1B1","SLC12A1","KCNJ1","BSND","CLCNKB",
    "MAGED2","FAM20A","AGXT","GRHPR","HOGA1","APRT","SLC7A9","SLC3A1","ADCY10",
    "SLC26A1","HNF4A","CTNS","CLDN14"),

  # short stature (non-dwarfism syndromic overlap)
  short_stature = c(
    "SHOX","NPR2","ACAN","FGFR3","PTPN11","NF1","IGF1R","STAT5B","GHR","GH1",
    "IHH","COL2A1","COL10A1","FBN1","BMP2","NPPC"))

pan <- rbindlist(lapply(names(panels), function(p)
  data.table(SYMBOL = panels[[p]], panel = p)))
pan <- unique(pan)

# panel weights (evidence strength for THIS phenotype)
pw <- c(MVA_syndrome = 5, chromosomal_instability = 3, cancer_predisposition = 3,
        primordial_dwarfism = 2, nephrocalcinosis = 2, short_stature = 1)
pan[, panel_weight := pw[panel]]

fwrite(pan, "gene_panels.txt", sep = "\t")
cat("panel genes: ", uniqueN(pan$SYMBOL), " across ", uniqueN(pan$panel), " panels\n")

#* gene -> which proband features it explains (for a phenotype-overlap score)
# features: CA=cancer/RMS, NC=nephrocalcinosis, GR=growth restriction/short stature/SGA,
#           FTT=failure to thrive/muscle, PERI=prematurity, MISC=parental miscarriage/aneuploidy
feat <- fread(text = 'SYMBOL	features
BUB1B	CA;NC;GR;FTT;PERI;MISC
CEP57	GR;MISC;CA
TRIP13	CA;GR;MISC
TP53	CA
DICER1	CA
SMARCB1	CA
NF1	CA;GR
PCNT	GR;PERI;MISC;CA
CEP152	GR;PERI;MISC
CENPJ	GR;PERI
ATR	GR;PERI;CA
BLM	CA;GR;PERI;MISC
NBN	CA;GR;MISC
RECQL4	CA;GR;NC
FANCA	CA;GR;PERI;MISC
FANCD2	CA;GR;MISC
BRCA2	CA;GR;MISC
LIG4	GR;PERI;CA;MISC
SLC34A1	NC;GR
SLC34A3	NC;GR
CYP24A1	NC
CASR	NC
CLDN16	NC;GR
CLDN19	NC
HOGA1	NC
AGXT	NC
GRHPR	NC
SHOX	GR
ACAN	GR
NPR2	GR
IGF1R	GR;PERI', header = TRUE)

fwrite(feat, "gene_feature_map.txt", sep = "\t")
cat("gene-feature map rows: ", nrow(feat), "\n")
