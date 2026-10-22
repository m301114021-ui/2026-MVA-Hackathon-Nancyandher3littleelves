# Stage 01 step 1 — normalise + hard-filter the proband VCF, write the genotype /
# allele-balance table and bcftools stats. Replaces the bash block of the old
# Main_MAV_Hackathon_01.QC.sh (same filters, same outputs).
import subprocess
import collections
from pathlib import Path

# working directory
working_dir = Path("/home/ianchen/Documents/Projects/MAV_Hackathon/MVA")
# source directory
source_dir  = working_dir / "data/raw"
source_dir2 = working_dir / "data/01_interim"

Input    = source_dir / "WGS_EX2312012_HGWCNDSX7.vcf.gz"
Fasta    = "/home/ianchen/.vep/homo_sapiens/105_GRCh38/Homo_sapiens.GRCh38.dna.toplevel.fa"
bcftools = "/home/ianchen/.conda/envs/bcftools/bin/bcftools"

Norm_VCF = source_dir2 / "proband.norm.pass.vcf.gz"


def run(cmd):
    print("+", cmd, flush=True)
    subprocess.run(cmd, shell=True, check=True, executable="/bin/bash")


#* [1] FILTER breakdown of the RAW vcf (feeds Fig.01_filter_breakdown in 01.QC_02)
res = subprocess.run([bcftools, "query", "-f", "%FILTER\n", str(Input)],
                     capture_output=True, text=True, check=True)
filt = collections.Counter(res.stdout.splitlines())
with open(source_dir2 / "raw_filter_breakdown.txt", "w") as f:
    f.write("count\tfilter\n")
    for k, n in filt.most_common():
        f.write(f"{n}\t{k}\n")
print("raw records: ", sum(filt.values()), " PASS: ", filt["PASS"])

#* [2] normalise + hard-filter to analysable proband variants
# - view -f PASS -r 1..22,X,Y   PASS sites on the standard chromosomes only
#   (drops chrM + decoys; the toplevel FASTA names MT not M -> -c w just warns)
# - norm -m -any -f <fasta>     split multiallelics, left-align, check REF
# - view -e 'GT="0/0"|...'      drop homozygous-ref / missing in the proband
# - view -i FMT/DP>=8 & FMT/GQ>=20   genotype-level quality gate
# - annotate --set-id           STABLE JOIN KEY: ID = CHROM_POS_REF_ALT, reused
#   by VEP (Uploaded_variation) and every downstream merge
regions = ",".join([str(i) for i in range(1, 23)] + ["X", "Y"])
run(f"set -o pipefail; "
    f"{bcftools} view -f PASS -r {regions} {Input} -Ou "
    f"| {bcftools} norm -m -any -f {Fasta} -c w -Ou 2> {source_dir2}/norm.err "
    f"| {bcftools} view -e 'GT=\"0/0\" | GT=\"./.\" | GT=\".\"' -Ou "
    f"| {bcftools} view -i 'FMT/DP>=8 & FMT/GQ>=20' -Ou "
    f"| {bcftools} annotate --set-id '%CHROM\\_%POS\\_%REF\\_%ALT' -Oz -o {Norm_VCF}")
run(f"{bcftools} index -f -t {Norm_VCF}")
run(f"{bcftools} stats {Norm_VCF} > {source_dir2}/proband.norm.pass.stats.txt")

#* [3] a lightweight allele-balance table for the mosaic-aneuploidy scan
# CHROM POS GT AD_REF AD_ALT DP GQ  -> feeds 01.QC_02/03 and stage 02/03
#! split AD into two columns — a bare "ref,alt" string is mis-parsed as a decimal
#! number by data.table::fread under a comma-decimal locale
run(f"{{ printf 'MarkerID\\tCHROM\\tPOS\\tREF\\tALT\\tGT\\tAD_REF\\tAD_ALT\\tDP\\tGQ\\n'; "
    f"{bcftools} query -f '%ID\\t%CHROM\\t%POS\\t%REF\\t%ALT\\t[%GT\\t%AD{{0}}\\t%AD{{1}}\\t%DP\\t%GQ]\\n' {Norm_VCF}; }} "
    f"> {source_dir2}/proband.ab_table.txt")

n = sum(1 for _ in open(source_dir2 / "proband.ab_table.txt")) - 1
print("analysable proband variants: ", n)
# expected from the Aug-2026 run: 4,682,895
