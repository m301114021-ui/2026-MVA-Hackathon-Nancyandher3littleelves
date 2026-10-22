# Stage 02 step 1 — Ensembl VEP (release 105, GRCh38, offline cache) on the
# analysable proband set (~4.7 M variants). Same flags as the Aug-2026 run.
#! VEP cache 105 carries gnomAD r2.1.1 EXOME AF (--af_gnomad / --max_af),
#! ClinVar 2021-06, dbSNP 154. No CADD / dbNSFP / SpliceAI plugin -> we rely on
#! SIFT, PolyPhen, consequence severity and ClinVar. gnomAD GENOME AF is NOT in
#! this cache; re-check the shortlist on gnomad.broadinstitute.org.
import os
import subprocess
from pathlib import Path

# working directory
working_dir = Path("/home/ianchen/Documents/Projects/MVA_Hackathon/MVA")
source_dir2 = working_dir / "data/01_interim"

Cache   = "/home/ianchen/.vep"
Fasta   = f"{Cache}/homo_sapiens/105_GRCh38/Homo_sapiens.GRCh38.dna.toplevel.fa"
VEP_ENV = "/home/ianchen/miniconda3/envs/vep/bin"
FORK    = 16

Input  = source_dir2 / "proband.norm.pass.vcf.gz"
Output = source_dir2 / "proband.vep.txt"
Log    = source_dir2 / "proband.vep.log"

#* [1] run VEP
# --tab positional output, one line per variant thanks to --pick
cmd = [f"{VEP_ENV}/vep",
       "--input_file", str(Input), "--output_file", str(Output),
       "--species", "homo_sapiens", "--assembly", "GRCh38",
       "--cache", "--offline", "--dir_cache", Cache, "--cache_version", "105",
       "--fasta", Fasta,
       "--pick", "--pick_order", "canonical,mane,biotype,rank,length",
       "--symbol", "--canonical", "--biotype", "--numbers", "--hgvs",
       "--sift", "b", "--polyphen", "b",
       "--check_existing", "--show_ref_allele",
       "--af_gnomad", "--max_af",
       "--no_stats", "--force_overwrite", "--tab",
       "--fork", str(FORK), "--buffer_size", "20000"]
env = dict(os.environ, PATH=f"{VEP_ENV}:{os.environ['PATH']}")
print("+", " ".join(cmd), flush=True)
with open(Log, "w") as log:
    subprocess.run(cmd, check=True, env=env, stdout=log, stderr=subprocess.STDOUT)

#* [2] sanity-check the VEP output before trusting it
# (the Perl "Possible precedence issue ... IndexedBase.pm" line is a harmless warning)
bad = [l for l in open(Log) if any(k in l for k in ("ERROR", "died", "MSG:"))]
if bad:
    raise SystemExit("VEP log contains errors:\n" + "".join(bad[-5:]))
rows = sum(1 for l in open(Output) if not l.startswith("#"))
n_in = sum(1 for _ in open(source_dir2 / "proband.ab_table.txt")) - 1
print(f"VEP output rows: {rows} / input variants: {n_in}")
if rows < 0.98 * n_in:
    raise SystemExit("VEP output looks short — VEP may have been killed.")
