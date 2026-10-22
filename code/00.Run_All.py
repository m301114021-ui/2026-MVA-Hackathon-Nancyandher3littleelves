# Driver: run every stage in order, from the raw VCF to the Track 1 submission CSV.
# All steps are Python or R (no shell scripts).
#   python3 code/00.Run_All.py                 # full run
#   python3 code/00.Run_All.py --from 02.Annotation_02.ParseVEP.R
# Each step overwrites its own outputs, so a partial re-run is safe.
# Launch detached so it survives an SSH / VPN drop:
#   setsid nohup python3 code/00.Run_All.py > data/01_interim/run_all.log 2>&1 &
import argparse
import os
import subprocess
import sys
import time
from pathlib import Path

working_dir = Path("/home/ianchen/Documents/Projects/MAV_Hackathon/MVA")
code_dir    = working_dir / "code"
Rscript     = "/home/ianchen/miniconda3/envs/R4.5.2/bin/Rscript"

STEPS = [
    "01.QC_01.Normalise.py",               # bcftools: PASS, norm, DP/GQ gate, AB table
    "01.QC_02.Stats.R",                    # FILTER breakdown, het/hom, sex check
    "02.Annotation_01.RunVEP.py",          # VEP 105 offline (~20-30 min)
    "02.Annotation_02.ParseVEP.R",         # VEP + genotype -> proband.annot.txt.gz
    "03.Prioritization_01.Filter.R",       # rarity + consequence filter
    "03.Prioritization_02.CallQC.py",      # genotype / cluster / tandem-repeat QC
    "03.Prioritization_03.gnomAD.py",      # gnomAD v4.1 joint AF (coordinates only)
    "03.Prioritization_04.InheritanceModels.R",
    "04.Phenotype_01.GenePanels.R",        # curated HPO gene panels
    "04.Phenotype_02.Rank.R",              # score & rank
    "05.Report_01.SubmissionCSV.R",        # hackathon Track 1 CSV
]


def log(msg):
    print(f"[{time.strftime('%F %T')}] {msg}", flush=True)


ap = argparse.ArgumentParser()
ap.add_argument("--from", dest="start", default=STEPS[0], choices=STEPS)
ap.add_argument("--to", dest="end", default=STEPS[-1], choices=STEPS)
args = ap.parse_args()

#* record tool versions with every run
VERSIONS = {
    "python":   [sys.executable, "--version"],
    "bcftools": ["/home/ianchen/.conda/envs/bcftools/bin/bcftools", "--version-only"],
    "vep":      ["/home/ianchen/miniconda3/envs/vep/bin/vep", "--help"],
    "R":        [Rscript, "-e", 'cat(R.version.string, "| data.table",'
                 ' as.character(packageVersion("data.table")), "\\n")'],
}
for tool, cmd in VERSIONS.items():
    out = subprocess.run(cmd, capture_output=True, text=True,
                         env=dict(os.environ, PATH=f"{Path(cmd[0]).parent}:{os.environ['PATH']}"))
    lines = (out.stdout + out.stderr).strip().splitlines()
    hit = [l for l in lines if "ensembl-vep" in l] or lines[:1]
    log(f"version {tool}: {' '.join(hit[0].split())}")

for step in STEPS[STEPS.index(args.start):STEPS.index(args.end) + 1]:
    exe = [sys.executable, "-u"] if step.endswith(".py") else [Rscript, "--no-save"]
    log(f"RUN {step}")
    t0 = time.time()
    if subprocess.run(exe + [str(code_dir / step)], cwd=working_dir).returncode != 0:
        log(f"!! FAILED {step}")
        sys.exit(1)
    log(f"OK  {step} ({time.time() - t0:.0f}s)")

log("DONE")
