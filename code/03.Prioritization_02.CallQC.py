# Stage 03 step 2 — call-level QC of the rare + impactful set, BEFORE the
# inheritance models, so mapping artefacts cannot pair up as fake compound hets.
#! Operational rules (documented choices, not an established standard):
#!   [a] genotype   DP >= 10, GQ >= 20, het alt-allele balance >= 0.20
#!                  (same cut-offs as gnomAD's "adj" genotype filter),
#!                  hom-alt alt-allele balance >= 0.80
#!   [b] cluster    >= 2 other analysable variants within +-50 bp -> local
#!                  misalignment (paralogous / segmentally duplicated sequence)
#!   [c] repeat     indel inside a short tandem repeat of the reference
#!                  (unit 1 bp tract >= 6 bp, or unit 2-6 bp >= 3 copies and >= 8 bp)
# Input  data/02_processed/proband.rare_coding.txt   (kept unchanged)
# Output data/02_processed/proband.rare_coding.qc.txt      all rows + QC flags
#        data/02_processed/proband.rare_coding.qcpass.txt  rows passing all QC
import bisect
import collections
import csv
from pathlib import Path

# working directory
working_dir = Path("/home/ianchen/Documents/Projects/MVA_Hackathon/MVA")
source_dir2 = working_dir / "data/01_interim"
source_dir3 = working_dir / "data/02_processed"

Fasta = "/home/ianchen/.vep/homo_sapiens/105_GRCh38/Homo_sapiens.GRCh38.dna.toplevel.fa"

HET_AB_MIN, HOM_AB_MIN = 0.20, 0.80
MIN_DP, MIN_GQ = 10, 20
CLUSTER_BP, CLUSTER_N = 50, 2


#* [1] minimal indexed FASTA reader (.fai: name, length, offset, linebases, linewidth)
fai = {}
for line in open(Fasta + ".fai"):
    name, length, offset, linebases, linewidth = line.split("\t")[:5]
    fai[name] = tuple(int(x) for x in (length, offset, linebases, linewidth))
fa = open(Fasta, "rb")


def fetch(chrom, start, end):
    """0-based half-open reference sequence."""
    length, offset, lb, lw = fai[chrom]
    start, end = max(0, start), min(length, end)
    fa.seek(offset + (start // lb) * lw + start % lb)
    raw = fa.read((end - start) + (end - start) // lb * (lw - lb) + lw)
    return raw.replace(b"\n", b"").replace(b"\r", b"")[: end - start].decode().upper()


def str_tract(seq):
    """Longest tandem tract (bp, copies, unit) starting at seq[0], unit 1-6 bp."""
    best = (0, 0, "")
    if not seq:
        return best
    for k in range(1, 7):
        unit, n = seq[:k], 0
        while seq[n * k:(n + 1) * k] == unit:
            n += 1
        if n * k > best[0]:
            best = (n * k, n, unit)
    return best


#* [2] DISTINCT positions of every analysable variant, for the cluster rule
# (split multiallelic records share a POS and must not count as neighbours)
pos = collections.defaultdict(set)
with open(source_dir2 / "proband.ab_table.txt") as f:
    next(f)
    for line in f:
        c, p = line.split("\t", 3)[1:3]
        pos[c].add(int(p))
pos = {c: sorted(v) for c, v in pos.items()}

#* [3] flag each rare + impactful variant
rows = list(csv.DictReader(open(source_dir3 / "proband.rare_coding.txt"), delimiter="\t"))
if not rows:
    raise SystemExit("proband.rare_coding.txt is empty — check stage 03_01")


def num(x):
    try:
        return float(x)
    except ValueError:
        return float("nan")


for r in rows:
    chrom, p, ref, alt = r["MarkerID"].split("_")
    p = int(p)
    ad_ref, ad_alt = num(r["AD_REF"]), num(r["AD_ALT"])
    ab = ad_alt / (ad_ref + ad_alt) if ad_ref + ad_alt > 0 else float("nan")

    # [a] genotype quality
    dp, gq = num(r["DP"]), num(r["GQ"])
    ab_min = HET_AB_MIN if r["zyg"] == "het" else HOM_AB_MIN
    r["qc_genotype"] = "PASS" if (dp >= MIN_DP and gq >= MIN_GQ and ab >= ab_min) else "FAIL"

    # [b] neighbours within +-50 bp (excluding itself)
    v = pos[chrom]
    n_near = bisect.bisect_right(v, p + CLUSTER_BP) - bisect.bisect_left(v, p - CLUSTER_BP) - 1
    r["n_within_50bp"] = n_near
    r["qc_cluster"] = "FAIL" if n_near >= CLUSTER_N else "PASS"

    # [c] indel in a short tandem repeat: scan the reference right after the
    # (left-aligned) anchor base; for an insertion also test the inserted sequence
    r["str_unit"], r["str_bp"] = "", 0
    r["qc_repeat"] = "PASS"
    if len(ref) != len(alt):
        flank = fetch(chrom, p, p + 60)                  # 0-based p = base after anchor
        ins = alt[1:] if len(alt) > len(ref) else ""
        bp, copies, unit = max(str_tract(flank), str_tract(ins + flank))
        r["str_unit"], r["str_bp"] = unit, bp
        if (len(unit) == 1 and bp >= 6) or (len(unit) >= 2 and copies >= 3 and bp >= 8):
            r["qc_repeat"] = "FAIL"

    r["ab"] = round(ab, 3)
    r["qc_pass"] = all(r[k] == "PASS" for k in ("qc_genotype", "qc_cluster", "qc_repeat"))

#* [4] write
cols = list(rows[0].keys())
with open(source_dir3 / "proband.rare_coding.qc.txt", "w", newline="") as f:
    w = csv.DictWriter(f, cols, delimiter="\t")
    w.writeheader()
    w.writerows(rows)
with open(source_dir3 / "proband.rare_coding.qcpass.txt", "w", newline="") as f:
    w = csv.DictWriter(f, cols, delimiter="\t")
    w.writeheader()
    w.writerows([r for r in rows if r["qc_pass"]])

print("rare + impactful input: ", len(rows))
for k in ("qc_genotype", "qc_cluster", "qc_repeat"):
    print(f"  {k} FAIL: ", sum(r[k] == "FAIL" for r in rows))
print("pass all QC: ", sum(r["qc_pass"] for r in rows))
