# Stage 03 step 3 — population frequency from gnomAD v4.1 (exomes + genomes, "joint")
# for the call-QC-passing set. The offline VEP 105 cache only carries gnomAD r2.1.1
# EXOME AF, so a variant outside the exome capture (or mis-represented there) looks
# "absent" and scores as ultra-rare; this step closes that gap.
#! Only variant coordinates (CHROM-POS-REF-ALT, GRCh38) are sent to the public
#! gnomAD GraphQL API — no genotype, depth or sample identifier.
#! Responses are cached in data/01_interim/gnomad_r4_cache.tsv (with query date),
#! so a re-run is offline and reproducible; delete the cache to re-query.
#! Operational rules (documented choices, not an established standard):
#!   - frequency used = max(v2.1.1 exome AF / popmax from VEP, v4.1 joint AF,
#!     v4.1 joint faf95 popmax); drop if > 1e-3 (same cut-off as 03.Prioritization_01)
#!   - drop a homozygous-alt candidate if gnomAD v4 has >= 1 homozygote
#!     (a severe paediatric recessive allele is not expected homozygous in gnomAD)
#!   - a site present in gnomAD v4 with joint AC = 0 (all genotypes filtered there)
#!     is flagged gnomad4_ac0 and its v4 frequency treated as unknown, not as 0
import csv
import json
import time
import urllib.error
import urllib.request
from pathlib import Path

# working directory
working_dir = Path("/home/ianchen/Documents/Projects/MAV_Hackathon/MVA")
source_dir2 = working_dir / "data/01_interim"
source_dir3 = working_dir / "data/02_processed"

API, DATASET = "https://gnomad.broadinstitute.org/api", "gnomad_r4"
BATCH, PAUSE = 20, 3                     # aliases per request, seconds between requests
AF_CUT = 1e-3
Cache = source_dir2 / "gnomad_r4_cache.tsv"
CACHE_COLS = ["variant_id", "found", "joint_ac", "joint_an", "joint_af",
              "joint_hom", "joint_faf95_popmax", "exome_af", "genome_af", "query_date"]


def af(d):
    return d["ac"] / d["an"] if d and d.get("an") else ""


def save_cache():
    with open(Cache, "w", newline="") as f:
        w = csv.DictWriter(f, CACHE_COLS, delimiter="\t")
        w.writeheader()
        w.writerows(cache[k] for k in sorted(cache))


def post(query):
    req = urllib.request.Request(API, data=json.dumps({"query": query}).encode(),
                                 headers={"Content-Type": "application/json"})
    for attempt in range(5):
        try:
            with urllib.request.urlopen(req, timeout=60) as resp:
                return json.load(resp)
        except urllib.error.HTTPError as e:
            if e.code != 429 or attempt == 4:
                raise
        except (urllib.error.URLError, TimeoutError):
            if attempt == 4:
                raise
        time.sleep(30 * (attempt + 1))                     # back off, then retry


#* [1] load the cache
cache = {}
if Cache.exists():
    for r in csv.DictReader(open(Cache), delimiter="\t"):
        cache[r["variant_id"]] = r

#* [2] query what is missing, BATCH variants per GraphQL request
rows = list(csv.DictReader(open(source_dir3 / "proband.rare_coding.qcpass.txt"), delimiter="\t"))
if not rows:
    raise SystemExit("proband.rare_coding.qcpass.txt is empty — check stage 03_02")
vid = {r["MarkerID"]: r["MarkerID"].replace("_", "-") for r in rows}
todo = sorted(set(vid.values()) - set(cache))
print("variants: ", len(vid), " cached: ", len(vid) - len(todo), " to query: ", len(todo))

fields = "ac an homozygote_count faf95 { popmax }"
for i in range(0, len(todo), BATCH):
    chunk = todo[i:i + BATCH]
    q = " ".join(f'v{j}: variant(variantId: "{v}", dataset: {DATASET}) '
                 f'{{ exome {{ {fields} }} genome {{ {fields} }} joint {{ {fields} }} }}'
                 for j, v in enumerate(chunk))
    res = post("{ " + q + " }")
    # only "Variant not found" may turn into an absent (null) record; anything else aborts
    other = [e for e in res.get("errors", []) if e.get("message") != "Variant not found"]
    if other or res.get("data") is None:
        save_cache()
        raise SystemExit(f"gnomAD API error: {other[:3]}")
    today = time.strftime("%F")
    for j, v in enumerate(chunk):
        d = res["data"].get(f"v{j}")                       # null = not in gnomAD v4
        jt = (d or {}).get("joint") or {}
        cache[v] = {
            "variant_id": v, "found": d is not None,
            "joint_ac": jt.get("ac", ""), "joint_an": jt.get("an", ""), "joint_af": af(jt),
            "joint_hom": jt.get("homozygote_count", ""),
            "joint_faf95_popmax": (jt.get("faf95") or {}).get("popmax", ""),
            "exome_af": af((d or {}).get("exome")), "genome_af": af((d or {}).get("genome")),
            "query_date": today}
    save_cache()                                           # incremental: survives a crash
    print(f"  queried {min(i + BATCH, len(todo))} / {len(todo)}", flush=True)
    time.sleep(PAUSE)

#* [3] attach gnomAD v4 and apply the frequency rules
def fnum(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


for r in rows:
    g = cache[vid[r["MarkerID"]]]
    ac0 = str(g["found"]) == "True" and fnum(g["joint_ac"]) == 0
    r["gnomad4_found"] = g["found"]
    r["gnomad4_ac0"] = ac0
    r["gnomad4_af"] = "" if ac0 else g["joint_af"]
    r["gnomad4_faf95_popmax"] = "" if ac0 else g["joint_faf95_popmax"]
    r["gnomad4_hom"] = g["joint_hom"]
    # af_use feeds the stage-04 rarity score: worst case of v2 (exome AF / popmax)
    # and v4 (joint AF, faf95 popmax); empty = not observed in either release
    afs = [fnum(x) for x in (r["af_use"], r["gnomad4_af"], r["gnomad4_faf95_popmax"])]
    afs = [x for x in afs if x is not None]
    r["af_use"] = max(afs) if afs else ""
    hom_hit = r["zyg"] == "hom_alt" and (fnum(g["joint_hom"]) or 0) >= 1
    r["qc_gnomad4"] = "FAIL" if (r["af_use"] != "" and r["af_use"] > AF_CUT) or hom_hit else "PASS"

cols = list(rows[0].keys())
with open(source_dir3 / "proband.rare_coding.gnomad.txt", "w", newline="") as f:
    w = csv.DictWriter(f, cols, delimiter="\t")
    w.writeheader()
    w.writerows(rows)
keep = [r for r in rows if r["qc_gnomad4"] == "PASS"]
with open(source_dir3 / "proband.rare_coding.analysis.txt", "w", newline="") as f:
    w = csv.DictWriter(f, cols, delimiter="\t")
    w.writeheader()
    w.writerows(keep)

print("in gnomAD v4: ", sum(str(r["gnomad4_found"]) == "True" for r in rows), " / ", len(rows),
      " (AC=0 sites: ", sum(r["gnomad4_ac0"] for r in rows), ")")
print("fail gnomAD v4 frequency / homozygote rule: ", len(rows) - len(keep))
print("analysis set: ", len(keep))
