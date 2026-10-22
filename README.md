# MVA Hackathon 2026 — Track 1 variant prioritisation

Team **Nancyandher3littleelves** · *Rare Disease, Real Kid: The MVA Hackathon 2026*, Track 1 (Variant Identification)

A single-sample WGS prioritisation pipeline that goes from the proband VCF to the hackathon submission CSV. Each step is a Python or R script, run in order by one driver. Annotation is offline (Ensembl VEP 105 cache). The only network call is a gnomAD v4.1 frequency lookup that sends variant IDs (chrom-pos-ref-alt) and nothing else, and its results are cached.

> **No data in this repository.** The hackathon dataset is gated and may not be redistributed. The repository holds code and documentation only. `.gitignore` is a whitelist (`code/`, `README.md`, `.gitignore`), so VCFs, annotations, candidate tables, the gnomAD cache and the submission CSV are never committed.

## Pipeline

| Step | Script | What it does |
|---|---|---|
| 01 | `01.QC_01.Normalise.py` | bcftools: keep PASS sites on chr1–22, X, Y; split multiallelics and left-normalise; drop hom-ref and missing genotypes; require `FMT/DP ≥ 8` and `GQ ≥ 20`; set ID = `CHROM_POS_REF_ALT`; write the genotype / allele-depth table |
| 01 | `01.QC_02.Stats.R` | FILTER breakdown, bcftools stats, het/hom ratio, sex check from chrX/chrY |
| 02 | `02.Annotation_01.RunVEP.py` | Ensembl VEP 105, offline cache, GRCh38, `--pick` (one consequence per variant), SIFT, PolyPhen, gnomAD r2.1.1 exome AF, ClinVar (cache release) |
| 02 | `02.Annotation_02.ParseVEP.R` | Parse VEP, rank consequence severity, flag ClinVar pathogenic (mixed benign/pathogenic assertions excluded), attach the genotype |
| 03 | `03.Prioritization_01.Filter.R` | Rarity (max AF ≤ 1e-3 or absent) + impactful consequence (HIGH/MODERATE, splice-region classes) or ClinVar pathogenic; protein-coding only |
| 03 | `03.Prioritization_02.CallQC.py` | Call-level QC: DP ≥ 10, GQ ≥ 20, het allele balance ≥ 0.20, hom allele balance ≥ 0.80; drop variants with ≥ 2 other variant positions within ±50 bp (misalignment clusters); drop indels inside short tandem repeats of the reference |
| 03 | `03.Prioritization_03.gnomAD.py` | gnomAD v4.1 joint (exome + genome) AF and faf95 popmax from the public GraphQL API, **variant IDs (chrom-pos-ref-alt) only**, cached; drop AF > 1e-3; drop hom-alt candidates that have any gnomAD homozygote |
| 03 | `03.Prioritization_04.InheritanceModels.R` | Proband-only models: homozygous, potential compound heterozygous (≥ 2 rare het variants in one gene, phase unknown), X-linked (GQ ≥ 30, DP ≥ 10), private heterozygous (absent from gnomAD v2 and v4; de novo status cannot be tested without parents) |
| 04 | `04.Phenotype_01.GenePanels.R` | Curated gene panels for the proband's HPO profile (MVA, chromosomal instability, cancer predisposition, primordial dwarfism, nephrocalcinosis, short stature) and a gene → phenotype-feature map |
| 04 | `04.Phenotype_02.Rank.R` | Additive score = consequence + rarity + ClinVar + panel weight + inheritance model + phenotype overlap |
| 05 | `05.Report_01.SubmissionCSV.R` | Group candidates into units (one compound-het pair per gene = one row), exclude highly polymorphic gene families, keep one row per gene, put panel genes first, compute EPCR, write up to 10 rows in the hackathon format |

`code/00.Run_All.py` runs all steps in order. At the start of each run it logs the Python, bcftools, VEP and R / data.table versions.

## Requirements

- bcftools 1.24 (conda env)
- Ensembl VEP 105 with the offline GRCh38 cache 105 and its toplevel FASTA (`~/.vep/homo_sapiens/105_GRCh38/`)
- R 4.5 with `data.table`, `dplyr`, `stringr`
- Python ≥ 3.9, standard library only
- Network access for the first run of the gnomAD step; later runs use the cache

The absolute paths at the top of each script (`working_dir`, `Fasta`, tool binaries) are machine-specific. Edit them before running anywhere else.

## Running

1. Request access to the gated dataset on Hugging Face (`SageBio/mva-hackathon-2026-data`).
2. Put `WGS_EX2312012_HGWCNDSX7.vcf.gz` and its `.tbi` in `data/raw/`.
3. Run:

```bash
python3 code/00.Run_All.py                       # full run (~20 min, mostly VEP)
python3 code/00.Run_All.py --from 03.Prioritization_02.CallQC.py   # partial re-run
```

The submission is written to `out/Nancyandher3littleelves_panelrank.csv`. Intermediate tables are in `data/01_interim/` and `data/02_processed/`.

## Design choices and limitations

- **All thresholds are operational choices, not established standards.** This covers the AF cut-off, the DP/GQ/allele-balance cut-offs, the cluster and tandem-repeat rules, the gene-family exclusion, the panel weights and the score weights.
- **The gene panels and the gene → feature map are curated by hand** for this proband's phenotype. They are not derived automatically from HPO.
- **EPCR is not calibrated.** It is a softmax of the integer score, forced to be non-increasing down the rows. Only the ranking is scored by the hackathon.
- **Compound-heterozygous pairs are unphased.** The VCF has no phasing information for distant variants and there are no parental samples, so a pair is a candidate *conditional on the variants being in trans*.
- **The pipeline only sees SNVs and small indels from the provided VCF.** It does not assess SVs, CNVs, repeat expansions, mitochondrial variants or low-level mosaicism, and it does not consider non-coding variants beyond the splice-region classes.
- **Variant-effect evidence is limited** to VEP consequence, SIFT, PolyPhen and the cached ClinVar release. No SpliceAI, REVEL, CADD or functional data were used.

## Data handling and AI use

- The dataset is used under the hackathon Data Use rules and will be deleted within 30 days of the hackathon close.
- Apart from the AI-assistant session noted below, the only external transfer is variant IDs (CHROM-POS-REF-ALT) sent to the public gnomAD API. No genotypes, depths or sample identifiers are sent.
- Code was developed with **Anthropic Claude Code (Claude Team/Enterprise plan, commercial terms, no training on customer content)**. Variant-level intermediate results were visible in the assistant session.

## Acknowledgement

This work was made possible through the Hackathon, organized by Sage Bionetworks in partnership with the MVA Society, Hugging Face, and BEACON (The Benchmarking, Evaluation, and Assessment Consortium for Science), with prize sponsorship from AWS and Anthropic. We are deeply grateful to the child and their family who generously contributed their data and their story to advance research into this rare disease. We acknowledge their trust in making this Hackathon possible.
