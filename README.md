# Analysis code for germline secondary findings in tumor-normal genome sequencing

This repository contains the analysis and visualization scripts used in:

> *Auditable return governance for germline secondary findings in tumor-normal
> genome sequencing: A single-institution audit of the Japan WGS Program.*

The scripts reproduce the three computational components of the study:

1. **SNV/indel pathogenic-variant extraction** (`01_snv_indel_pathogenic/`)
2. **CNV calling and ACMG class 4/5 prioritization** (`02_cnv_acmg/`)
3. **Genetic-ancestry inference by PCA projection** (`03_ancestry_pca/`)

> **No participant-level data are included in this repository.** Input genomes,
> VCFs, BAMs, sample manifests, and individual identifiers are *not* distributed
> here. All file-system paths in the scripts are de-identified placeholders
> (`/path/to/...`) and must be edited before use. Access to the underlying data
> is subject to institutional data-use and privacy restrictions; see the
> Data Availability statement of the manuscript.

These scripts were run on an HPC cluster under a Grid Engine (SGE) scheduler.
The `#$` directives and resource requests are retained as a record of how the
analyses were executed; adapt them to your own scheduler/environment.

---

## Repository layout

```
cancerWGS-germline-secondary-findings/   # repository root
├── README.md
├── LICENSE
├── .gitignore
├── 01_snv_indel_pathogenic/
│   ├── variant_filter_parallel_v4.py   # per-sample P/LP SNV/indel extraction (array task)
│   ├── run_variant_filter.sh           # SGE array-job wrapper
│   ├── merge_results.py                # merge per-batch tables into a cohort table
│   └── merge_results.sh                # SGE wrapper for the merge step
├── 02_cnv_acmg/
│   ├── cnvkit_germline_array_v4.sh     # CNVkit WGS calling (array task, one per sample)
│   ├── convert_cnvkit_to_bed_v2.py     # CNVkit .call.cns -> AnnotSV BED
│   ├── run_annotsv_all_samples_v2.sh   # AnnotSV annotation (assigns ACMG/ClinGen class)
│   ├── filter_cnvs_v3.py               # ACMG class 4/5 + autosomal + segdup<70% selection
│   └── install_annotsv_local.sh        # local AnnotSV v3.3.6 install helper
└── 03_ancestry_pca/
    ├── run_01_prepare_loadings_05232026.py        # gnomAD v3.1 loadings -> PLINK2 inputs
    ├── run_02_extract_subset_vcfs_05232026_v2.sh  # subset normal VCFs to loading SNPs (array)
    ├── run_03_merge_and_plink_05232026_v2.sh      # merge subsets -> PLINK2 pgen
    ├── run_04_pca_projection_05232026.sh          # PLINK2 --score projection
    ├── run_05_ancestry_inference_05232026_v2.py   # offset correction + gnomAD RF (ONNX) inference
    ├── run_06_visualize_05232026.py               # PCA / ancestry plots
    ├── run_07_overlay_with_reference_05232026_v2.py  # overlay cohort on HGDP+1KG reference
    └── run_08_eas_substructure_05232026.py        # EAS sub-structure / Japanese-cluster check
```

Filenames retain their internal version suffixes (`_v4`, `_v2`, `_v7`) so they
map one-to-one to the exact scripts that produced the reported results.

---

## Pipeline 1 — SNV/indel pathogenic-variant extraction

**Input.** Per-sample ANNOVAR-annotated tables of the *normal* DeepVariant calls
(`*.deepvariant.vcf.gz.hg38_multianno.txt`), an ACMG/return-gene list
(`SF_106Genes.txt`), the ClinVar `variant_summary_GRCh38.txt.gz` table, and a
sample manifest.

**Procedure.** For each record the script keeps `FILTER == PASS`, `GQ ≥ 20`,
`DP ≥ 10`, then flags Pathogenic/Likely-pathogenic by either ClinVar
(excluding *Conflicting* classifications) or `InterVar_automated`, and annotates
ClinVar VariationID / clinical significance / review-status stars and SF-gene
membership. ToMMo 38KJPN frequency is carried as an annotation column and applied
as a downstream filter (it is **not** an automated cut-off in this script).
All columns are resolved by header **name**, not by fixed position.

**Output.** Per-sample filtered tables, per-batch statistics, and a merged
cohort table (`all_samples_pathogenic_variants_<timestamp>.txt`).

**Run.**
```bash
# array job (set -t to ceil(N_samples / batch_size))
qsub 01_snv_indel_pathogenic/run_variant_filter.sh
# then merge
qsub 01_snv_indel_pathogenic/merge_results.sh
```

## Pipeline 2 — CNV calling and ACMG class 4/5 prioritization

**Input.** Per-sample normal BAMs, a pooled CNVkit reference (`*.cnn`), a CNVkit
container image (Apptainer/Singularity), AnnotSV v3.3.6 with its annotation
databases, and a merged segmental-duplication track (`segdup_hg38_merged.txt`).

**Procedure.** CNVkit (WGS mode) → CBS segmentation (threshold 0.1) → call
(ploidy 2) per sample; convert to BED; annotate with AnnotSV (GRCh38), which
assigns the ACMG/ClinGen SV classification; then `filter_cnvs_v3.py` keeps only
AnnotSV **ACMG class 4 (likely pathogenic) / class 5 (pathogenic)** events,
excludes sex-chromosome calls, and removes CNVs overlapping segmental
duplications by ≥ 70 %. This yields the high-priority CNV workload reported in
the manuscript. The `ACMG_class` column is selected by name; the retained
classes are configurable with `--acmg-classes` (default `4,5`).

> **Note — filter version.** `filter_cnvs_v3.py` matches the CNV criteria stated
> in the manuscript (ACMG class 4/5, autosomal, segdup < 70 %). A later
> `filter_cnvs_v4.py` additionally removes common CNVs by ToMMo JCNVv1 (54KJPN)
> population frequency; that extra filter is **not** part of the manuscript's
> 82-event CNV definition and is therefore not included here.

## Pipeline 3 — Genetic-ancestry inference (PCA projection)

**Input.** Per-sample normal DeepVariant VCFs and public gnomAD v3.1 resources:
PCA loadings (`gnomad.v3.1.pca_loadings.tsv.gz`), the random-forest ancestry
classifier (`gnomad.v3.1.RF_fit.onnx`), and HGDP+1KG sample metadata
(`gnomad_meta_v1.tsv`).

**Procedure.** Build PLINK2 score inputs from the loadings → subset each sample's
VCF to the loading SNPs → merge → project with `PLINK2 --score`
(variance-standardized, Hail `pc_project`–equivalent, with an analytic offset
correction for sites absent from the merge) → classify with the gnomAD RF model
→ visualize and check East-Asian sub-structure relative to the HGDP+1KG Japanese
cluster.

> **Note — cohort count.** Array sizes and some hard-coded counts in the PCA
> scripts (e.g. `-t 1-248`, "n = 247") reflect the sample manifest at the time
> the ancestry analysis was run and may include samples beyond the final
> analytic cohort (N = 197). They do not change the projection logic. Adjust to
> your cohort when reproducing.

---

## Dependencies

The analyses were run on Linux (HPC, Grid Engine). Versions used:

| Tool | Version | Used in |
|------|---------|---------|
| Python | 3.x | all pipelines |
| pandas, numpy | — | 1, 2, 3 |
| intervaltree | — | 2 (`filter_cnvs_v3.py`) |
| scipy | — | 3 (`run_08`) |
| onnxruntime | — | 3 (`run_05`) |
| matplotlib | — | 3 (visualization) |
| DeepVariant | (upstream caller) | input VCFs |
| ANNOVAR + InterVar | — | 1 (annotation/classification, upstream) |
| ClinVar `variant_summary` (GRCh38) | 2021-10 anchor + later refresh | 1 |
| CNVkit | 0.9 | 2 |
| AnnotSV | 3.3.6 | 2 |
| bedtools / samtools | 2.31.1 / 1.19 | 2 |
| htslib / bcftools (bgzip, tabix) | 1.19 | 3 |
| PLINK2 | — | 3 |
| gnomAD v3.1 loadings / RF model / HGDP+1KG metadata | v3.1 | 3 |
| ToMMo 38KJPN allele frequencies | 38KJPN | 1 (annotation) |

CNVkit is executed from a container image; AnnotSV can be installed locally with
`02_cnv_acmg/install_annotsv_local.sh`.

---

## De-identification

Before publication, the following were removed or replaced with placeholders:
internal file-system paths, cluster/queue/group names, usernames, host names,
and example sample identifiers. Scripts that emit per-sample tables or labelled
figures (e.g. `03_ancestry_pca/run_08_*`) should be published **without** their
output directories, since those outputs can contain sample IDs.

## License

Released under the MIT License (see `LICENSE`).

## Citation

Please cite the associated manuscript. For questions about the code, contact the
corresponding author.
