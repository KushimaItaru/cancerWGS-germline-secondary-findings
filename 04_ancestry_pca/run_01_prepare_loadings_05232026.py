#!/usr/bin/env python3
# ========================================================================
# Script: run_01_prepare_loadings_05232026.py
# Pipeline 4 (ancestry / PCA), step 1.
# 処理内容:
#   - gnomAD v3.1 PCA loadings TSV (locus, alleles, loadings配列, pca_af)
#     を読み込み、PCAプロジェクション用に変換
#   - 出力1: sites.bed (chr, start, end) - 76,399 SNP位置 (0-based BED)
#   - 出力2: sites.tsv (chr, pos, REF, ALT, SNPID=chr:pos:REF:ALT) - VCF抽出用
#   - 出力3: loadings_plink.tsv (PLINK2 --score 用)
#       header: SNPID  A1(=ALT, 効果allele)  A2(=REF)  PC1..PC16
#   - 出力4: pca_af.afreq (PLINK2 --read-freq 用)
#       header: #CHROM  ID  REF  ALT  PROVISIONAL_REF?  ALT_FREQS  OBS_CT
#   - 30PCのうち最初の16PCのみ使用 (gnomAD RF分類器が16次元入力)
#   - 実行時間を計測してログ出力
#
# De-identification note:
#   WORK_DIR などのパスはプレースホルダ ("/path/to/...") です。実行前に編集してください。
#   gnomAD リファレンス（loadings/RF/metadata）は公開リソースです。
#   参加者個別データは本リポジトリに含まれません（README.md 参照）。
# ========================================================================

import gzip
import json
import time
from pathlib import Path

start_time = time.time()

# ===== パス設定 =====
WORK_DIR = Path("/path/to/PCA_project")
LOADINGS_TSV = WORK_DIR / "reference" / "gnomad.v3.1.pca_loadings.tsv.gz"

OUT_BED = WORK_DIR / "sitelist" / "sites.bed"
OUT_SITES_TSV = WORK_DIR / "sitelist" / "sites.tsv"
OUT_LOADINGS = WORK_DIR / "sitelist" / "loadings_plink.tsv"
OUT_AFREQ = WORK_DIR / "sitelist" / "pca_af.afreq"
OUT_LOG = WORK_DIR / "logs" / "01_prepare_loadings_05232026.log"

N_PCS = 16  # gnomAD RF分類器が要求する次元

# 出力ディレクトリ作成
for p in (OUT_BED, OUT_SITES_TSV, OUT_LOADINGS, OUT_AFREQ, OUT_LOG):
    p.parent.mkdir(parents=True, exist_ok=True)

log_lines = []

def log(msg):
    line = f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {msg}"
    print(line)
    log_lines.append(line)

log(f"Reading loadings from {LOADINGS_TSV}")

# ===== loadings.tsv.gz を読み込み =====
# header: locus  alleles  loadings  pca_af
# locus形式: "chr1:930204"
# alleles形式: '["G","A"]'  (REF=最初, ALT=2番目)
# loadings形式: '[v1,v2,...,v30]' (30PC)
# pca_af: ALT allele frequency in training set

# Header column resolution (動的にカラム位置を取る)
with gzip.open(LOADINGS_TSV, "rt") as f:
    header_line = f.readline().rstrip("\n")
    header = header_line.split("\t")
    col_locus = header.index("locus")
    col_alleles = header.index("alleles")
    col_loadings = header.index("loadings")
    col_af = header.index("pca_af")
    log(f"Header columns: locus={col_locus}, alleles={col_alleles}, loadings={col_loadings}, pca_af={col_af}")

    n_total = 0
    n_kept = 0

    bed_fh = open(OUT_BED, "w")
    sites_fh = open(OUT_SITES_TSV, "w")
    sites_fh.write("CHROM\tPOS\tREF\tALT\tSNPID\n")
    load_fh = open(OUT_LOADINGS, "w")
    load_fh.write("SNPID\tA1\tA2\t" + "\t".join([f"PC{i+1}" for i in range(N_PCS)]) + "\n")
    af_fh = open(OUT_AFREQ, "w")
    af_fh.write("#CHROM\tID\tREF\tALT\tPROVISIONAL_REF?\tALT_FREQS\tOBS_CT\n")

    for line in f:
        n_total += 1
        parts = line.rstrip("\n").split("\t")
        locus = parts[col_locus]
        alleles_json = parts[col_alleles]
        loadings_json = parts[col_loadings]
        pca_af = float(parts[col_af])

        # locus → chr, pos
        chrom, pos_str = locus.split(":")
        pos = int(pos_str)

        # alleles JSON: ["G","A"]
        alleles = json.loads(alleles_json)
        if len(alleles) != 2:
            continue
        ref, alt = alleles[0], alleles[1]

        # bi-allelic SNP のみ (insertion/deletionも理論上含むがgnomAD loadingsはほぼSNV)
        # SNVの定義: REFもALTも1塩基
        if len(ref) != 1 or len(alt) != 1:
            continue
        if ref not in "ACGT" or alt not in "ACGT":
            continue

        # loadings array
        loadings = json.loads(loadings_json)
        if len(loadings) < N_PCS:
            continue
        loadings_pc = loadings[:N_PCS]

        # SNPID: chr:pos:REF:ALT (PLINK2標準)
        snpid = f"{chrom}:{pos}:{ref}:{alt}"

        # BED (0-based, half-open)
        bed_fh.write(f"{chrom}\t{pos-1}\t{pos}\n")

        # sites TSV
        sites_fh.write(f"{chrom}\t{pos}\t{ref}\t{alt}\t{snpid}\n")

        # loadings PLINK2 score format
        # A1 = ALT (effect allele), A2 = REF
        loadings_str = "\t".join([f"{v:.10g}" for v in loadings_pc])
        load_fh.write(f"{snpid}\t{alt}\t{ref}\t{loadings_str}\n")

        # AF file (.afreq PLINK2 format)
        # PROVISIONAL_REF? は 'N' (not provisional)
        af_fh.write(f"{chrom}\t{snpid}\t{ref}\t{alt}\tN\t{pca_af:.10g}\t8200\n")
        # OBS_CT = 8200 はHGDP+1KG callset (4150 samples * 2 alleles) 目安

        n_kept += 1

    bed_fh.close()
    sites_fh.close()
    load_fh.close()
    af_fh.close()

log(f"Total loadings lines processed: {n_total}")
log(f"SNPs kept (bi-allelic SNV): {n_kept}")
log(f"Filtered out: {n_total - n_kept}")
log("")
log(f"Output files:")
log(f"  BED:      {OUT_BED}")
log(f"  Sites:    {OUT_SITES_TSV}")
log(f"  Loadings: {OUT_LOADINGS}")
log(f"  AF:       {OUT_AFREQ}")

# BED をソート + マージ (bcftools view -R 用)
import subprocess
sorted_bed = OUT_BED.with_suffix(".sorted.bed")
log(f"Sorting BED → {sorted_bed}")
subprocess.run(
    f"sort -k1,1V -k2,2n {OUT_BED} > {sorted_bed} && mv {sorted_bed} {OUT_BED}",
    shell=True, check=True
)

elapsed = time.time() - start_time
log(f"=== Total elapsed time: {elapsed:.1f} sec ({elapsed/60:.2f} min) ===")

# ログ保存
with open(OUT_LOG, "w") as f:
    f.write("\n".join(log_lines) + "\n")
print(f"Log saved: {OUT_LOG}")
