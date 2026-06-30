#!/usr/bin/env python3
# ========================================================================
# Script: run_05_ancestry_inference_05232026_v2.py
# Pipeline 4 (ancestry / PCA), step 5.
# 処理内容:
#   - PLINK2 --score の出力 (pca_scores.sscore) を読み込み
#   - merge済みデータに含まれていない SNP (=全サンプル 0/0 想定) の
#     スコア寄与を解析的に補正 (offset_k)
#     offset_k = -Σ_{i∈missing} sqrt(2*af_i/(1-af_i)) * loading_ik
#   - PC_k = (PC_k_SUM + offset_k) / sqrt(N_total) で Hail pc_project に一致
#       N_total = 76,399 (gnomAD loadings 全体)
#   - gnomAD RF (ONNX) で 10 ancestry クラス確率を予測
#   - 出力: pca/ancestry_assignments.tsv
#   - 動的列指定: sscore のヘッダから PC_SUM 列名を取得、loadings_plink.tsvから
#     PC1-PC16の loading 列、pca_af.afreqから ALT_FREQS を取得
#   - 実行時間を計測
#
# v2変更点:
#   - merge済みデータの SNPs だけでなく、未merge SNPs の
#     constant offset 補正を追加 (全サンプル 0/0 想定)
#   - normalization を sqrt(N_used) → sqrt(N_total=76399) に修正
#
# De-identification note:
#   WORK_DIR などのパスはプレースホルダ ("/path/to/...") です。実行前に編集してください。
#   gnomAD リファレンス（loadings/RF/metadata）は公開リソースです。
#   参加者個別データは本リポジトリに含まれません（README.md 参照）。
# ========================================================================

import sys
import time
import math
from pathlib import Path
import numpy as np
import pandas as pd
import onnxruntime as rt

start_time = time.time()

# ===== Path =====
WORK_DIR = Path("/path/to/PCA_project")
SSCORE = WORK_DIR / "pca" / "pca_scores.sscore"
LOADINGS = WORK_DIR / "sitelist" / "loadings_plink.tsv"
AFREQ = WORK_DIR / "sitelist" / "pca_af.afreq"
PVAR = WORK_DIR / "merged" / "cancer_wgs.pvar"
ONNX_MODEL = WORK_DIR / "reference" / "gnomad.v3.1.RF_fit.onnx"
OUT_TSV = WORK_DIR / "pca" / "ancestry_assignments.tsv"
LOG_FILE = WORK_DIR / "logs" / "05_ancestry_inference_05232026_v2.log"

OUT_TSV.parent.mkdir(parents=True, exist_ok=True)
LOG_FILE.parent.mkdir(parents=True, exist_ok=True)

log_lines = []
def log(msg):
    line = f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {msg}"
    print(line)
    log_lines.append(line)

# ===== Step 1: PLINK2 sscore 読み込み =====
log(f"Loading sscore: {SSCORE}")
sscore = pd.read_csv(SSCORE, sep="\t")
log(f"sscore shape: {sscore.shape}, columns: {list(sscore.columns)[:5]}...")

iid_col = "#IID" if "#IID" in sscore.columns else "IID"

# PC_SUM 列を動的取得 (PC1_SUM..PC16_SUM)
pc_sum_cols = [c for c in sscore.columns if c.startswith("PC") and c.endswith("_SUM")]
pc_sum_cols = sorted(pc_sum_cols, key=lambda x: int(x.replace("PC","").replace("_SUM","")))
log(f"PC SUM columns: {pc_sum_cols}")
N_PCS = len(pc_sum_cols)
assert N_PCS == 16, f"Expected 16 PCs, got {N_PCS}"

# ===== Step 2: loadings 読み込み =====
log(f"Loading loadings: {LOADINGS}")
loadings = pd.read_csv(LOADINGS, sep="\t")
log(f"loadings shape: {loadings.shape}")
N_total_snps = len(loadings)
log(f"Total gnomAD loadings SNPs: {N_total_snps}")

# PC列名動的取得
pc_load_cols = [c for c in loadings.columns if c.startswith("PC") and not c.endswith("_SUM")]
pc_load_cols = sorted(pc_load_cols, key=lambda x: int(x.replace("PC","")))
log(f"loadings PC cols: {pc_load_cols[:5]}... (n={len(pc_load_cols)})")

# ===== Step 3: AF 読み込み =====
log(f"Loading AF: {AFREQ}")
af = pd.read_csv(AFREQ, sep="\t")
# 列名: '#CHROM', 'ID', 'REF', 'ALT', 'PROVISIONAL_REF?', 'ALT_FREQS', 'OBS_CT'
af_id_col = "ID"
af_freq_col = "ALT_FREQS"
log(f"AF columns: {list(af.columns)}, using {af_id_col} and {af_freq_col}")

# ===== Step 4: merge済みデータの SNP ID リスト =====
log(f"Reading pvar: {PVAR}")
# Skip ##header, but #CHROM line is the real header
pvar = pd.read_csv(PVAR, sep="\t", comment="#",
                   names=["CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO"])
log(f"pvar rows: {len(pvar)}")
merged_ids = set(pvar["ID"].astype(str).values)

# ===== Step 5: 未merge SNPの offset計算 =====
loadings["in_merge"] = loadings["SNPID"].isin(merged_ids)
n_in = loadings["in_merge"].sum()
n_out = (~loadings["in_merge"]).sum()
log(f"SNPs in merge: {n_in}")
log(f"SNPs NOT in merge (need offset correction): {n_out}")

# af を loadings にマージ
af_merge = af[[af_id_col, af_freq_col]].rename(columns={af_id_col:"SNPID", af_freq_col:"AF"})
loadings = loadings.merge(af_merge, on="SNPID", how="left")
log(f"Loadings with AF: {loadings['AF'].notna().sum()}/{len(loadings)}")

# 未mergeサンプルだけ抽出
missing = loadings[~loadings["in_merge"]].copy()
# 0/0 想定でのcontribution: -2*af / sqrt(2*af*(1-af)) * loading_ik = -sqrt(2*af/(1-af)) * loading_ik
# (af → 0/1 の境界は安全に扱う; af=0 or af=1の場合はskip)
missing = missing[(missing["AF"] > 0) & (missing["AF"] < 1)]
log(f"Missing with valid AF (>0, <1): {len(missing)}")

# offset per PC
af_vec = missing["AF"].values
factor = -np.sqrt(2 * af_vec / (1 - af_vec))  # shape (n_missing,)
offsets = {}
for i, col in enumerate(pc_load_cols):
    load_vec = missing[col].values
    offsets[col] = float(np.sum(factor * load_vec))
    log(f"  offset {col}: {offsets[col]:.4f}")

# ===== Step 6: 補正後の PC スコア計算 =====
N_norm = math.sqrt(N_total_snps)
log(f"Normalization factor: sqrt({N_total_snps}) = {N_norm:.4f}")

# 注意: PLINK2 sscore の PC_k_SUM は merge内 SNP の variance-standardized 寄与の合計
# Hail pc_project 相当の PC = (PC_k_SUM + offset_k) / sqrt(N_total)
df = pd.DataFrame()
df["SAMPLE_ID"] = sscore[iid_col].astype(str)
for k in range(1, N_PCS+1):
    sum_col = f"PC{k}_SUM"
    load_col = f"PC{k}"
    pc_partial = sscore[sum_col].values
    offset = offsets[load_col]
    pc_full = (pc_partial + offset) / N_norm
    df[f"PC{k}"] = pc_full

log(f"PC ranges:")
for k in range(1, 9):
    pc = df[f"PC{k}"].values
    log(f"  PC{k}: min={pc.min():.4f}, max={pc.max():.4f}, mean={pc.mean():.4f}, sd={pc.std():.4f}")

# ===== Step 7: ONNX RF 推論 =====
log(f"Loading ONNX model: {ONNX_MODEL}")
sess = rt.InferenceSession(str(ONNX_MODEL))
input_name = sess.get_inputs()[0].name

X = df[[f"PC{k}" for k in range(1, N_PCS+1)]].values.astype(np.float32)
log(f"Input shape for RF: {X.shape}")
labels, probs = sess.run(None, {input_name: X})
log(f"Inference done. labels[:5]: {labels[:5]}")

classes = sorted(probs[0].keys())
log(f"Ancestry classes: {classes}")
prob_matrix = np.array([[p[c] for c in classes] for p in probs])

df["pred_ancestry"] = labels
df["max_prob"] = prob_matrix.max(axis=1)
for i, cls in enumerate(classes):
    df[f"prob_{cls}"] = prob_matrix[:, i]

df.to_csv(OUT_TSV, sep="\t", index=False)
log(f"Saved → {OUT_TSV}")
log(f"Sample count: {len(df)}")

# Summary
log("=== Ancestry distribution ===")
for k, v in df["pred_ancestry"].value_counts().items():
    log(f"  {k}: {v} ({100*v/len(df):.1f}%)")
low = (df["max_prob"] < 0.9).sum()
log(f"Low confidence (max_prob<0.9): {low}")

elapsed = time.time() - start_time
log(f"=== Total elapsed: {elapsed:.1f} sec ===")
with open(LOG_FILE, "w") as f:
    f.write("\n".join(log_lines) + "\n")
log(f"Log saved: {LOG_FILE}")
