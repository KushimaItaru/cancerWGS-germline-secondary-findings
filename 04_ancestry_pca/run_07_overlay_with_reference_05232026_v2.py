#!/usr/bin/env python3
# ========================================================================
# Script: run_07_overlay_with_reference_05232026_v2.py
# Pipeline 4 (ancestry / PCA), step 7.
# 処理内容:
#   - gnomAD metadata の population_inference.pca_scores 列 (JSON配列)から
#     HGDP+1KG samplesの projection scale PC scoresを取得
#       (これは cancer WGS scoring と同じスケール)
#   - sample metadata の project_meta.project_pop でsuper-pop label付与
#   - cancer WGS samples を ancestry_assignments.tsv から読み込み
#   - 同一PC空間にプロット (PC1-PC2, PC2-PC3, PC3-PC4)
#       reference は色付きdots、cancer WGS は黒xマーカ
#   - 数値検証: cancer の PC1-4 が EAS reference の範囲と重なるか
#
# v2変更点:
#   - GLOBAL_pc_scores.txt.bgz (original PCA scale) ではなく
#     population_inference.pca_scores (projection scale) を使用
#   - これで projection結果と完全に同じスケールでプロット可能
#
# De-identification note:
#   WORK_DIR などのパスはプレースホルダ ("/path/to/...") です。実行前に編集してください。
#   gnomAD リファレンス（metadata）は公開リソースです。
#   参加者個別データは本リポジトリに含まれません（README.md 参照）。
# ========================================================================

import sys, time, json
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

start_time = time.time()

WORK_DIR = Path("/path/to/PCA_project")
REF_META = WORK_DIR / "reference" / "gnomad_meta_v1.tsv"
ANC_TSV  = WORK_DIR / "pca" / "ancestry_assignments.tsv"
OUT_DIR  = WORK_DIR / "pca"
LOG_FILE = WORK_DIR / "logs" / "07_overlay_v2_05232026.log"

OUT_DIR.mkdir(parents=True, exist_ok=True)
LOG_FILE.parent.mkdir(parents=True, exist_ok=True)

log_lines = []
def log(msg):
    line = f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {msg}"
    print(line)
    log_lines.append(line)

# ===== Step 1: metadata からprojection scale PC scoresを取得 =====
log(f"Loading metadata: {REF_META}")
header = pd.read_csv(REF_META, sep="\t", nrows=0).columns.tolist()

# 動的列指定
sample_col = "s"
pop_col = [c for c in header if c == "project_meta.project_pop"][0]
infer_pop_col = "population_inference.pop"
pca_col = "population_inference.pca_scores"

use_cols = [sample_col, pop_col, infer_pop_col, pca_col]
log(f"Using columns: {use_cols}")

meta = pd.read_csv(REF_META, sep="\t", usecols=use_cols)
log(f"Metadata loaded: {len(meta)} samples")

# Parse pca_scores JSON array
def parse_pca_scores(s):
    if pd.isna(s) or s == "NA":
        return None
    try:
        arr = json.loads(s)
        return arr
    except:
        return None

log("Parsing PCA scores JSON...")
parsed = meta[pca_col].apply(parse_pca_scores)
valid_mask = parsed.notna()
log(f"Valid PCA scores: {valid_mask.sum()} / {len(meta)}")

# Determine PC count from first valid entry
pc_lens = parsed[valid_mask].apply(len).value_counts()
log(f"PC count distribution: {pc_lens.to_dict()}")
N_PCS = pc_lens.index[0]
log(f"Using {N_PCS} PCs")

# Build df_ref
df_ref = pd.DataFrame()
df_ref["SAMPLE_ID"] = meta[sample_col]
for k in range(N_PCS):
    df_ref[f"PC{k+1}"] = parsed.apply(lambda x: x[k] if x is not None else np.nan)

# Add pop labels
df_ref["project_pop"] = meta[pop_col].values
df_ref["inferred_pop"] = meta[infer_pop_col].values
df_ref = df_ref[valid_mask.values].reset_index(drop=True)
log(f"Reference with valid scores: {len(df_ref)}")

# Use project_pop as primary label
df_ref["pop_label"] = df_ref["project_pop"].fillna(df_ref["inferred_pop"]).fillna("unknown")
log(f"Pop label distribution:")
for k, v in df_ref["pop_label"].value_counts().items():
    log(f"  {k}: {v}")

# ===== Step 2: Cancer samples 読み込み =====
df_can = pd.read_csv(ANC_TSV, sep="\t")
log(f"\nCancer WGS shape: {df_can.shape}")
log(f"Cancer pred_ancestry: {df_can['pred_ancestry'].value_counts().to_dict()}")

# ===== Step 3: 数値検証 (projection scale) =====
log(f"\n=== Sanity check: cancer vs EAS reference (projection scale) ===")
eas_ref = df_ref[df_ref["pop_label"].str.lower()=="eas"]
log(f"EAS reference samples: {len(eas_ref)}")
for k in [1,2,3,4,5,6,7,8]:
    pc = f"PC{k}"
    log(f"  {pc}:")
    log(f"    EAS ref:  range [{eas_ref[pc].min():.4f}, {eas_ref[pc].max():.4f}]  mean={eas_ref[pc].mean():.4f} ± {eas_ref[pc].std():.4f}")
    log(f"    Cancer:   range [{df_can[pc].min():.4f}, {df_can[pc].max():.4f}]  mean={df_can[pc].mean():.4f} ± {df_can[pc].std():.4f}")

# ===== Step 4: 重ね描き =====
pop_colors = {
    "afr": "#1f77b4", "amr": "#d62728", "ami": "#9467bd",
    "asj": "#ff7f0e", "eas": "#2ca02c", "csa": "#17becf",
    "sas": "#17becf", "eur": "#e377c2", "nfe": "#e377c2",
    "fin": "#bcbd22", "mid": "#8c564b", "oce": "#7f7f7f",
    "oth": "#bbbbbb", "unknown": "#cccccc",
}

def plot_pc(ax, pcx, pcy, ref_df, can_df, title):
    # reference: small colored dots
    for pop, grp in ref_df.groupby("pop_label"):
        ax.scatter(grp[pcx], grp[pcy],
                   c=pop_colors.get(pop, "gray"),
                   s=14, alpha=0.5,
                   label=f"{pop} ref (n={len(grp)})",
                   edgecolors="none")
    # cancer WGS: large black X
    ax.scatter(can_df[pcx], can_df[pcy],
               c="black", marker="x", s=35, alpha=0.7,
               linewidths=1.0,
               label=f"Cancer WGS (n={len(can_df)})")
    ax.set_xlabel(pcx, fontsize=12)
    ax.set_ylabel(pcy, fontsize=12)
    ax.set_title(title, fontsize=12)
    ax.legend(fontsize=7, loc="best", markerscale=1.5)
    ax.grid(True, alpha=0.3)

# Figure 1: PC1-PC2, PC2-PC3, PC3-PC4 (overall)
fig, axes = plt.subplots(1, 3, figsize=(24, 7.5))
plot_pc(axes[0], "PC1", "PC2", df_ref, df_can, "PC1 vs PC2")
plot_pc(axes[1], "PC2", "PC3", df_ref, df_can, "PC2 vs PC3")
plot_pc(axes[2], "PC3", "PC4", df_ref, df_can, "PC3 vs PC4")
plt.suptitle(
    f"Cancer WGS (n={len(df_can)}) projected onto gnomAD HGDP+1KG (n={len(df_ref)}) PCA space [projection scale]",
    fontsize=14)
plt.tight_layout()
out_png1 = OUT_DIR / "pca_overlay_v2_PC1234_05232026.png"
plt.savefig(out_png1, dpi=130, bbox_inches="tight")
plt.close()
log(f"\nSaved: {out_png1}")

# Figure 2: focused PC1-PC2
fig, ax = plt.subplots(figsize=(11, 9))
plot_pc(ax, "PC1", "PC2", df_ref, df_can,
        "PC1 vs PC2 — Cancer WGS on HGDP+1KG reference (projection scale)")
plt.tight_layout()
out_png2 = OUT_DIR / "pca_overlay_v2_PC1_PC2_05232026.png"
plt.savefig(out_png2, dpi=150, bbox_inches="tight")
plt.close()
log(f"Saved: {out_png2}")

# Figure 3: zoom on EAS region
fig, ax = plt.subplots(figsize=(11, 9))
plot_pc(ax, "PC1", "PC2", df_ref, df_can, "PC1 vs PC2 (zoomed on EAS region)")
# zoom: cancer samples region + a bit
pc1_lo = min(df_can["PC1"].min(), eas_ref["PC1"].min()) - 0.01
pc1_hi = max(df_can["PC1"].max(), eas_ref["PC1"].max()) + 0.01
pc2_lo = min(df_can["PC2"].min(), eas_ref["PC2"].min()) - 0.01
pc2_hi = max(df_can["PC2"].max(), eas_ref["PC2"].max()) + 0.01
ax.set_xlim(pc1_lo, pc1_hi)
ax.set_ylim(pc2_lo, pc2_hi)
plt.tight_layout()
out_png3 = OUT_DIR / "pca_overlay_v2_PC1_PC2_zoom_eas_05232026.png"
plt.savefig(out_png3, dpi=150, bbox_inches="tight")
plt.close()
log(f"Saved: {out_png3}")

elapsed = time.time() - start_time
log(f"=== Total elapsed: {elapsed:.1f} sec ===")
with open(LOG_FILE, "w") as f:
    f.write("\n".join(log_lines) + "\n")
log(f"Log saved: {LOG_FILE}")
