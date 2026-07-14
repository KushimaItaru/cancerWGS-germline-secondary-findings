#!/usr/bin/env python3
# ========================================================================
# Script: run_06_visualize_05232026.py
# Pipeline 4 (ancestry / PCA), step 6.
# 処理内容:
#   - ancestry_assignments.tsv から PC1-PC16 + 予測ancestry を読み込み
#   - 可視化:
#     1. PC1 vs PC2 散布図 (ancestry色分け)
#     2. PC2 vs PC3 散布図
#     3. PC3 vs PC4 散布図
#     4. ancestry分布 棒グラフ
#     5. max_prob ヒストグラム (確信度分布)
#   - サンプルIDの中で max_prob < 0.9 のものを別マーカで強調 (outlier候補)
#   - 出力: pca/pca_plot_*.png
#   - 動的列指定: PC列名・確率列名を columns から抽出
#   - 実行時間を計測
#
# De-identification note:
#   WORK_DIR などのパスはプレースホルダ ("/path/to/...") です。実行前に編集してください。
#   参加者個別データは本リポジトリに含まれません（README.md 参照）。
# ========================================================================

import time
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

start_time = time.time()

# ===== Path =====
WORK_DIR = Path("/path/to/PCA_project")
ANC_TSV = WORK_DIR / "pca" / "ancestry_assignments.tsv"
OUT_DIR = WORK_DIR / "pca"
LOG_FILE = WORK_DIR / "logs" / "06_visualize_05232026.log"

OUT_DIR.mkdir(parents=True, exist_ok=True)
LOG_FILE.parent.mkdir(parents=True, exist_ok=True)

log_lines = []
def log(msg):
    line = f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {msg}"
    print(line)
    log_lines.append(line)

# ===== 読み込み =====
df = pd.read_csv(ANC_TSV, sep="\t")
log(f"Loaded: {ANC_TSV}, shape={df.shape}")
log(f"Columns: {list(df.columns)}")

# PC列を動的取得 (PC1, PC2, ..., PC16)
pc_cols = sorted([c for c in df.columns if c.startswith("PC") and c[2:].isdigit()],
                 key=lambda x: int(x[2:]))
log(f"PC columns: {pc_cols}")

# 確率列を動的取得 (prob_*)
prob_cols = [c for c in df.columns if c.startswith("prob_")]
log(f"Probability columns: {prob_cols}")

# ancestry → 色 (gnomAD配色を参考)
ancestry_colors = {
    "afr": "#1f77b4",
    "ami": "#9467bd",
    "amr": "#d62728",
    "asj": "#ff7f0e",
    "eas": "#2ca02c",
    "fin": "#bcbd22",
    "mid": "#e377c2",
    "nfe": "#7f7f7f",
    "oth": "#8c564b",
    "sas": "#17becf",
}

# ===== 散布図描画関数 =====
def scatter_pc(ax, pcx, pcy, df, title):
    # 高確信度サンプル
    high = df[df["max_prob"] >= 0.9]
    low = df[df["max_prob"] < 0.9]

    # 高確信度: 各ancestryで色分け
    for anc, grp in high.groupby("pred_ancestry"):
        ax.scatter(grp[pcx], grp[pcy],
                   c=ancestry_colors.get(anc, "gray"),
                   label=f"{anc} (n={len(grp)})",
                   s=40, alpha=0.7, edgecolors="white", linewidths=0.4)
    # 低確信度: 黒の×マーカ
    if len(low) > 0:
        ax.scatter(low[pcx], low[pcy],
                   c="black", marker="x",
                   label=f"max_prob<0.9 (n={len(low)})",
                   s=60, alpha=0.9, linewidths=1.5)
    ax.set_xlabel(pcx)
    ax.set_ylabel(pcy)
    ax.set_title(title)
    ax.legend(fontsize=8, loc="best")
    ax.grid(True, alpha=0.3)

# Figure 1: PC1-PC2, PC2-PC3, PC3-PC4
fig, axes = plt.subplots(1, 3, figsize=(21, 6))
scatter_pc(axes[0], "PC1", "PC2", df, "PC1 vs PC2")
scatter_pc(axes[1], "PC2", "PC3", df, "PC2 vs PC3")
scatter_pc(axes[2], "PC3", "PC4", df, "PC3 vs PC4")
plt.suptitle(f"Cancer WGS PCA Projection onto gnomAD v3.1 (n={len(df)})", fontsize=14)
plt.tight_layout()
out_png1 = OUT_DIR / "pca_plot_PC1234_05232026.png"
plt.savefig(out_png1, dpi=150, bbox_inches="tight")
plt.close()
log(f"Saved: {out_png1}")

# Figure 2: ancestry分布バー
fig, ax = plt.subplots(figsize=(10, 6))
counts = df["pred_ancestry"].value_counts().sort_index()
colors = [ancestry_colors.get(a, "gray") for a in counts.index]
bars = ax.bar(counts.index, counts.values, color=colors, edgecolor="black")
for bar, v in zip(bars, counts.values):
    ax.text(bar.get_x() + bar.get_width()/2, v + 1, f"{v}",
            ha="center", va="bottom", fontsize=10)
ax.set_xlabel("Predicted Ancestry")
ax.set_ylabel("Number of samples")
ax.set_title(f"Ancestry Distribution (n={len(df)})")
ax.grid(True, axis="y", alpha=0.3)
plt.tight_layout()
out_png2 = OUT_DIR / "ancestry_distribution_05232026.png"
plt.savefig(out_png2, dpi=150, bbox_inches="tight")
plt.close()
log(f"Saved: {out_png2}")

# Figure 3: max_prob ヒストグラム
fig, ax = plt.subplots(figsize=(8, 5))
ax.hist(df["max_prob"], bins=40, color="steelblue", edgecolor="white")
ax.axvline(0.9, color="red", linestyle="--", label="threshold=0.9")
ax.set_xlabel("Max ancestry probability")
ax.set_ylabel("Sample count")
ax.set_title(f"Max ancestry probability distribution (n={len(df)})")
ax.legend()
ax.grid(True, alpha=0.3)
plt.tight_layout()
out_png3 = OUT_DIR / "max_prob_distribution_05232026.png"
plt.savefig(out_png3, dpi=150, bbox_inches="tight")
plt.close()
log(f"Saved: {out_png3}")

# Figure 4: PC heatmap (各PCがancestryをどれだけ分離するか)
fig, axes = plt.subplots(4, 4, figsize=(16, 16))
for i in range(16):
    ax = axes[i // 4][i % 4]
    pc = f"PC{i+1}"
    for anc, grp in df.groupby("pred_ancestry"):
        ax.hist(grp[pc], bins=30, alpha=0.5,
                color=ancestry_colors.get(anc, "gray"), label=anc, density=True)
    ax.set_xlabel(pc)
    ax.set_title(pc, fontsize=10)
    if i == 0:
        ax.legend(fontsize=7, loc="best")
plt.suptitle("PC distributions by predicted ancestry", fontsize=14)
plt.tight_layout()
out_png4 = OUT_DIR / "pc_histograms_by_ancestry_05232026.png"
plt.savefig(out_png4, dpi=120, bbox_inches="tight")
plt.close()
log(f"Saved: {out_png4}")

# ===== サマリ =====
log("")
log("=== Summary ===")
log(f"Total samples: {len(df)}")
log("Ancestry distribution:")
for k, v in df["pred_ancestry"].value_counts().items():
    log(f"  {k}: {v} ({100*v/len(df):.1f}%)")
log(f"Low confidence samples (max_prob<0.9): {(df['max_prob']<0.9).sum()}")

elapsed = time.time() - start_time
log(f"=== Total elapsed time: {elapsed:.1f} sec ===")

with open(LOG_FILE, "w") as f:
    f.write("\n".join(log_lines) + "\n")
log(f"Log saved: {LOG_FILE}")
