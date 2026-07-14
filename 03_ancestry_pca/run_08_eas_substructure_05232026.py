#!/usr/bin/env python3
# ========================================================================
# Script: run_08_eas_substructure_05232026.py
# Pipeline 4 (ancestry / PCA), step 8.
# 処理内容:
#   - gnomAD HGDP+1KG metadataのproject_pop=eas samplesのsubpopを取得
#   - 主要EAS sub-popにグルーピング:
#       Japanese (jpt + HGDP japanese)
#       Han_Chinese (chb + chs + han)
#       SE_Asian (khv + cdx + cambodian + dai + lahu + miaozu + she + yizu ...)
#       Northern_Asian (yakut + daur + hezhen + mongola + oroqen + xibo + uygur ...)
#       Other_EAS: その他
#   - EAS reference を sub-pop別に色分け、cancer WGSを重ね描き
#   - cancer WGS samplesのうち、Japanese clusterから外れるサンプルを特定:
#       PC1-PC8空間でJapanese centroidからの距離を計算
#       Mahalanobis距離 > 閾値 のサンプルをサンプルIDラベル付きで強調表示
#   - non-Japanese候補サンプルのリストを TSV 保存
#   - 動的列指定: PC列・subpop列はheader動的取得
#   - 実行時間計測
#
# De-identification note:
#   WORK_DIR などのパスはプレースホルダ ("/path/to/...") です。実行前に編集してください。
#   gnomAD リファレンス（metadata）は公開リソースです。
#   出力 TSV/図には sample_id が含まれ得るため、公開時は出力物（pca/ 配下）を
#   リポジトリに含めないでください（スクリプトのみ公開）。
#   参加者個別データは本リポジトリに含まれません（README.md 参照）。
# ========================================================================

import sys, time, json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.spatial.distance import mahalanobis
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

start_time = time.time()

WORK_DIR = Path("/path/to/PCA_project")
REF_META = WORK_DIR / "reference" / "gnomad_meta_v1.tsv"
ANC_TSV  = WORK_DIR / "pca" / "ancestry_assignments.tsv"
OUT_DIR  = WORK_DIR / "pca"
LOG_FILE = WORK_DIR / "logs" / "08_eas_substructure_05232026.log"

OUT_DIR.mkdir(parents=True, exist_ok=True)
LOG_FILE.parent.mkdir(parents=True, exist_ok=True)

log_lines = []
def log(msg):
    line = f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {msg}"
    print(line)
    log_lines.append(line)

# ===== Step 1: メタデータ読み込み (subpop, project_pop, pca_scores) =====
log(f"Loading metadata: {REF_META}")
header = pd.read_csv(REF_META, sep="\t", nrows=0).columns.tolist()

sample_col = "s"
pop_col = "project_meta.project_pop"
subpop_col = "project_meta.project_subpop"
desc_col = "project_meta.subpop_description"
pca_col = "population_inference.pca_scores"

# Verify all exist
for c in (sample_col, pop_col, subpop_col, desc_col, pca_col):
    if c not in header:
        log(f"WARNING: column '{c}' not in metadata header")

use_cols = [sample_col, pop_col, subpop_col, desc_col, pca_col]
meta = pd.read_csv(REF_META, sep="\t", usecols=use_cols)
log(f"Metadata: {len(meta)} samples")

# Parse pca_scores
def parse_scores(s):
    if pd.isna(s) or s == "NA": return None
    try: return json.loads(s)
    except: return None
parsed = meta[pca_col].apply(parse_scores)
valid = parsed.notna()
log(f"Valid PCA scores: {valid.sum()}")

N_PCS = 16
df_ref = pd.DataFrame({"SAMPLE_ID": meta[sample_col]})
for k in range(N_PCS):
    df_ref[f"PC{k+1}"] = parsed.apply(lambda x: x[k] if x is not None else np.nan)
df_ref["pop"] = meta[pop_col]
df_ref["subpop"] = meta[subpop_col]
df_ref["desc"] = meta[desc_col]
df_ref = df_ref[valid.values].reset_index(drop=True)

# ===== Step 2: EAS sub-pop grouping =====
df_eas = df_ref[df_ref["pop"]=="eas"].copy()
log(f"EAS reference: {len(df_eas)}")

# Sub-pop grouping
def group_subpop(sp):
    sp_l = (sp or "").lower()
    if sp_l in ("jpt", "japanese"):
        return "Japanese"
    if sp_l in ("chb", "chs", "han"):
        return "Han_Chinese"
    if sp_l in ("khv", "cdx", "cambodian", "dai", "lahu", "miaozu", "she", "yizu",
                "naxi", "tujia", "tu"):
        return "SE_Asian"
    if sp_l in ("yakut", "daur", "hezhen", "mongola", "oroqen", "xibo", "uygur"):
        return "Northern_Asian"
    return f"Other_EAS({sp_l})"

df_eas["group"] = df_eas["subpop"].apply(group_subpop)
log(f"EAS group distribution:")
for k, v in df_eas["group"].value_counts().items():
    log(f"  {k}: {v}")

# ===== Step 3: Cancer WGS 読み込み =====
df_can = pd.read_csv(ANC_TSV, sep="\t")
log(f"Cancer WGS: {len(df_can)}")

# ===== Step 4: Japanese centroid からの Mahalanobis distance =====
df_jpt = df_eas[df_eas["group"]=="Japanese"].copy()
log(f"Japanese reference: {len(df_jpt)} (JPT + HGDP Japanese)")

# Use PC1-PC8 for distance (higher PCs are noisier within EAS)
PC_USE = [f"PC{k}" for k in range(1, 9)]
jpt_pc = df_jpt[PC_USE].values
jpt_mean = jpt_pc.mean(axis=0)
jpt_cov = np.cov(jpt_pc.T)
jpt_cov_inv = np.linalg.pinv(jpt_cov)
log(f"Japanese centroid (PC1-PC8): {jpt_mean.round(4)}")

def maha(vec):
    diff = vec - jpt_mean
    return float(np.sqrt(max(diff @ jpt_cov_inv @ diff, 0)))

# Reference EAS Mahalanobis (for setting threshold)
df_eas["maha_jpt"] = df_eas[PC_USE].apply(lambda r: maha(r.values), axis=1)
log(f"\nMahalanobis distance from Japanese centroid (reference EAS by group):")
for k, grp in df_eas.groupby("group"):
    log(f"  {k}: median={grp['maha_jpt'].median():.2f}, max={grp['maha_jpt'].max():.2f}")

# Cancer WGS Mahalanobis
df_can["maha_jpt"] = df_can[PC_USE].apply(lambda r: maha(r.values), axis=1)
log(f"\nCancer WGS Mahalanobis from Japanese centroid:")
log(f"  median: {df_can['maha_jpt'].median():.2f}")
log(f"  90%ile: {df_can['maha_jpt'].quantile(0.9):.2f}")
log(f"  95%ile: {df_can['maha_jpt'].quantile(0.95):.2f}")
log(f"  max:    {df_can['maha_jpt'].max():.2f}")

# Threshold: Japanese reference 99%ile (Chi-square distribution with df=8: 99% = 20.1)
# But for robustness, use empirical Japanese ref 99%ile
jpt_99 = df_jpt[PC_USE].apply(lambda r: maha(r.values), axis=1).quantile(0.99)
log(f"\nJapanese ref 99%ile Mahalanobis: {jpt_99:.2f}")
THRESH = jpt_99
log(f"Threshold for 'non-Japanese candidate': {THRESH:.2f}")

df_can["is_non_japanese_candidate"] = df_can["maha_jpt"] > THRESH
n_candidates = df_can["is_non_japanese_candidate"].sum()
log(f"\nNon-Japanese candidates among cancer WGS cohort: {n_candidates}")

# ===== Step 5: 最も近いEAS sub-popを推定 (kNN) =====
# For each candidate, find nearest reference sub-pop (in PC1-PC8 space)
from scipy.spatial import cKDTree
eas_pc = df_eas[PC_USE].values
tree = cKDTree(eas_pc)

def nearest_subpop(vec, k=20):
    """Return majority group among k nearest neighbors in EAS reference."""
    dists, idxs = tree.query(vec, k=k)
    nearest_groups = df_eas.iloc[idxs]["group"].values
    # majority vote
    vals, counts = np.unique(nearest_groups, return_counts=True)
    return vals[counts.argmax()], counts.max(), dict(zip(vals, counts))

cand_results = []
for _, row in df_can[df_can["is_non_japanese_candidate"]].iterrows():
    nearest, count, dist = nearest_subpop(row[PC_USE].values, k=20)
    cand_results.append({
        "SAMPLE_ID": row["SAMPLE_ID"],
        "maha_jpt": row["maha_jpt"],
        "nearest_subpop_k20": nearest,
        "k20_votes": count,
        "k20_detail": str(dist),
        **{pc: row[pc] for pc in PC_USE},
    })

cand_df = pd.DataFrame(cand_results).sort_values("maha_jpt", ascending=False)
out_cand = OUT_DIR / "non_japanese_candidates_05232026.tsv"
cand_df.to_csv(out_cand, sep="\t", index=False)
log(f"\nSaved candidate list: {out_cand}")
log(f"\n=== Non-Japanese candidate samples ===")
log(cand_df[["SAMPLE_ID","maha_jpt","nearest_subpop_k20","k20_votes"]].to_string(index=False))

# Also save full ranking
df_can_sorted = df_can.sort_values("maha_jpt", ascending=False)
out_rank = OUT_DIR / "cancer_wgs_mahalanobis_ranked_05232026.tsv"
df_can_sorted[["SAMPLE_ID","maha_jpt","is_non_japanese_candidate"] + PC_USE].to_csv(out_rank, sep="\t", index=False)
log(f"Full ranking saved: {out_rank}")

# ===== Step 6: 可視化 =====
group_colors = {
    "Japanese":      "#d62728",   # red
    "Han_Chinese":   "#2ca02c",   # green
    "SE_Asian":      "#ff7f0e",   # orange
    "Northern_Asian":"#9467bd",   # purple
    "Other_EAS":     "#7f7f7f",
}

def plot_eas(ax, pcx, pcy, ref_df, can_df, candidates, title, with_labels=False):
    # Reference EAS by group
    for grp, sub in ref_df.groupby("group"):
        # If group name contains parentheses, use Other_EAS color
        col = group_colors.get(grp, group_colors["Other_EAS"])
        ax.scatter(sub[pcx], sub[pcy], c=col, s=18, alpha=0.55,
                   label=f"{grp} (n={len(sub)})", edgecolors="none")
    # Cancer WGS - main (Japanese-like)
    main = can_df[~can_df["is_non_japanese_candidate"]]
    ax.scatter(main[pcx], main[pcy], c="black", marker="x", s=22, alpha=0.5,
               linewidths=0.7,
               label=f"Cancer WGS Japanese-like (n={len(main)})")
    # Cancer WGS - candidates
    if len(candidates) > 0:
        ax.scatter(candidates[pcx], candidates[pcy], c="#e377c2", marker="D",
                   s=85, edgecolors="black", linewidths=1.5,
                   label=f"Cancer WGS non-Japanese candidates (n={len(candidates)})")
        if with_labels:
            for _, r in candidates.iterrows():
                ax.annotate(r["SAMPLE_ID"], (r[pcx], r[pcy]),
                            xytext=(6, 6), textcoords="offset points",
                            fontsize=8, color="black",
                            bbox=dict(boxstyle="round,pad=0.15", fc="white",
                                      ec="black", alpha=0.7))
    ax.set_xlabel(pcx); ax.set_ylabel(pcy)
    ax.set_title(title)
    ax.legend(fontsize=8, loc="best", markerscale=1.2)
    ax.grid(True, alpha=0.3)

# group projection — merge Other_EAS(*) into single "Other_EAS"
df_eas_plot = df_eas.copy()
df_eas_plot["group"] = df_eas_plot["group"].apply(
    lambda g: "Other_EAS" if g.startswith("Other_EAS") else g)

candidates = df_can[df_can["is_non_japanese_candidate"]].copy()

# Figure 1: full PC1-PC2-PC3-PC4 panels
fig, axes = plt.subplots(1, 3, figsize=(24, 7.5))
plot_eas(axes[0], "PC1", "PC2", df_eas_plot, df_can, candidates, "PC1 vs PC2 (EAS only)")
plot_eas(axes[1], "PC2", "PC3", df_eas_plot, df_can, candidates, "PC2 vs PC3 (EAS only)")
plot_eas(axes[2], "PC3", "PC4", df_eas_plot, df_can, candidates, "PC3 vs PC4 (EAS only)")
plt.suptitle(
    f"EAS sub-population structure: Cancer WGS (n={len(df_can)}) overlaid on HGDP+1KG EAS reference",
    fontsize=13)
plt.tight_layout()
out1 = OUT_DIR / "pca_eas_substructure_PC1234_05232026.png"
plt.savefig(out1, dpi=140, bbox_inches="tight")
plt.close()
log(f"\nSaved: {out1}")

# Figure 2: zoomed with sample ID labels
fig, axes = plt.subplots(1, 2, figsize=(20, 9))
plot_eas(axes[0], "PC3", "PC4", df_eas_plot, df_can, candidates,
         "PC3 vs PC4 (with non-Japanese candidate IDs)", with_labels=True)
plot_eas(axes[1], "PC4", "PC5", df_eas_plot, df_can, candidates,
         "PC4 vs PC5 (with non-Japanese candidate IDs)", with_labels=True)
plt.suptitle(
    f"Non-Japanese candidates labeled (Mahalanobis > {THRESH:.1f} from Japanese centroid in PC1-PC8)",
    fontsize=13)
plt.tight_layout()
out2 = OUT_DIR / "pca_eas_substructure_labeled_05232026.png"
plt.savefig(out2, dpi=150, bbox_inches="tight")
plt.close()
log(f"Saved: {out2}")

# Figure 3: Mahalanobis distance distribution
fig, ax = plt.subplots(figsize=(11, 6))
for grp, sub in df_eas_plot.groupby("group"):
    col = group_colors.get(grp, group_colors["Other_EAS"])
    ax.hist(sub["maha_jpt"], bins=40, alpha=0.5, color=col, label=f"{grp} (n={len(sub)})",
            density=True)
ax.hist(df_can["maha_jpt"], bins=40, alpha=0.6, color="black",
        label=f"Cancer WGS (n={len(df_can)})", density=True, edgecolor="white")
ax.axvline(THRESH, color="red", linestyle="--", label=f"threshold={THRESH:.1f}")
ax.set_xlabel("Mahalanobis distance from Japanese centroid (PC1-PC8)")
ax.set_ylabel("Density")
ax.set_title("Distribution of Mahalanobis distances from Japanese reference centroid")
ax.set_xlim(0, min(50, ax.get_xlim()[1]))
ax.legend(fontsize=9)
ax.grid(True, alpha=0.3)
plt.tight_layout()
out3 = OUT_DIR / "pca_mahalanobis_jpt_05232026.png"
plt.savefig(out3, dpi=140, bbox_inches="tight")
plt.close()
log(f"Saved: {out3}")

elapsed = time.time() - start_time
log(f"=== Total elapsed: {elapsed:.1f} sec ===")
with open(LOG_FILE, "w") as f:
    f.write("\n".join(log_lines) + "\n")
log(f"Log saved: {LOG_FILE}")
