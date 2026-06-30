#!/bin/bash
# ========================================================================
# Script: run_04_pca_projection_05232026.sh
# Pipeline 4 (ancestry / PCA), step 4.
# 処理内容:
#   - PLINK2 --score を使い、gnomAD v3.1 loadings × cancer WGS genotype で
#     PCスコアを計算 (project方式)
#   - --read-freq で gnomAD のpca_af を allele freq として利用
#     → 'variance-standardize' で (g - 2*af) / sqrt(2*af*(1-af)) 正規化
#     → loading との内積 → PCスコア
#   - cols=+scoresums で各PC合計を出力 (Hail pc_project相当)
#   - 出力: pca/pca_scores.sscore (各サンプルのPC1-PC16)
#   - 実行時間を計測
#
# 数式メモ (Hail pc_project equivalent):
#   PC_k = sum_i [ (g_i - 2*af_i) / sqrt(N * 2*af_i * (1-af_i)) * loading_ik ]
#   PLINK2 --score variance-standardize は (g - 2af) / sqrt(2af(1-af)) なので、
#   結果を sqrt(N_SNP) で割れば Hail pc_project に一致
#
# De-identification note:
#   パスはプレースホルダ ("/path/to/..." / "~/miniconda3") です。実行前に編集してください。
#   参加者個別データは本リポジトリに含まれません（README.md 参照）。
# ========================================================================

#$ -S /bin/bash
#$ -cwd
# (optional) specify your queue if required, e.g.  #$ -q your.q
#$ -l s_vmem=16G,mem_req=16G

#$ -j y
#$ -o /path/to/PCA_project/logs/
#$ -N pca_score

set -e
set -o pipefail

START_SEC=$(date +%s)

# ===== Path（編集してください） =====
WORK_DIR=/path/to/PCA_project
PFILE_PREFIX=${WORK_DIR}/merged/cancer_wgs
LOADINGS=${WORK_DIR}/sitelist/loadings_plink.tsv
AFREQ=${WORK_DIR}/sitelist/pca_af.afreq
OUT_PREFIX=${WORK_DIR}/pca/pca_scores
LOG_DIR=${WORK_DIR}/logs

mkdir -p ${WORK_DIR}/pca

# PLINK2 (conda)
source ~/miniconda3/etc/profile.d/conda.sh
conda activate pca
PLINK2=~/miniconda3/envs/pca/bin/plink2

echo "[$(date)] PLINK2 version: $(${PLINK2} --version | head -1)"

# ===== PLINK2 --score 列指定 =====
# loadings_plink.tsv layout:
#   1: SNPID
#   2: A1 (ALT, effect allele)
#   3: A2 (REF)
#   4-19: PC1..PC16
SCORE_COL_START=4
SCORE_COL_END=$((SCORE_COL_START + 16 - 1))   # = 19

echo "[$(date)] Running PLINK2 --score with PC columns ${SCORE_COL_START}-${SCORE_COL_END}"

${PLINK2} \
    --pfile ${PFILE_PREFIX} \
    --read-freq ${AFREQ} \
    --score ${LOADINGS} 1 2 header-read no-mean-imputation variance-standardize cols=+scoresums,+denom,-fid \
    --score-col-nums ${SCORE_COL_START}-${SCORE_COL_END} \
    --threads 4 \
    --memory 12000 \
    --out ${OUT_PREFIX} 2>&1 | tee ${LOG_DIR}/04_pca_score.log

echo "[$(date)] --score done."
ls -lh ${OUT_PREFIX}* 2>&1

# サマリ
N_USED=$(awk '/--score: / {for(i=1;i<=NF;i++) if($i=="processed,") print $(i-1)}' ${LOG_DIR}/04_pca_score.log | head -1)
echo "[$(date)] SNPs used for scoring (from log): ${N_USED:-N/A}"

# ===== 実行時間 =====
END_SEC=$(date +%s)
ELAPSED=$((END_SEC - START_SEC))
echo "[$(date)] Total elapsed: ${ELAPSED} sec"

exit 0
