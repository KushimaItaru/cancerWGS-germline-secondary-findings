#!/bin/bash
# ========================================================================
# Script: run_03_merge_and_plink_05232026.sh
# Pipeline 4 (ancestry / PCA), step 3.
# 処理内容:
#   - vcf_subset/ 内の subset.vcf.gz を bcftools merge で結合
#   - merge オプション -m none で各SNPの multi-allelic化を避け、bi-allelic維持
#   - 結果VCFを PLINK2 pgen 形式に変換 (--vcf-half-call missing)
#   - SNPサイトリスト (sites.tsv) と extract 結合数の一致を確認
#   - 実行時間を計測
#
# De-identification note:
#   パスはプレースホルダ ("/path/to/..." / "~/miniconda3") です。実行前に編集してください。
#   参加者個別データは本リポジトリに含まれません（README.md 参照）。
# ========================================================================

#$ -S /bin/bash
#$ -cwd

#$ -l s_vmem=32G,mem_req=32G,d_rt=14400,s_rt=14400

#$ -j y
#$ -o /path/to/PCA_project/logs/
#$ -N pca_merge

set -e
set -o pipefail

START_SEC=$(date +%s)

# ===== Path（編集してください） =====
WORK_DIR=/path/to/PCA_project
SUBSET_DIR=${WORK_DIR}/vcf_subset
MERGE_DIR=${WORK_DIR}/merged
LOG_DIR=${WORK_DIR}/logs
SAMPLEINFO=/path/to/sample_manifest.txt

mkdir -p ${MERGE_DIR}

# bcftools 1.19（パスは編集してください）
export LD_LIBRARY_PATH=/path/to/htslib-1.19:${LD_LIBRARY_PATH}
BCFTOOLS=/path/to/bcftools
TABIX=/path/to/htslib-1.19/tabix

# PLINK2 (conda)
source ~/miniconda3/etc/profile.d/conda.sh
conda activate pca
PLINK2=~/miniconda3/envs/pca/bin/plink2

echo "[$(date)] PLINK2: $(${PLINK2} --version | head -1)"
echo "[$(date)] bcftools: $(${BCFTOOLS} --version | head -1)"

# ===== Step 1: subset VCFリスト作成 (動的にanalysis_status==1のサンプルIDを取得) =====
HEADER=$(head -1 ${SAMPLEINFO})
col_sampleID=$(echo "$HEADER" | tr '\t' '\n' | awk '{print NR, $0}' | awk '$2=="sampleID"{print $1}')
col_status=$(echo "$HEADER" | tr '\t' '\n' | awk '{print NR, $0}' | awk '$2=="analysis_status"{print $1}')

VCF_LIST=${MERGE_DIR}/vcf_list.txt
> ${VCF_LIST}
N_FOUND=0
N_MISSING=0
MISSING_SAMPLES=()

while IFS= read -r SAMPLE_ID; do
    VCF=${SUBSET_DIR}/${SAMPLE_ID}.subset.vcf.gz
    if [ -f "$VCF" ] && [ -f "${VCF}.tbi" ]; then
        echo "${VCF}" >> ${VCF_LIST}
        N_FOUND=$((N_FOUND+1))
    else
        N_MISSING=$((N_MISSING+1))
        MISSING_SAMPLES+=("$SAMPLE_ID")
    fi
done < <(awk -F'\t' -v s=$col_status 'NR>1 && $s==1 {print $1}' ${SAMPLEINFO} | sort -u)

echo "[$(date)] subset VCFs found:   $N_FOUND"
echo "[$(date)] subset VCFs missing: $N_MISSING"
if [ $N_MISSING -gt 0 ]; then
    echo "[$(date)] WARNING: missing samples (showing first 10):"
    printf '%s\n' "${MISSING_SAMPLES[@]:0:10}"
fi

# ===== Step 2: bcftools merge =====
MERGED_VCF=${MERGE_DIR}/cancer_wgs_merged.vcf.gz

echo "[$(date)] Step 2: bcftools merge → ${MERGED_VCF}"
${BCFTOOLS} merge \
    -l ${VCF_LIST} \
    -0 \
    -m none \
    --threads 4 \
    -O z \
    -o ${MERGED_VCF} 2>> ${LOG_DIR}/03_merge.log

${TABIX} -f -p vcf ${MERGED_VCF}

N_SITES_MERGED=$(${BCFTOOLS} view -H ${MERGED_VCF} | wc -l)
N_SAMPLES_MERGED=$(${BCFTOOLS} query -l ${MERGED_VCF} | wc -l)
echo "[$(date)] merged VCF: ${N_SITES_MERGED} sites x ${N_SAMPLES_MERGED} samples"

# ===== Step 3: PLINK2 pgen 変換 =====
PFILE_PREFIX=${MERGE_DIR}/cancer_wgs

echo "[$(date)] Step 3: PLINK2 pgen conversion → ${PFILE_PREFIX}"
${PLINK2} \
    --vcf ${MERGED_VCF} \
    --vcf-half-call missing \
    --set-all-var-ids '@:#:$r:$a' \
    --new-id-max-allele-len 10 \
    --make-pgen \
    --threads 4 \
    --memory 16000 \
    --out ${PFILE_PREFIX} 2>&1 | tee ${LOG_DIR}/03_plink2_convert.log

echo "[$(date)] PLINK2 pgen done"
ls -lh ${PFILE_PREFIX}.* 2>&1

# ===== Step 4: サイトリスト整合性確認 =====
N_LOAD_SITES=$(awk 'NR>1' ${WORK_DIR}/sitelist/sites.tsv | wc -l)
N_PFILE_SITES=$(wc -l < ${PFILE_PREFIX}.pvar | awk '{print $1-1}')  # ヘッダ行除く

echo "[$(date)] Site count: loadings=${N_LOAD_SITES}, pgen=${N_PFILE_SITES}"

# ===== 実行時間 =====
END_SEC=$(date +%s)
ELAPSED=$((END_SEC - START_SEC))
echo "[$(date)] Total elapsed: ${ELAPSED} sec ($(awk -v t=$ELAPSED 'BEGIN{printf "%.2f", t/60}') min)"

exit 0
