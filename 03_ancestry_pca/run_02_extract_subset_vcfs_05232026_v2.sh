#!/bin/bash
# ========================================================================
# Script: run_02_extract_subset_vcfs_05232026.sh
# Pipeline 4 (ancestry / PCA), step 2.
# 処理内容:
#   - SGE array job で各サンプルの normal DeepVariant VCF から
#     gnomAD loadingsの76,399 SNP位置のみを抽出
#   - bcftools view -R sites.bed -f PASS,RefCall でPASS/RefCall両方取得
#       (RefCallはhomozygous-ref相当、欠損ジェノタイプを最小化)
#   - bcftools view -m2 -M2 -v snps で bi-allelic SNV のみ
#   - bcftools annotate --set-id '%CHROM:%POS:%REF:%FIRST_ALT' でID統一
#   - bcftools reheader でサンプル名を UUID から sample_id に変更
#   - 結果: ${WORK_DIR}/vcf_subset/{sample_id}.subset.vcf.gz
#   - 動的列指定: manifest のヘッダから列位置を取得
#   - 実行時間を計測
#
# De-identification note:
#   パスはプレースホルダ ("/path/to/...") です。実行前に編集してください。
#   array サイズ (-t) は実行時の manifest 行数を反映します（最終解析コホートと
#   一致しない場合があります。README の Pipeline 4 注記を参照）。
#   参加者個別データは本リポジトリに含まれません。
# ========================================================================

#$ -S /bin/bash
#$ -cwd

#$ -l s_vmem=4G,mem_req=4G,d_rt=3600,s_rt=3600

#$ -t 1-248
#$ -tc 30
#$ -j y
#$ -o /path/to/PCA_project/logs/
#$ -N pca_extract

set -e
set -o pipefail

START_SEC=$(date +%s)

# ===== Path（編集してください） =====
WORK_DIR=/path/to/PCA_project
SAMPLEINFO=/path/to/sample_manifest.txt
BED=${WORK_DIR}/sitelist/sites.bed
VCF_BASE=/path/to/wgs_results
OUT_DIR=${WORK_DIR}/vcf_subset
LOG_DIR=${WORK_DIR}/logs

mkdir -p ${OUT_DIR}

# bcftools / tabix (user-local htslib/bcftools 1.19; パスは編集してください)
export LD_LIBRARY_PATH=/path/to/htslib-1.19:${LD_LIBRARY_PATH}
BCFTOOLS=/path/to/bcftools
TABIX=/path/to/htslib-1.19/tabix
BGZIP=/path/to/htslib-1.19/bgzip

# ===== 動的列指定: ヘッダから列番号取得 =====
# header: sampleID sex age tumor_sample_name tumor_uuid normal_sample_name normal_uuid analysis_status running_status shared_date bamFile
HEADER=$(head -1 ${SAMPLEINFO})

# スペース区切りで列名→列番号 (1-based for awk)
col_sampleID=$(echo "$HEADER" | tr '\t' '\n' | awk '{print NR, $0}' | awk '$2=="sampleID"{print $1}')
col_tumor_name=$(echo "$HEADER" | tr '\t' '\n' | awk '{print NR, $0}' | awk '$2=="tumor_sample_name"{print $1}')
col_normal_uuid=$(echo "$HEADER" | tr '\t' '\n' | awk '{print NR, $0}' | awk '$2=="normal_uuid"{print $1}')
col_status=$(echo "$HEADER" | tr '\t' '\n' | awk '{print NR, $0}' | awk '$2=="analysis_status"{print $1}')

echo "[$(date)] column positions: sampleID=$col_sampleID, tumor_name=$col_tumor_name, normal_uuid=$col_normal_uuid, status=$col_status"

# ===== analysis_status==1 の行を抽出してアレイ ID 番目をとる =====
LINE=$(awk -F'\t' -v s=$col_status -v t=$SGE_TASK_ID 'NR>1 && $s==1 {n++; if(n==t){print; exit}}' ${SAMPLEINFO})

if [ -z "$LINE" ]; then
    echo "ERROR: no row for SGE_TASK_ID=$SGE_TASK_ID"
    exit 1
fi

SAMPLE_ID=$(echo "$LINE" | awk -F'\t' -v c=$col_sampleID '{print $c}')
TUMOR_NAME=$(echo "$LINE" | awk -F'\t' -v c=$col_tumor_name '{print $c}')
NORMAL_UUID=$(echo "$LINE" | awk -F'\t' -v c=$col_normal_uuid '{print $c}')

echo "[$(date)] Task=$SGE_TASK_ID sample=$SAMPLE_ID tumor_name=$TUMOR_NAME normal_uuid=$NORMAL_UUID"

# ===== 入力VCFパス =====
IN_VCF=${VCF_BASE}/${TUMOR_NAME}/deepvariant/${NORMAL_UUID}/${NORMAL_UUID}.deepvariant.vcf.gz
OUT_VCF=${OUT_DIR}/${SAMPLE_ID}.subset.vcf.gz

if [ ! -f "$IN_VCF" ]; then
    echo "ERROR: input VCF not found: $IN_VCF"
    exit 2
fi

# tbi インデックス確認 (なければ作成)
if [ ! -f "${IN_VCF}.tbi" ]; then
    echo "[$(date)] tbi missing, creating..."
    ${TABIX} -p vcf ${IN_VCF}
fi

# ===== bcftools pipeline =====
# 1. view -R で region 抽出 (PASS/RefCall両方OK)
# 2. view -m2 -M2 -v snps で bi-allelic SNV のみ
# 3. annotate でID統一
# 4. reheader でサンプル名を sample_id に
TMP_VCF=${OUT_DIR}/${SAMPLE_ID}.tmp.vcf.gz
SAMPLE_RENAME=${OUT_DIR}/${SAMPLE_ID}.rename.txt
echo "${SAMPLE_ID}" > ${SAMPLE_RENAME}

echo "[$(date)] Step 1: bcftools view -R BED + bi-allelic SNV filter"
${BCFTOOLS} view \
    -R ${BED} \
    -f PASS,RefCall \
    -m2 -M2 -v snps \
    ${IN_VCF} \
    -O u 2>> ${LOG_DIR}/02_extract_${SAMPLE_ID}.log \
| ${BCFTOOLS} annotate --set-id '%CHROM\:%POS\:%REF\:%FIRST_ALT' \
    -O u 2>> ${LOG_DIR}/02_extract_${SAMPLE_ID}.log \
| ${BCFTOOLS} reheader -s ${SAMPLE_RENAME} 2>> ${LOG_DIR}/02_extract_${SAMPLE_ID}.log \
> ${TMP_VCF}.bcf

echo "[$(date)] Step 2: convert to vcf.gz + tabix"
${BCFTOOLS} view ${TMP_VCF}.bcf -O z -o ${OUT_VCF} 2>> ${LOG_DIR}/02_extract_${SAMPLE_ID}.log
${TABIX} -f -p vcf ${OUT_VCF}

# 中間ファイル削除
rm -f ${TMP_VCF}.bcf ${SAMPLE_RENAME}

# 結果サマリ
N_VARIANTS=$(${BCFTOOLS} view -H ${OUT_VCF} | wc -l)
echo "[$(date)] Sample=${SAMPLE_ID} extracted ${N_VARIANTS} variants"

# ===== 実行時間 =====
END_SEC=$(date +%s)
ELAPSED=$((END_SEC - START_SEC))
echo "[$(date)] Sample=${SAMPLE_ID} done. Elapsed: ${ELAPSED} sec"

exit 0
