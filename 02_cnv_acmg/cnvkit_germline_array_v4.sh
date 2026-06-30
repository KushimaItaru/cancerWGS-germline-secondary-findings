#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -o logs/cnvkit_array_$JOB_ID.$TASK_ID.log
#$ -e logs/cnvkit_array_$JOB_ID.$TASK_ID.err
#$ -l s_vmem=128G
#$ -l mem_req=128G
#$ -pe def_slot 4
#$ -l d_rt=14400
#$ -l s_rt=14400
#$ -tc 20
#$ -N cnvkit_germline_full

# cnvkit_germline_array_v4.sh
# Pipeline 2 (CNV), step 1: per-sample germline CNV calling with CNVkit (WGS mode)
#   against a pooled reference, run as an SGE array job (one task per sample).
#   Steps: batch coverage -> CBS segmentation (threshold 0.1) -> call (ploidy 2).
# - 全サンプルに対して CNVkit を並列実行する(全サンプル再計算版)
# - CNVkit のパラメータ(参照・コンテナ・batch/segment/call)は従来版と同一
# - 前提(作業ディレクトリで prepare_sample_list を実行し sample_list.txt を作成済み)
# - 投入時の指示:
#     cd /path/to/cnvkit_workdir
#     N=$(wc -l < sample_list.txt)
#     qsub -t 1-${N} cnvkit_germline_array_v4.sh
#
# De-identification note:
#   パス(コンテナ/参照/bind mount)はプレースホルダ ("/path/to/...") です。実行前に編集してください。
#   sample_list.txt は 2列TSV(<sample_id>\t<bam_path>)。array サイズは投入時に指定します。
#   参加者個別データは本リポジトリに含まれません(README.md 参照)。

STAGE_START=$(date +%s)

echo "============================================"
echo "=== CNVkit Germline Analysis (full recompute) ==="
echo "============================================"
echo "Task ID: ${SGE_TASK_ID}"
echo "Job ID: ${JOB_ID}"
echo "Host: $(hostname)"
echo "Working Dir: $(pwd)"
echo "Job started at: $(date)"
echo "Allocated slots: ${NSLOTS:-1}"
echo ""

# モジュールロード
module use /usr/local/package/modulefiles/
module load apptainer/

# 基本設定(パスは環境に合わせて編集してください)
CNVKIT_IMAGE="/path/to/cnvkit_latest.sif"
POOLED_REFERENCE="/path/to/cnvkit_reference/pooled_reference_20samples.cnn"
BASE_OUTPUT_DIR="./results"
SAMPLE_LIST="./sample_list.txt"   # 2列TSV: <sample_id>\t<bam_path>
LOG_DIR="./logs"

# 必須ファイル確認
echo "[INFO] 必須ファイル確認..."
for f in "${CNVKIT_IMAGE}" "${POOLED_REFERENCE}" "${SAMPLE_LIST}"; do
    if [ ! -f "${f}" ]; then
        echo "ERROR: 必須ファイルが見つかりません: ${f}"
        exit 1
    fi
    echo "  OK ${f}"
done
echo ""

# ログ・出力ディレクトリ作成
mkdir -p "${LOG_DIR}"
mkdir -p "${BASE_OUTPUT_DIR}"

# 現在のタスクのサンプル情報を取得
IFS=$'\t' read -r SAMPLE_ID SAMPLE_BAM < <(sed -n "${SGE_TASK_ID}p" "${SAMPLE_LIST}")

if [ -z "${SAMPLE_ID}" ] || [ -z "${SAMPLE_BAM}" ]; then
    echo "ERROR: タスクID ${SGE_TASK_ID} に対応するサンプルが見つかりません"
    echo "       サンプルリスト行数: $(wc -l < ${SAMPLE_LIST})"
    exit 1
fi

echo "[INFO] Processing sample: ${SAMPLE_ID}"
echo "[INFO] BAM file: ${SAMPLE_BAM}"

# BAMファイルの存在確認
if [ ! -f "${SAMPLE_BAM}" ]; then
    echo "ERROR: BAM file not found: ${SAMPLE_BAM}"
    BAM_DIR=$(dirname "${SAMPLE_BAM}")
    if [ -d "${BAM_DIR}" ]; then
        echo "DEBUG: ディレクトリ内の BAM ファイル一覧:"
        ls -la "${BAM_DIR}" | grep ".bam" | head -5
    fi
    exit 1
fi

# サンプル固有の出力ディレクトリ
OUTPUT_DIR="${BASE_OUTPUT_DIR}/${SAMPLE_ID}"
mkdir -p "${OUTPUT_DIR}"
cd "${OUTPUT_DIR}" || { echo "ERROR: ${OUTPUT_DIR} に移動できません"; exit 1; }

# サンプル名(BAMファイル名から)
SAMPLE_NAME=$(basename "${SAMPLE_BAM}" .bam)

# ステップ1: カバレッジ計算
STEP1_START=$(date +%s)
echo ""
echo "[$(date)] Step 1: Coverage calculation for ${SAMPLE_ID}..."

apptainer exec --bind /path/to/data:/path/to/data "${CNVKIT_IMAGE}" \
    cnvkit.py batch \
    "${SAMPLE_BAM}" \
    --reference "${POOLED_REFERENCE}" \
    --output-dir . \
    --method wgs \
    --segment-method none \
    --count-reads \
    --drop-low-coverage \
    -p 4

if [ $? -ne 0 ]; then
    echo "ERROR: Step 1 (Coverage calculation) failed for ${SAMPLE_ID}"
    exit 1
fi

STEP1_END=$(date +%s)
echo "[INFO] Step 1 completed in $((STEP1_END - STEP1_START)) seconds"

# ステップ2: CBS閾値0.1でセグメンテーション
STEP2_START=$(date +%s)
echo ""
echo "[$(date)] Step 2: Segmentation with CBS threshold 0.1..."

CNR_FILE="${SAMPLE_NAME}.cnr"
CNS_FILE="${SAMPLE_NAME}_cbs01.cns"

apptainer exec --bind /path/to/data:/path/to/data "${CNVKIT_IMAGE}" \
    cnvkit.py segment \
    "${CNR_FILE}" \
    --method cbs \
    --threshold 0.1 \
    --drop-low-coverage \
    -p 4 \
    -o "${CNS_FILE}"

if [ $? -ne 0 ]; then
    echo "ERROR: Step 2 (Segmentation) failed for ${SAMPLE_ID}"
    exit 1
fi

STEP2_END=$(date +%s)
echo "[INFO] Step 2 completed in $((STEP2_END - STEP2_START)) seconds"

# ステップ3: CNVコール
STEP3_START=$(date +%s)
echo ""
echo "[$(date)] Step 3: Calling CNVs..."

CALL_FILE="${SAMPLE_NAME}_cbs01.call.cns"

apptainer exec --bind /path/to/data:/path/to/data "${CNVKIT_IMAGE}" \
    cnvkit.py call \
    "${CNS_FILE}" \
    --ploidy 2 \
    --thresholds=-1.1,-0.85,0.5,0.7 \
    -o "${CALL_FILE}"

if [ $? -ne 0 ]; then
    echo "ERROR: Step 3 (CNV calling) failed for ${SAMPLE_ID}"
    exit 1
fi

STEP3_END=$(date +%s)
echo "[INFO] Step 3 completed in $((STEP3_END - STEP3_START)) seconds"

# 結果のサマリー
echo ""
echo "============================================"
echo "=== CNV Summary for ${SAMPLE_ID} ==="
echo "============================================"
TOTAL_CNVS=$(awk 'NR>1 && $6!=2' "${CALL_FILE}" | wc -l)
DELETIONS=$(awk 'NR>1 && $6<2' "${CALL_FILE}" | wc -l)
DUPLICATIONS=$(awk 'NR>1 && $6>2' "${CALL_FILE}" | wc -l)
HOMO_DEL=$(awk 'NR>1 && $6==0' "${CALL_FILE}" | wc -l)
HET_DEL=$(awk 'NR>1 && $6==1' "${CALL_FILE}" | wc -l)
DUP=$(awk 'NR>1 && $6==3' "${CALL_FILE}" | wc -l)
AMP=$(awk 'NR>1 && $6>=4' "${CALL_FILE}" | wc -l)
SEGMENTS=$(tail -n +2 "${CNS_FILE}" | wc -l)

echo "Total CNVs detected: ${TOTAL_CNVS}"
echo "  - Homo Deletions (CN=0): ${HOMO_DEL}"
echo "  - Het Deletions (CN=1): ${HET_DEL}"
echo "  - Duplications (CN=3): ${DUP}"
echo "  - Amplifications (CN>=4): ${AMP}"
echo "Segments created: ${SEGMENTS}"

# 成功フラグファイル作成
touch SUCCESS

STAGE_END=$(date +%s)
TOTAL_ELAPSED=$((STAGE_END - STAGE_START))

echo ""
echo "============================================"
echo "=== Job Completion ==="
echo "============================================"
echo "Sample: ${SAMPLE_ID}"
echo "Task ID: ${SGE_TASK_ID}"
echo "Job ID: ${JOB_ID}"
echo "Job completed at: $(date)"
echo "Total processing time: ${TOTAL_ELAPSED} seconds ($((TOTAL_ELAPSED / 60)) min)"
echo "Results saved in: $(pwd)"
echo "============================================"
