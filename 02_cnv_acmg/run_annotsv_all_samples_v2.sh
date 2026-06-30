#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -l s_vmem=32G,mem_req=32G
#$ -l d_rt=86400
#$ -l s_rt=86400
#$ -N AnnotSV_all_samples
#$ -o logs/annotsv_$JOB_ID.log
#$ -e logs/annotsv_$JOB_ID.err

# run_annotsv_all_samples_v2.sh
# Pipeline 2 (CNV), step 3: annotate the merged CNV BED with AnnotSV (GRCh38),
#   which assigns the ACMG/ClinGen SV classification used by the next step
#   (filter_cnvs_v3.py keeps ACMG class 4/5).
# - 統合BED(all_samples_cnv.bed)にAnnotSVでアノテーションを付与
# - v1 -> v2 の変更点:
#   * パスを現作業ディレクトリ基準($(pwd))に変更
#   * ログファイル名に JOB_ID を含める形式
#   * 実行時間計測を追加
#   * ウォールタイム上限を明示指定 d_rt=86400(24時間)
# - 処理内容:
#   1. bedtools, samtools モジュール環境準備
#   2. AnnotSV を統合BEDに対して実行
#   3. アノテーション済みTSVを出力
#
# De-identification note:
#   ANNOTSV のパスはプレースホルダ ("/path/to/AnnotSV") です。実行前に編集してください。
#   入出力は作業ディレクトリ($(pwd))基準です。
#   参加者個別データは本リポジトリに含まれません(README.md 参照)。

START_TIME=$(date +%s)

mkdir -p logs

echo "============================================"
echo "  AnnotSV all samples (v2)"
echo "============================================"
echo "開始時刻: $(date)"
echo "JOB_ID: ${JOB_ID}"
echo "ホスト名: $(hostname)"
echo "作業ディレクトリ: $(pwd)"
echo ""

# モジュール環境設定
source /etc/profile.d/modules.sh
module use /usr/local/package/modulefiles/
module load bedtools/2.31.1
module load samtools/1.19

# モジュール確認
echo "[INFO] ロード済みモジュール:"
module list 2>&1
echo ""
echo "[INFO] ツールのパス:"
echo "  bedtools: $(which bedtools)"
echo "  samtools: $(which samtools)"
echo ""

# AnnotSV のパスを設定(インストール済みディレクトリを参照。環境に合わせて編集してください)
# 注: 複数の解析で同じ AnnotSV を共有可能
export ANNOTSV=/path/to/AnnotSV

if [ ! -d "${ANNOTSV}" ]; then
    echo "ERROR: AnnotSV ディレクトリが見つかりません: ${ANNOTSV}"
    echo "       AnnotSV をインストール済みかご確認ください(install_annotsv_local.sh)"
    exit 1
fi

if [ ! -x "${ANNOTSV}/bin/AnnotSV" ]; then
    echo "ERROR: AnnotSV 実行ファイルが見つかりません: ${ANNOTSV}/bin/AnnotSV"
    exit 1
fi

echo "[INFO] ANNOTSV: ${ANNOTSV}"

# 入出力パス(現ディレクトリ基準)
INPUT_FILE="$(pwd)/cnv_allSample/all_samples_cnv.bed"
OUTPUT_FILE="$(pwd)/cnv_allSample/all_samples_cnv_annotated.tsv"

# 入力ファイル確認
if [ ! -f "${INPUT_FILE}" ]; then
    echo "ERROR: 入力BEDファイルが見つかりません: ${INPUT_FILE}"
    echo "       先に convert_cnvkit_to_bed_v2.py を実行してください"
    exit 1
fi

echo ""
echo "[INFO] 入力ファイル: ${INPUT_FILE}"
echo "[INFO] 入力ファイル行数: $(wc -l < ${INPUT_FILE})"
echo "[INFO] 出力ファイル: ${OUTPUT_FILE}"
echo ""

# AnnotSV実行
echo "[INFO] AnnotSV を実行..."
echo "============================================"

"${ANNOTSV}/bin/AnnotSV" \
    -SVinputFile "${INPUT_FILE}" \
    -SVinputInfo 1 \
    -outputFile "${OUTPUT_FILE}" \
    -genomeBuild GRCh38 \
    -svtBEDcol 4 \
    -annotationMode full

EXIT_CODE=$?

echo "============================================"
echo ""

if [ ${EXIT_CODE} -eq 0 ]; then
    echo "[INFO] AnnotSV completed successfully"
    if [ -f "${OUTPUT_FILE}" ]; then
        echo "[INFO] 出力ファイル行数: $(wc -l < ${OUTPUT_FILE})"
        echo "[INFO] 出力ファイルサイズ: $(ls -lh ${OUTPUT_FILE} | awk '{print $5}')"
    fi
else
    echo "[ERROR] AnnotSV failed with exit code ${EXIT_CODE}"
fi

END_TIME=$(date +%s)
ELAPSED=$((END_TIME - START_TIME))

echo ""
echo "============================================"
echo "終了時刻: $(date)"
echo "処理時間: ${ELAPSED} 秒 ($((ELAPSED / 60)) 分)"
echo "============================================"

exit ${EXIT_CODE}
