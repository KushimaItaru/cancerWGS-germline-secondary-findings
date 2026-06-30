#!/bin/bash
# =============================================================================
# merge_results.sh
# Pipeline 1 (SNV/indel): SGE wrapper that runs merge_results.py over all batches.
#   Set WORK_DIR and NUM_BATCHES to match your run. No participant data is
#   included (see README.md).
# =============================================================================
#$ -S /bin/bash
#$ -cwd
#$ -l s_vmem=16G
#$ -l mem_req=16G
#$ -o logs/merge_results.log
#$ -e logs/merge_results.err
#$ -N merge_results

echo "全バッチの結果をマージします..."
echo "開始時刻: $(date)"

# 作業ディレクトリの設定（編集してください）
WORK_DIR="/path/to/SF_analysis"
cd $WORK_DIR || { echo "作業ディレクトリに移動できません"; exit 1; }

# 出力ディレクトリの設定
OUTPUT_DIR="./filtered_variants"

# バッチ数（0-16の17バッチ）
NUM_BATCHES=17

# マージスクリプトの実行
python merge_results.py --output-dir $OUTPUT_DIR --num-batches $NUM_BATCHES

echo "マージ処理が完了しました"
echo "終了時刻: $(date)"
