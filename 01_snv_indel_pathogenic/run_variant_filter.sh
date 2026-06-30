#!/bin/bash
# =============================================================================
# run_variant_filter.sh
# Pipeline 1 (SNV/indel): SGE array-job wrapper for variant_filter_parallel_v4.py.
#   Each array task processes one batch of samples (batch_size default = 12).
#   Adjust the array range (-t) to ceil(N_samples / batch_size).
#   This is the exact as-run grid-engine script; resource directives are kept
#   for reproducibility. No participant data is included (see README.md).
# =============================================================================
#$ -S /bin/bash
#$ -cwd
#$ -l s_vmem=24G
#$ -l mem_req=24G
#$ -pe def_slot 1
#$ -o logs/variant_filter_$JOB_ID_$SGE_TASK_ID.log
#$ -e logs/variant_filter_$JOB_ID_$SGE_TASK_ID.err
#$ -t 1-17  # タスクIDは1から始まる（1-17で17バッチ）
#$ -N variant_filter

# ログディレクトリが存在しない場合は作成
mkdir -p logs

# SGE_TASK_IDは1から始まるが、バッチインデックスは0から始まるため調整
BATCH_IDX=$((SGE_TASK_ID - 1))

echo "開始時刻: $(date)"
echo "ホスト名: $(hostname)"
echo "SGE_TASK_ID: $SGE_TASK_ID"
echo "バッチインデックス: $BATCH_IDX"

# バリアントフィルタリング実行
python variant_filter_parallel_v4.py --batch-index $BATCH_IDX

echo "終了時刻: $(date)"
