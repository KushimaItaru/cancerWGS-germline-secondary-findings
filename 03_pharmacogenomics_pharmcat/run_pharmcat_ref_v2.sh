#!/bin/bash
# =============================================================================
# run_pharmcat_ref_v2.sh
# Pipeline 3 (PharmCAT), step 2: run PharmCAT (v2.15.5) on every preprocessed VCF
#   produced by step 1, generating per-sample report/phenotype JSON + HTML.
#   Pin pharmcat-2.15.5-all.jar; edit OUTPUT_DIR and the jar path. No participant
#   data is included (see README.md).
# =============================================================================

#$ -S /bin/bash
#$ -cwd
#$ -l s_vmem=8G,mem_req=8G
#$ -o pharmcat_ref_fixed.out
#$ -e pharmcat_ref_fixed.err

# 必要なJavaモジュールをロード
module use /usr/local/package/modulefiles/
module load java/17.0.9.8.1

# Javaのメモリ設定を調整
unset JAVA_TOOL_OPTIONS
export _JAVA_OPTIONS="-Xmx4g -Xms2g -XX:ParallelGCThreads=2 -XX:ConcGCThreads=2 -XX:G1ConcRefinementThreads=2 -Djava.util.concurrent.ForkJoinPool.common.parallelism=2"

# 出力ディレクトリを設定（編集してください）
OUTPUT_DIR="/path/to/pharmcat/output"

echo "開始時刻: $(date)"

# 出力ディレクトリ内のbgzファイルを処理
echo "出力ディレクトリ内のbgzファイルを検索しています..."
for VCF_FILE in ${OUTPUT_DIR}/*.preprocessed.vcf.bgz; do
    # ファイル名からサンプル情報を抽出
    FILENAME=$(basename "$VCF_FILE")
    SAMPLE_PREFIX=${FILENAME%.preprocessed.vcf.bgz}  # AAA__BBB形式

    # ファイル名からUUID部分を抽出（サンプル名として使用）
    UUID=$(echo $SAMPLE_PREFIX | cut -d'__' -f2)

    echo "処理中: $SAMPLE_PREFIX"
    echo "使用するサンプル名: $UUID"

    # PharmCATを実行（jar のパスは編集してください）
    echo "PharmCATを実行中..."
    java -jar /path/to/pharmcat-2.15.5-all.jar \
      -vcf "$VCF_FILE" \
      -o ${OUTPUT_DIR} \
      -s "$UUID" \
      -bf "${SAMPLE_PREFIX}"

    # 処理結果の確認
    if [ -f "${OUTPUT_DIR}/${SAMPLE_PREFIX}.report.html" ]; then
        echo "成功: レポートファイルが生成されました"
    else
        echo "警告: レポートファイルが見つかりません"
    fi

    echo "完了: $SAMPLE_PREFIX"
    echo "-------------------------------"
done

echo "全サンプルの処理が完了しました"
echo "終了時刻: $(date)"
