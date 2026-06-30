#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# convert_cnvkit_to_bed_v2.py
# Pipeline 2 (CNV), step 2: CNVkit .call.cns -> AnnotSV BED, merged across samples.
# - CNVkitの.call.cnsファイルを AnnotSV 用 BED 形式に変換し、全サンプルを統合
# - v1 -> v2 の変更点:
#   * 実行時間計測を追加
#   * 引数化(--sample-list, --results-dir, --output-dir)
#   * 進捗表示の改善(N/total 形式)
#   * エラーログをタイムスタンプ付きで詳細化
#   * 統合BED出力前のソート処理を追加(AnnotSV互換性向上)
# - 処理内容:
#   1. sample_list.txt から処理対象サンプルを取得
#   2. 各サンプルの results/<sample>/*_cbs01.call.cns を読み込み
#   3. CN != 2 のCNVのみ抽出
#   4. DEL/DUP の SV_type を判定
#   5. AnnotSV 用 BED 形式(svtype は4列目)に変換
#   6. 全サンプル分を all_samples_cnv.bed に統合
#
# De-identification note:
#   入出力は引数(--sample-list / --results-dir / --output-dir)で指定します
#   （ハードコードされた個人パスはありません）。
#   参加者個別データは本リポジトリに含まれません(README.md 参照)。

import argparse
import glob
import os
import sys
import time
from datetime import datetime
from pathlib import Path

import pandas as pd


def now_str() -> str:
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def parse_arguments():
    parser = argparse.ArgumentParser(description="CNVkitの結果をAnnotSV用BED形式に変換")
    parser.add_argument("--sample-list", default="./sample_list.txt",
                        help="サンプルリストファイル(tab区切り: sample_id, bam_path)")
    parser.add_argument("--results-dir", default="./results",
                        help="CNVkit結果ディレクトリ")
    parser.add_argument("--output-dir", default="./cnv_allSample",
                        help="出力ディレクトリ")
    return parser.parse_args()


def convert_cnvkit_to_bed(cns_file: str, sample_id: str):
    """CNVkitの.call.cnsファイルをAnnotSV用BED形式に変換"""
    try:
        df = pd.read_csv(cns_file, sep='\t')

        # CNVのみ抽出(CN != 2)
        cnv_df = df[df['cn'] != 2].copy()

        if len(cnv_df) == 0:
            return None

        # SVタイプ判定
        cnv_df['svtype'] = cnv_df['cn'].apply(lambda x: 'DEL' if x < 2 else 'DUP')

        # AnnotSV用BED形式(4列目に svtype が来る形式: -svtBEDcol 4 と整合)
        bed_df = pd.DataFrame({
            'chr': cnv_df['chromosome'],
            'start': cnv_df['start'],
            'end': cnv_df['end'],
            'svtype': cnv_df['svtype'],
            'sample_id': sample_id,
            'cn': cnv_df['cn'],
            'log2': cnv_df['log2'],
            'probes': cnv_df['probes'],
            'gene': cnv_df['gene']
        })

        return bed_df

    except Exception as e:
        print(f"  [ERROR] {sample_id}: {e}", file=sys.stderr)
        return None


def main():
    overall_start = time.time()
    print("============================================")
    print("  convert_cnvkit_to_bed_v2.py")
    print("============================================")
    print(f"開始時刻: {now_str()}")

    args = parse_arguments()

    SAMPLE_LIST = args.sample_list
    CNVKIT_RESULTS = args.results_dir
    OUTPUT_DIR = args.output_dir

    # 出力ディレクトリ作成(race condition対応)
    Path(OUTPUT_DIR).mkdir(parents=True, exist_ok=True)

    print(f"サンプルリスト: {SAMPLE_LIST}")
    print(f"結果ディレクトリ: {CNVKIT_RESULTS}")
    print(f"出力ディレクトリ: {OUTPUT_DIR}")
    print()

    # 必須ファイル確認
    if not os.path.exists(SAMPLE_LIST):
        print(f"[ERROR] サンプルリストが見つかりません: {SAMPLE_LIST}")
        sys.exit(1)
    if not os.path.isdir(CNVKIT_RESULTS):
        print(f"[ERROR] 結果ディレクトリが見つかりません: {CNVKIT_RESULTS}")
        sys.exit(1)

    # サンプルリスト読み込み
    sample_ids = []
    with open(SAMPLE_LIST, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split('\t')
            if len(parts) >= 1:
                sample_ids.append(parts[0])

    total_samples = len(sample_ids)
    print(f"処理対象サンプル数: {total_samples}")
    print()

    # 処理統計
    processed_samples = 0
    samples_with_cnvs = 0
    total_cnvs = 0
    no_cnv_samples = []
    failed_samples = []
    all_beds = []

    for i, sample_id in enumerate(sample_ids, start=1):
        if i % 20 == 0 or i == total_samples:
            print(f"[INFO] Processing {i}/{total_samples} samples...")

        cnvkit_dir = os.path.join(CNVKIT_RESULTS, sample_id)
        cns_files = glob.glob(os.path.join(cnvkit_dir, "*cbs01.call.cns"))

        if not cns_files:
            failed_samples.append(sample_id)
            continue

        cns_file = cns_files[0]
        bed_df = convert_cnvkit_to_bed(cns_file, sample_id)

        if bed_df is not None:
            # 個別BED保存
            output_file = os.path.join(OUTPUT_DIR, f"{sample_id}_cnv.bed")
            bed_df.to_csv(output_file, sep='\t', index=False, header=False)
            all_beds.append(bed_df)
            processed_samples += 1
            samples_with_cnvs += 1
            total_cnvs += len(bed_df)
        else:
            no_cnv_samples.append(sample_id)
            processed_samples += 1

    # 全サンプル統合
    if all_beds:
        print()
        print("[INFO] 全サンプルのCNVを統合中...")
        merged_bed = pd.concat(all_beds, ignore_index=True)

        # 染色体・位置でソート(AnnotSV互換性のため)
        # chr1, chr2, ..., chrX, chrY の順序を保つために染色体を文字列ソートしない
        def chr_sort_key(chrom):
            c = str(chrom).replace('chr', '')
            try:
                return (0, int(c))
            except ValueError:
                # X, Y, MT などは後ろへ
                order = {'X': 100, 'Y': 101, 'M': 102, 'MT': 103}
                return (0, order.get(c, 999))

        merged_bed['_sort_key'] = merged_bed['chr'].apply(chr_sort_key)
        merged_bed = merged_bed.sort_values(by=['_sort_key', 'start']).drop(columns=['_sort_key'])

        merged_file = os.path.join(OUTPUT_DIR, "all_samples_cnv.bed")
        merged_bed.to_csv(merged_file, sep='\t', index=False, header=False)
        print(f"[INFO] 統合ファイル保存: {merged_file}")
        print(f"[INFO] 統合BEDの行数: {len(merged_bed):,}")

    # サマリー出力
    print()
    print("============================================")
    print("  Processing Summary")
    print("============================================")
    print(f"Total samples in list: {total_samples}")
    print(f"Successfully processed: {processed_samples}")
    print(f"Samples with CNVs: {samples_with_cnvs}")
    print(f"Samples without CNVs: {len(no_cnv_samples)}")
    print(f"Failed to process: {len(failed_samples)}")
    print(f"Total CNVs detected: {total_cnvs}")
    if samples_with_cnvs > 0:
        print(f"Average CNVs per sample: {total_cnvs/samples_with_cnvs:.1f}")

    # 失敗・CNVなしサンプルをファイル出力
    if no_cnv_samples:
        no_cnv_file = os.path.join(OUTPUT_DIR, "samples_without_cnvs.txt")
        with open(no_cnv_file, 'w') as f:
            for s in no_cnv_samples:
                f.write(f"{s}\n")
        print(f"[INFO] CNVなしサンプル: {no_cnv_file}")

    if failed_samples:
        failed_file = os.path.join(OUTPUT_DIR, "failed_samples.txt")
        with open(failed_file, 'w') as f:
            for s in failed_samples:
                f.write(f"{s}\n")
        print(f"[INFO] 失敗サンプル: {failed_file}")

    # 統計ファイル
    stats_file = os.path.join(OUTPUT_DIR, "conversion_stats.txt")
    with open(stats_file, 'w') as f:
        f.write("=== CNVkit to BED Conversion Statistics ===\n")
        f.write(f"Generated: {now_str()}\n")
        f.write(f"Total samples: {total_samples}\n")
        f.write(f"Processed: {processed_samples}\n")
        f.write(f"With CNVs: {samples_with_cnvs}\n")
        f.write(f"Without CNVs: {len(no_cnv_samples)}\n")
        f.write(f"Failed: {len(failed_samples)}\n")
        f.write(f"Total CNVs: {total_cnvs}\n")
        if samples_with_cnvs > 0:
            f.write(f"Average CNVs per sample: {total_cnvs/samples_with_cnvs:.1f}\n")

    overall_end = time.time()
    elapsed = overall_end - overall_start

    print()
    print(f"[INFO] 統計ファイル: {stats_file}")
    print(f"[INFO] 終了時刻: {now_str()}")
    print(f"[INFO] 処理時間: {elapsed:.2f}秒 ({elapsed/60:.2f}分)")
    print("============================================")


if __name__ == "__main__":
    main()
