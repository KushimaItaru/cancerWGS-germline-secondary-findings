#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# filter_cnvs_v3.py
# Pipeline 2 (CNV), step 4: ACMG class 4/5 selection + positional filtering.
#   This is the step that yields the high-priority CNV workload reported in the
#   manuscript (AnnotSV ACMG class 4/5, autosomal, segmental-duplication overlap
#   < 70%).
# - AnnotSV アノテーション済みCNVをフィルタリング
# - v2 -> v3 の変更点:
#   * ACMG_class による絞り込みを新規追加:
#     - ACMG_class が 4(likely pathogenic) または 5(pathogenic) のCNVのみ通過
#     - 残すクラスは --acmg-classes で変更可能(既定 "4,5")
#     - 対象列は列名(既定 "ACMG_class")で動的に取得(列ズレ防止)。--acmg-column で変更可
#     - 欠損(NA/空/"."/数値化不可)は非該当として除外
#   * 計算効率化: 先にACMG_classをベクトル化で絞り込み、性染色体/Segdup判定は該当行のみ実施
#   * 統計レポートに「ACMG_class非該当で除外」を追加
# - フィルタリング処理(3条件すべてを満たすCNVのみ残す):
#   1. ACMG_class が 4 または 5                      ★v3で追加
#   2. 性染色体 (X, Y) のCNVを除外
#   3. Segmental duplication領域と70%以上重なるCNVを除外
#   AnnotSV の全列を保持したまま出力
# - 実行時間を計測してログに記録
#
# De-identification note:
#   入出力は --input / --segdup / --output で指定します（ハードコードされた個人パスはありません）。
#   参加者個別データは本リポジトリに含まれません（README.md 参照）。

import argparse
import os
import sys
import time
from datetime import datetime

# 環境準備済み前提でインポート
try:
    import pandas as pd
    from intervaltree import IntervalTree, Interval
except ImportError as e:
    print(f"[ERROR] 必要なパッケージが見つかりません: {e}", file=sys.stderr)
    print("[INFO] 以下を実行してください:", file=sys.stderr)
    print("  pip install pandas intervaltree", file=sys.stderr)
    sys.exit(1)


def now_str() -> str:
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def parse_arguments():
    parser = argparse.ArgumentParser(description="CNVをフィルタリング(ACMG_class 4/5・性染色体・Segdup領域)")
    parser.add_argument("--input", required=True,
                        help="AnnotSV出力TSVファイル(all_samples_cnv_annotated.tsv)")
    parser.add_argument("--segdup", required=True,
                        help="Segdup領域ファイル(segdup_hg38_merged.txt)")
    parser.add_argument("--output", required=True,
                        help="フィルタ済み出力ファイル")
    parser.add_argument("--overlap-threshold", type=float, default=70.0,
                        help="Segdup重複率の閾値(%%) (デフォルト: 70)")
    parser.add_argument("--exclude-sex-chr", action="store_true", default=True,
                        help="性染色体(X,Y)のCNVを除外(デフォルト: 有効)")
    parser.add_argument("--acmg-classes", default="4,5",
                        help="残すACMG_classをカンマ区切りで指定(デフォルト: 4,5)")
    parser.add_argument("--acmg-column", default="ACMG_class",
                        help="ACMG分類の列名(デフォルト: ACMG_class)")
    return parser.parse_args()


def parse_acmg_classes(spec: str):
    """'4,5' のような指定を整数集合に変換"""
    allowed = set()
    for x in str(spec).split(','):
        x = x.strip()
        if x == "":
            continue
        try:
            allowed.add(int(x))
        except ValueError:
            print(f"[WARN] ACMG_class指定を無視(数値化不可): {x}", file=sys.stderr)
    return allowed


def acmg_pass(value, allowed_classes) -> bool:
    """ACMG_class が許可クラス(既定:{4,5})に含まれれば True
    欠損(NA/空/'.'/数値化不可)は False(除外)。'4'や'4.0'などの表記揺れにも対応。"""
    s = str(value).strip()
    if s in ("", "NA", "nan", "NaN", "None", "."):
        return False
    try:
        return int(float(s)) in allowed_classes
    except ValueError:
        return False


def load_segdup_regions(segdup_file: str) -> dict:
    """Segmental duplication領域を IntervalTree にロード"""
    chr_trees = {}
    print(f"[INFO] Segdup領域を読み込み: {segdup_file}")

    line_count = 0
    with open(segdup_file, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 3:
                chr_name = parts[0].replace('chr', '')
                try:
                    start = int(parts[1])
                    end = int(parts[2])
                except ValueError:
                    continue

                if chr_name not in chr_trees:
                    chr_trees[chr_name] = IntervalTree()
                chr_trees[chr_name].add(Interval(start, end))
                line_count += 1

    print(f"[INFO] Segdup領域: {len(chr_trees)} 染色体, {line_count} 区間")
    return chr_trees


def calculate_overlap_percentage(cnv_start: int, cnv_end: int, segdup_tree: IntervalTree) -> float:
    """CNVとsegdup領域の重複率(%)を計算"""
    if not segdup_tree:
        return 0.0

    cnv_length = cnv_end - cnv_start
    if cnv_length == 0:
        return 0.0

    overlapping = segdup_tree.overlap(cnv_start, cnv_end)
    if not overlapping:
        return 0.0

    # マージして重複領域の合計を計算
    merged_intervals = []
    for interval in sorted(overlapping):
        overlap_start = max(cnv_start, interval.begin)
        overlap_end = min(cnv_end, interval.end)
        if merged_intervals and overlap_start <= merged_intervals[-1][1]:
            merged_intervals[-1] = (merged_intervals[-1][0], max(merged_intervals[-1][1], overlap_end))
        else:
            merged_intervals.append((overlap_start, overlap_end))

    overlap_length = sum(end - start for start, end in merged_intervals)
    return (overlap_length / cnv_length) * 100.0


def filter_cnvs(input_file: str, segdup_file: str, output_file: str,
                overlap_threshold: float = 70.0, exclude_sex_chr: bool = True,
                acmg_column: str = "ACMG_class", allowed_classes=frozenset({4, 5})) -> None:
    """CNVをフィルタリング"""

    print(f"\n[INFO] 入力ファイル読み込み: {input_file}")
    df = pd.read_csv(input_file, sep='\t', low_memory=False)
    initial_count = len(df)
    print(f"[INFO] 初期CNV数: {initial_count:,}")
    print(f"[INFO] 列名(先頭10): {list(df.columns[:10])}")

    # ACMG_class 列の存在確認(列名で動的取得)
    if acmg_column not in df.columns:
        print(f"[ERROR] ACMG列 '{acmg_column}' が見つかりません。", file=sys.stderr)
        print(f"[INFO] 利用可能な列(一部): {list(df.columns)}", file=sys.stderr)
        sys.exit(1)

    # --- ACMG_class 絞り込み (v3で追加。ベクトル化で高速に該当行のみ抽出) ---
    print(f"\n[INFO] ACMG_class フィルター: 残すクラス = {sorted(allowed_classes)} (列: {acmg_column})")
    acmg_mask = df[acmg_column].apply(lambda v: acmg_pass(v, allowed_classes))
    acmg_excluded_count = int((~acmg_mask).sum())
    df = df[acmg_mask].reset_index(drop=True)
    print(f"[INFO] ACMG_class 該当: {len(df):,} 件 (非該当で除外: {acmg_excluded_count:,})")

    # Segdup領域ロード
    segdup_trees = load_segdup_regions(segdup_file)

    # フィルタリング(ACMG該当行に対してのみ 性染色体/Segdup を判定)
    sex_chr_count = 0
    segdup_overlap_count = 0
    filtered_indices = []

    print(f"\n[INFO] フィルタリング中...")
    print(f"  性染色体除外: {'ON' if exclude_sex_chr else 'OFF'}")
    print(f"  Segdup重複率閾値: {overlap_threshold}%")

    total_rows = len(df)
    progress_interval = max(1, total_rows // 10)

    for idx in range(total_rows):
        if idx > 0 and idx % progress_interval == 0:
            print(f"  処理中... {int(idx / total_rows * 100)}%")

        row = df.iloc[idx]

        # 性染色体チェック
        chr_name = str(row['SV_chrom'])
        if exclude_sex_chr and chr_name in ['X', 'Y', 'chrX', 'chrY']:
            sex_chr_count += 1
            continue

        # Segdup重複チェック
        chr_name_clean = chr_name.replace('chr', '')
        try:
            cnv_start = int(row['SV_start'])
            cnv_end = int(row['SV_end'])
        except (ValueError, KeyError):
            continue

        if chr_name_clean in segdup_trees:
            overlap_pct = calculate_overlap_percentage(
                cnv_start, cnv_end, segdup_trees[chr_name_clean]
            )
            if overlap_pct >= overlap_threshold:
                segdup_overlap_count += 1
                continue

        filtered_indices.append(idx)

    print(f"  処理中... 100%")

    # フィルタ適用
    df_filtered = df.iloc[filtered_indices]

    # 保存
    print(f"\n[INFO] 結果を保存中: {output_file}")
    df_filtered.to_csv(output_file, sep='\t', index=False)

    # 統計
    print()
    print("=" * 50)
    print("  フィルタリング結果")
    print("=" * 50)
    print(f"初期CNV数:                          {initial_count:,}")
    print(f"ACMG_class非該当({sorted(allowed_classes)}以外)で除外: {acmg_excluded_count:,}")
    print(f"性染色体CNVで除外:                  {sex_chr_count:,}")
    print(f"Segdup重複(>={overlap_threshold}%)で除外:         {segdup_overlap_count:,}")
    print(f"除外されたCNV総数:                  {acmg_excluded_count + sex_chr_count + segdup_overlap_count:,}")
    print(f"残ったCNV数:                        {len(df_filtered):,}")
    if initial_count > 0:
        print(f"保持率:                             {(len(df_filtered) / initial_count * 100):.2f}%")
    print("=" * 50)
    print(f"\n出力ファイル: {output_file}")


def main():
    overall_start = time.time()
    print("============================================")
    print("  filter_cnvs_v3.py")
    print("============================================")
    print(f"開始時刻: {now_str()}")

    args = parse_arguments()
    allowed_classes = parse_acmg_classes(args.acmg_classes)
    if not allowed_classes:
        print("[ERROR] --acmg-classes が空です(例: --acmg-classes 4,5)", file=sys.stderr)
        sys.exit(1)

    print(f"入力: {args.input}")
    print(f"Segdup: {args.segdup}")
    print(f"出力: {args.output}")
    print(f"Segdup閾値: {args.overlap_threshold}%")
    print(f"ACMG_class 残すクラス: {sorted(allowed_classes)} (列: {args.acmg_column})")

    # ファイル存在確認
    for label, path in [("入力", args.input), ("Segdup", args.segdup)]:
        if not os.path.exists(path):
            print(f"[ERROR] {label}ファイルが見つかりません: {path}")
            sys.exit(1)

    # 出力ディレクトリ作成
    out_dir = os.path.dirname(args.output)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    filter_cnvs(
        input_file=args.input,
        segdup_file=args.segdup,
        output_file=args.output,
        overlap_threshold=args.overlap_threshold,
        exclude_sex_chr=args.exclude_sex_chr,
        acmg_column=args.acmg_column,
        allowed_classes=allowed_classes,
    )

    overall_end = time.time()
    elapsed = overall_end - overall_start
    print()
    print(f"[INFO] 終了時刻: {now_str()}")
    print(f"[INFO] 処理時間: {elapsed:.2f}秒 ({elapsed/60:.2f}分)")


if __name__ == "__main__":
    main()
