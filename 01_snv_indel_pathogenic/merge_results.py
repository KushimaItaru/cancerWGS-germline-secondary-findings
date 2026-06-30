#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# =============================================================================
# merge_results.py
# Pipeline 1 (SNV/indel): merge per-batch filtered tables into one cohort-level
#   table of P/LP germline SNV/indel records, plus a merged per-sample statistics
#   table. Builds a compact front-matter column block (SampleID, variant key,
#   gene, SF-gene flag, AAChange, ClinVar VariationID/CLNSIG/review, InterVar,
#   FORMAT) and appends the remaining ANNOVAR columns. Columns are resolved by
#   header NAME (no fixed indices). No participant data is included (see README.md).
# =============================================================================

import os
import sys
import glob
import argparse
import time

def parse_arguments():
    """コマンドライン引数を解析"""
    parser = argparse.ArgumentParser(description='複数バッチの結果を統合')
    parser.add_argument('--output-dir', type=str, default='./filtered_variants',
                        help='出力ディレクトリ')
    parser.add_argument('--num-batches', type=int, required=True,
                        help='バッチ数')
    return parser.parse_args()

def merge_output_files(output_dir, num_batches):
    """複数のバッチ結果を1つに統合"""
    print(f"\n{num_batches}バッチの結果を統合します...")

    # バッチ情報ファイルを確認
    batch_info_files = []
    for i in range(num_batches):
        batch_info_file = os.path.join(output_dir, f"batch_{i}_info.txt")
        if os.path.exists(batch_info_file):
            batch_info_files.append(batch_info_file)
        else:
            print(f"警告: バッチ情報ファイルが見つかりません: {batch_info_file}")

    if not batch_info_files:
        print("統合するバッチファイルがありません。終了します。")
        return

    # 出力ファイルリストを収集
    output_files = []
    for info_file in batch_info_files:
        with open(info_file, 'r') as f:
            for line in f:
                file_path = line.strip()
                if os.path.exists(file_path):
                    output_files.append(file_path)
                else:
                    print(f"警告: 出力ファイルが見つかりません: {file_path}")

    if not output_files:
        print("統合するファイルがありません。終了します。")
        return

    # 日時を含むファイル名を生成
    timestamp = time.strftime("%Y%m%d_%H%M%S")
    merged_output_file = os.path.join(output_dir, f"all_samples_pathogenic_variants_{timestamp}.txt")

    # 最初のファイルのヘッダーを取得し、インデックスを特定
    with open(output_files[0], 'r') as f:
        header = f.readline().strip()
        header_fields = header.split('\t')

        # 必要な列のインデックスを特定
        try:
            sample_id_idx = header_fields.index("SampleID")
            gene_idx = header_fields.index("Gene.refGene")
            sf_genes_idx = header_fields.index("SF_106Genes")
            aa_change_idx = header_fields.index("AAChange.refGene")
            clinvar_id_idx = header_fields.index("clinvarID")
            clinvar_clnsig_idx = header_fields.index("clinvar_CLNSIG")
            review_status_idx = header_fields.index("review_status")
            intervar_idx = header_fields.index("InterVar_automated")
            otherinfo13_idx = header_fields.index("Otherinfo13")
            chr_idx = header_fields.index("Chr")
            start_idx = header_fields.index("Start")
            end_idx = header_fields.index("End")
            ref_idx = header_fields.index("Ref")
            alt_idx = header_fields.index("Alt")
        except ValueError as e:
            print(f"エラー: 必要な列が見つかりません: {e}")
            print(f"ヘッダー行: {header}")
            return

    # 新しいヘッダー行を作成（先頭の特定列）
    new_header = ["SampleID", "variants", "Gene.refGene", "SF_106Genes", "AAChange.refGene",
                  "clinvarID/clinvar_CLNSIG/review_status", "InterVar", "GT:GQ:DP:AD:VAF:PL"]

    # 特別な処理を必要とする列のインデックスをセット
    special_indices = {
        sample_id_idx, gene_idx, sf_genes_idx, aa_change_idx,
        clinvar_id_idx, clinvar_clnsig_idx, review_status_idx,
        intervar_idx, otherinfo13_idx,
        chr_idx, start_idx, end_idx, ref_idx, alt_idx
    }

    # 残りの列を追加（特別処理が必要な列は除外）
    for i, column in enumerate(header_fields):
        if i not in special_indices and column not in new_header:
            new_header.append(column)

    # 統合ファイルを作成
    with open(merged_output_file, 'w') as outfile:
        # ヘッダーを書き込む
        outfile.write('\t'.join(new_header) + '\n')

        # 各ファイルの内容を追加
        total_lines = 0
        for file_path in output_files:
            with open(file_path, 'r') as infile:
                # ヘッダー行をスキップ
                next(infile)

                # 各行を統合ファイルに書き込む
                file_lines = 0
                for line in infile:
                    fields = line.strip().split('\t')

                    # variants列の作成（Chr:Start-End:Ref:Alt）
                    variants_value = f"{fields[chr_idx]}:{fields[start_idx]}-{fields[end_idx]}:{fields[ref_idx]}:{fields[alt_idx]}"

                    # clinvarID/clinvar_CLNSIG/review_statusの作成
                    clinvar_combined = f"{fields[clinvar_id_idx]}/{fields[clinvar_clnsig_idx]}/{fields[review_status_idx]}"

                    # 新しい行を作成
                    new_row = []

                    # 新しいヘッダーの順序に従って値を追加
                    for column in new_header:
                        if column == "SampleID":
                            new_row.append(fields[sample_id_idx])
                        elif column == "variants":
                            new_row.append(variants_value)
                        elif column == "Gene.refGene":
                            new_row.append(fields[gene_idx])
                        elif column == "SF_106Genes":
                            new_row.append(fields[sf_genes_idx])
                        elif column == "AAChange.refGene":
                            new_row.append(fields[aa_change_idx])
                        elif column == "clinvarID/clinvar_CLNSIG/review_status":
                            new_row.append(clinvar_combined)
                        elif column == "InterVar":
                            new_row.append(fields[intervar_idx])
                        elif column == "GT:GQ:DP:AD:VAF:PL":
                            new_row.append(fields[otherinfo13_idx])
                        else:
                            # その他の列（元のヘッダーからインデックスを取得）
                            try:
                                orig_idx = header_fields.index(column)
                                new_row.append(fields[orig_idx])
                            except ValueError:
                                # 元のヘッダーに存在しない列の場合
                                new_row.append(".")

                    # 行を書き込む
                    outfile.write('\t'.join(new_row) + '\n')

                    total_lines += 1
                    file_lines += 1

                print(f"ファイル {file_path}: {file_lines}行")

    print(f"統合が完了しました: {len(output_files)}ファイルから{total_lines}行をマージしました")
    print(f"統合ファイル: {merged_output_file}")

def merge_stats_files(output_dir, num_batches):
    """各バッチの統計情報を統合"""
    print(f"\n各バッチの統計情報を統合します...")

    # 統計ファイルを確認
    stats_files = []
    for i in range(num_batches):
        stats_file = os.path.join(output_dir, f"batch_{i}_stats.txt")
        if os.path.exists(stats_file):
            stats_files.append(stats_file)

    if not stats_files:
        print("統合する統計ファイルがありません。")
        return

    # 日時を含むファイル名を生成
    timestamp = time.strftime("%Y%m%d_%H%M%S")
    merged_stats_file = os.path.join(output_dir, f"all_samples_stats_{timestamp}.txt")

    # 統計情報を読み込む
    all_stats = []
    for stats_file in stats_files:
        with open(stats_file, 'r') as f:
            header = next(f).strip()  # ヘッダー行をスキップ
            for line in f:
                fields = line.strip().split('\t')
                if len(fields) >= 8:  # GQとDPのフィールドが追加されたので8に変更
                    all_stats.append({
                        "sample_id": fields[0],
                        "total_variants": int(fields[1]),
                        "pass_variants": int(fields[2]),
                        "gq_filtered_variants": int(fields[3]),
                        "dp_filtered_variants": int(fields[4]),
                        "rare_variants": int(fields[5]),
                        "output_variants": int(fields[6]),
                        "acmg_yes_variants": int(fields[7])
                    })

    # 統合統計ファイルを作成
    with open(merged_stats_file, 'w') as f:
        f.write("SampleID\t総バリアント数\tPASSバリアント\tGQフィルター\tDPフィルター\t稀なバリアント\t病原性バリアント\tACMG遺伝子\n")

        for stat in all_stats:
            f.write(f"{stat['sample_id']}\t{stat['total_variants']}\t{stat['pass_variants']}\t{stat['gq_filtered_variants']}\t{stat['dp_filtered_variants']}\t{stat['rare_variants']}\t{stat['output_variants']}\t{stat['acmg_yes_variants']}\n")

    # 全体の集計
    total_samples = len(all_stats)
    total_variants = sum(stat['total_variants'] for stat in all_stats)
    total_pass = sum(stat['pass_variants'] for stat in all_stats)
    total_gq_filtered = sum(stat['gq_filtered_variants'] for stat in all_stats)
    total_dp_filtered = sum(stat['dp_filtered_variants'] for stat in all_stats)
    total_rare = sum(stat['rare_variants'] for stat in all_stats)
    total_pathogenic = sum(stat['output_variants'] for stat in all_stats)
    total_acmg = sum(stat['acmg_yes_variants'] for stat in all_stats)

    print(f"統計情報の統合が完了しました: {total_samples}サンプル")
    print(f"総バリアント数: {total_variants}")
    print(f"PASS品質バリアント数: {total_pass} ({total_pass/total_variants*100:.2f}%)")
    print(f"GQフィルター通過数 (GQ>=20): {total_gq_filtered} ({total_gq_filtered/total_pass*100:.2f}%)")
    print(f"DPフィルター通過数 (DP>=10): {total_dp_filtered} ({total_dp_filtered/total_gq_filtered*100:.2f}%)")
    print(f"稀なバリアント数: {total_rare} ({total_rare/total_dp_filtered*100:.2f}%)")
    print(f"病原性バリアント数: {total_pathogenic}")
    print(f"ACMG遺伝子のバリアント数: {total_acmg}")
    print(f"統合統計ファイル: {merged_stats_file}")

def main():
    args = parse_arguments()

    print(f"=== 複数バッチの結果統合を開始します ===")
    print(f"注意: ClinVarのreview_status情報（星評価）を含めて統合します")

    # ディレクトリ確認
    if not os.path.exists(args.output_dir):
        print(f"エラー: 出力ディレクトリが存在しません: {args.output_dir}")
        return

    # 結果ファイルの統合
    merge_output_files(args.output_dir, args.num_batches)

    # 統計ファイルの統合
    merge_stats_files(args.output_dir, args.num_batches)

    print("\n=== 統合処理が完了しました ===")

if __name__ == "__main__":
    main()
