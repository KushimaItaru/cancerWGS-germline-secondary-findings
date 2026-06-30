#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# =============================================================================
# variant_filter_parallel_v4.py
# Pipeline 1 (SNV/indel): extract Pathogenic / Likely-pathogenic germline
#   SNV/indel records from per-sample ANNOVAR-annotated tumour-normal WGS output.
#
# What it does:
#   - Loads a sample manifest, an ACMG/return-gene list (SF_106Genes.txt), and a
#     ClinVar variant_summary table (AlleleID -> VariationID / CLNSIG / review status).
#   - For each sample's *.hg38_multianno.txt (DeepVariant + ANNOVAR / InterVar):
#       * keep FILTER == PASS, GQ >= 20, DP >= 10 (FORMAT field parsed by name)
#       * call P/LP by ClinVar (excluding "Conflicting") OR InterVar_automated
#       * annotate ClinVar VariationID / CLNSIG / review-status stars and SF-gene membership
#   - Writes per-sample filtered tables + per-batch statistics (array-job design).
#   - All columns are resolved dynamically by header NAME (no fixed indices).
#
# De-identification note:
#   File-system paths below are placeholders ("/path/to/..."). Edit them (or pass
#   --uuid-list / --clinvar-file / --output-dir / --acmg-gene-file) before running.
#   No participant-level data is included in this repository. See README.md.
# =============================================================================

import os
import re
import gzip
import sys
import argparse
import time

def parse_arguments():
    """コマンドライン引数を解析"""
    parser = argparse.ArgumentParser(description='バリアントフィルタリングを複数サンプルに対して実行')
    parser.add_argument('--batch-index', type=int, required=True, help='処理するバッチのインデックス（0から始まる）')
    parser.add_argument('--batch-size', type=int, default=12, help='バッチあたりのサンプル数（デフォルト：12）')
    parser.add_argument('--uuid-list', type=str, default='/path/to/sample_manifest.txt',
                        help='UUID情報を含むファイルのパス')
    parser.add_argument('--output-dir', type=str, default='./filtered_variants',
                        help='出力ディレクトリ')
    parser.add_argument('--acmg-gene-file', type=str, default='SF_106Genes.txt',
                        help='ACMG遺伝子リストファイル')
    parser.add_argument('--clinvar-file', type=str, default='/path/to/clinvar/variant_summary_GRCh38.txt.gz',
                        help='ClinVarデータファイル (variant_summary_GRCh38.txt.gz)')
    return parser.parse_args()

def load_uuid_list(uuid_list_file):
    """uuid_listファイルからサンプル情報を読み込む"""
    sample_info = []
    try:
        with open(uuid_list_file, 'r') as f:
            header_found = False
            for line in f:
                if line.startswith("tumor_sample_name"):
                    header_found = True
                    continue

                if header_found and line.strip() and not line.startswith("#"):
                    fields = line.strip().split()
                    if len(fields) >= 4:  # tumor_sample_name, tumor_uuid, normal_sample_name, normal_uuid
                        sample_info.append({
                            "tumor_sample_name": fields[0],
                            "tumor_uuid": fields[1],
                            "normal_sample_name": fields[2],
                            "normal_uuid": fields[3]
                        })
        print(f"uuid_listから{len(sample_info)}サンプルの情報を読み込みました")
    except Exception as e:
        print(f"エラー: uuid_listファイルの読み込み中にエラーが発生しました: {e}")
        sys.exit(1)

    return sample_info

def load_acmg_genes(gene_file):
    """SF 遺伝子リストを読み込む"""
    acmg_genes = set()
    try:
        with open(gene_file, 'r') as f:
            for line in f:
                gene = line.strip()
                if gene:  # 空行をスキップ
                    acmg_genes.add(gene)
        print(f"SF 遺伝子リスト: {len(acmg_genes)}遺伝子を読み込みました")
    except Exception as e:
        print(f"エラー: 遺伝子リストファイルの読み込み中にエラーが発生しました: {e}")
        sys.exit(1)

    return acmg_genes

def get_review_stars(review_status):
    """ReviewStatusを星評価に変換する関数"""
    if review_status == "practice guideline":
        return "★★★★"
    elif review_status == "reviewed by expert panel":
        return "★★★"
    elif review_status == "criteria provided, multiple submitters, no conflicts":
        return "★★"
    elif review_status == "criteria provided, single submitter":
        return "★"
    elif "conflicting" in review_status:
        return "★0"
    else:
        return "★0"  # 評価無し、または不明な場合

def load_clinvar_data(clinvar_file):
    """variant_summary_GRCh38.txtからAlleleID、VariationID、ClinicalSignificance、ReviewStatus情報を読み込む"""
    clinvar_data = {}

    # gzipファイルを開く
    print("ClinVarデータを読み込んでいます...")
    line_count = 0
    entry_count = 0

    try:
        with gzip.open(clinvar_file, 'rt') as f:
            # ヘッダー行を読み取り、必要な列のインデックスを特定
            header = f.readline().strip()
            if not header.startswith('#'):
                header = '#' + header  # #が省略されている場合に追加

            header_fields = header.strip('#').split('\t')

            try:
                alleleid_idx = header_fields.index("AlleleID")
                variationid_idx = header_fields.index("VariationID")
                clinsig_idx = header_fields.index("ClinicalSignificance")
                review_idx = header_fields.index("ReviewStatus")
            except ValueError as e:
                print(f"エラー: 必要な列が見つかりません: {e}")
                print(f"ヘッダー行: {header}")
                sys.exit(1)

            # データ行を処理
            for line in f:
                line_count += 1
                if line_count % 500000 == 0:
                    print(f"  {line_count}行を処理中... {entry_count}エントリ読み込み完了")

                fields = line.strip().split('\t')
                if len(fields) <= max(alleleid_idx, variationid_idx, clinsig_idx, review_idx):
                    continue

                try:
                    alleleid = fields[alleleid_idx]
                    variationid = fields[variationid_idx]
                    clinsig = fields[clinsig_idx]
                    review_status = fields[review_idx]
                    review_stars = get_review_stars(review_status)

                    # 空でないエントリのみ追加
                    if alleleid and variationid and clinsig:
                        clinvar_data[alleleid] = (variationid, clinsig, review_stars)
                        entry_count += 1
                except Exception as e:
                    # 特定の行の処理中にエラーが発生した場合はスキップして続行
                    print(f"警告: 行 {line_count} の処理中にエラーが発生しました: {e}")
                    continue

    except Exception as e:
        print(f"エラー: ClinVarデータの読み込み中にエラーが発生しました: {e}")
        sys.exit(1)

    print(f"ClinVarデータの読み込みが完了しました: {entry_count}エントリを読み込みました")
    return clinvar_data

def extract_pathogenic_variants(input_txt, output_txt, acmg_genes, clinvar_data, sample_id):
    """バリアントを抽出し加工する関数"""

    # カウンタの初期化
    total_variants = 0
    pass_variants = 0
    gq_filtered_variants = 0  # GQフィルターを通過したバリアント数
    dp_filtered_variants = 0  # DPフィルターを通過したバリアント数
    rare_variants = 0
    clinvar_matches = 0
    clinvar_potential_pathogenic = 0  # 矛盾する分類を含む病原性バリアント
    clinvar_pathogenic = 0  # 矛盾する分類を除外した病原性バリアント
    conflicting_variants = 0
    intervar_pathogenic = 0
    output_variants = 0
    acmg_yes_variants = 0

    # 出力ディレクトリを作成（存在しない場合）
    output_dir = os.path.dirname(output_txt)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)
        print(f"出力ディレクトリを作成しました: {output_dir}")

    # 病原性判定用の正規表現パターン
    pathogenic_pattern = re.compile(r'pathogenic|Pathogenic|likely_pathogenic|Likely_Pathogenic', re.IGNORECASE)
    conflicting_pattern = re.compile(r'Conflicting_classifications_of_pathogenicity|conflicting|Conflicting', re.IGNORECASE)

    print(f"\n1. サンプル {sample_id} の入力ファイルの読み込みを開始します...")

    try:
        with open(input_txt, 'r') as infile, open(output_txt, 'w') as outfile:
            # ヘッダー行を読み取り、インデックスを特定
            header = infile.readline().strip()
            header_fields = header.split('\t')

            try:
                intervar_idx = header_fields.index("InterVar_automated")
            except ValueError:
                intervar_idx = 82  # デフォルト値

            try:
                gene_idx = header_fields.index("Gene.refGene")
            except ValueError:
                print("Gene.refGene列が見つかりません。スクリプトを終了します。")
                return

            try:
                otherinfo10_idx = header_fields.index("Otherinfo10")
            except ValueError:
                print("Otherinfo10列が見つかりません。全ての列名を表示します：")
                print(header_fields)
                return

            try:
                otherinfo13_idx = header_fields.index("Otherinfo13")
            except ValueError:
                print("Otherinfo13列が見つかりません。全ての列名を表示します：")
                print(header_fields)
                return

            try:
                chr_idx = header_fields.index("Chr")
            except ValueError:
                chr_idx = 0  # デフォルト値

            try:
                start_idx = header_fields.index("Start")
            except ValueError:
                start_idx = 1  # デフォルト値

            try:
                ref_idx = header_fields.index("Ref")
            except ValueError:
                ref_idx = 3  # デフォルト値

            try:
                alt_idx = header_fields.index("Alt")
            except ValueError:
                alt_idx = 4  # デフォルト値

            try:
                clinvar_idx = header_fields.index("clinvar_20211010")
            except ValueError:
                print("clinvar_20211010列が見つかりません。全ての列名を表示します：")
                print(header_fields)
                return

            try:
                tommo_idx = header_fields.index("ToMMo")
            except ValueError:
                print("ToMMo列が見つかりません。全ての列名を表示します：")
                print(header_fields)
                return

            # 新しいヘッダー行を作成（SampleID列を追加）
            new_header_fields = ["SampleID"] + header_fields.copy()
            new_header_fields.insert(gene_idx + 2, "SF_106Genes")  # SampleIDを追加したので+2
            new_header_fields.insert(clinvar_idx + 3, "clinvarID")  # SampleIDとSF_106Genesを追加したので+3
            new_header_fields.insert(clinvar_idx + 4, "clinvar_CLNSIG")  # SampleID、SF_106Genes、clinvarIDを追加したので+4
            new_header_fields.insert(clinvar_idx + 5, "review_status")  # 新しい列を追加
            new_header = '\t'.join(new_header_fields)
            outfile.write(new_header + '\n')

            # 行カウント用
            line_count = 0
            progress_step = 1000000  # 100万行ごとに進捗を表示

            # 各行を処理
            for line in infile:
                line_count += 1
                total_variants += 1  # 総バリアント数をカウント

                # 進捗表示
                if line_count % progress_step == 0:
                    print(f"  {line_count}行を処理中... ({total_variants}バリアント)")

                fields = line.strip().split('\t')
                if len(fields) <= max(intervar_idx, gene_idx, otherinfo10_idx, otherinfo13_idx, chr_idx, start_idx, ref_idx, alt_idx, clinvar_idx, tommo_idx):
                    continue

                intervar_value = fields[intervar_idx]
                otherinfo10_value = fields[otherinfo10_idx]
                otherinfo13_value = fields[otherinfo13_idx]
                tommo_value = fields[tommo_idx]

                # フィルタリングステップ2: 品質フィルター - 完全一致に変更
                if otherinfo10_value != "PASS":
                    continue

                pass_variants += 1

                # 新規フィルタリングステップ: GQ>=20 と DP>=10
                # Otherinfo13の形式は "GT:GQ:DP:AD:VAF:PL" の想定
                otherinfo13_parts = otherinfo13_value.split(':')

                # フォーマットチェックと欠損値チェック
                if len(otherinfo13_parts) < 3 or otherinfo13_parts[1] == '.' or otherinfo13_parts[2] == '.':
                    continue  # GQまたはDPが欠損している場合はスキップ

                try:
                    gq_value = int(otherinfo13_parts[1])
                    dp_value = int(otherinfo13_parts[2])
                except ValueError:
                    continue  # 数値変換に失敗した場合はスキップ

                # GQフィルター
                if gq_value < 20:
                    continue

                gq_filtered_variants += 1

                # DPフィルター
                if dp_value < 10:
                    continue

                dp_filtered_variants += 1

                # フィルタリングステップ3: ToMMo頻度フィルターは削除
                # 代わりに全てのバリアントを通過させる
                is_rare = True  # 常にTrueを設定
                rare_variants += 1

                # ステップ4: アノテーション付与（ClinVarデータとマッチング）
                clinvar_id = fields[clinvar_idx]  # clinvar_20211010列からALLELEIDを取得

                # ClinVarデータを取得（VariationID、ClinicalSignificance、ReviewStatus）
                clnsig = "Not_in_ClinVar"
                clinvar_id_value = "."
                review_stars = "★0"  # デフォルト値
                if clinvar_id and clinvar_id != ".":
                    clinvar_info = clinvar_data.get(clinvar_id, (clinvar_id, "Not_in_ClinVar", "★0"))
                    if len(clinvar_info) == 3:  # 更新版のデータ形式
                        clinvar_id_value, clnsig, review_stars = clinvar_info
                    else:  # 旧バージョンのデータ形式（後方互換性）
                        clinvar_id_value, clnsig = clinvar_info

                if clnsig != "Not_in_ClinVar":
                    clinvar_matches += 1

                # フィルタリングステップ5: 病原性フィルター

                # ClinVarでpathogenicとされているかチェック
                is_clinvar_pathogenic = pathogenic_pattern.search(clnsig) is not None

                # Conflicting_classifications_of_pathogenicityをチェック
                is_conflicting = conflicting_pattern.search(clnsig) is not None

                # 潜在的な病原性バリアント(除外前)をカウント
                if is_clinvar_pathogenic:
                    clinvar_potential_pathogenic += 1

                    # 矛盾する分類のチェック
                    if is_conflicting:
                        conflicting_variants += 1
                        continue  # 矛盾する分類を持つバリアントはスキップ
                    else:
                        # 矛盾がなく、実際に病原性のあるバリアントをカウント
                        clinvar_pathogenic += 1

                # InterVarで病原性とされているかチェック
                is_intervar_pathogenic = pathogenic_pattern.search(intervar_value) is not None
                if is_intervar_pathogenic:
                    intervar_pathogenic += 1

                # 病原性条件: InterVarかClinVarのどちらかが病原性を示す（フレームシフト条件を削除）
                is_pathogenic = is_intervar_pathogenic or (is_clinvar_pathogenic and not is_conflicting)

                if not is_pathogenic:
                    continue

                output_variants += 1  # 出力バリアント数をカウント

                # ステップ6: ACMG情報の付与
                gene_value = fields[gene_idx]

                # 遺伝子がACMG SF v3.2リストに含まれるか確認
                genes = gene_value.split(';')
                is_acmg = "No"
                for gene in genes:
                    # 一部の特殊なフォーマットを処理
                    if '=' in gene:
                        continue

                    # 遺伝子名に余分な情報があれば除去
                    gene = gene.split('(')[0].strip()

                    if gene in acmg_genes:
                        is_acmg = "Yes"
                        acmg_yes_variants += 1  # ACMG Yes バリアントをカウント
                        break

                # 新しい行を作成（SampleID列と新しい列を追加）
                new_fields = [sample_id] + fields.copy()
                new_fields.insert(gene_idx + 2, is_acmg)  # SampleIDを追加したので+2
                new_fields.insert(clinvar_idx + 3, clinvar_id_value)  # SampleIDとSF_106Genesを追加したので+3
                new_fields.insert(clinvar_idx + 4, clnsig)  # SampleID、SF_106Genes、clinvarIDを追加したので+4
                new_fields.insert(clinvar_idx + 5, review_stars)  # レビュー星評価を追加
                new_line = '\t'.join(new_fields)
                outfile.write(new_line + '\n')
    except Exception as e:
        print(f"エラー: バリアント抽出中にエラーが発生しました: {e}")
        return None

    # 各ステップの結果を表示
    print(f"\n----- サンプル {sample_id} のフィルタリング結果 -----")
    print(f"1. 入力ファイルの読み込み: {total_variants}バリアント")
    print(f"2. 品質フィルター (Otherinfo10がPASSに完全一致): {pass_variants}バリアント ({pass_variants/total_variants*100:.2f}%)")
    print(f"3. GQ品質フィルター (GQ >= 20): {gq_filtered_variants}バリアント ({gq_filtered_variants/pass_variants*100:.2f}% of PASS)")
    print(f"4. DP品質フィルター (DP >= 10): {dp_filtered_variants}バリアント ({dp_filtered_variants/gq_filtered_variants*100:.2f}% of GQ filtered)")
    print(f"5. 頻度フィルター: ToMMoフィルターは削除されました (全てのバリアントを通過)")
    print(f"6. 病原性判定:")
    print(f"   - ClinVar: {clinvar_pathogenic}件の病原性バリアント（他に矛盾する分類{conflicting_variants}件は除外）")
    print(f"   - InterVar: {intervar_pathogenic}件の病原性バリアント")
    print(f"7. 最終的な病原性バリアント: {output_variants}件")
    print(f"8. SF_106Genes遺伝子に含まれるバリアント: {acmg_yes_variants}件")

    # 統計情報を保存
    stats_file = output_txt + ".stats"
    with open(stats_file, 'w') as f:
        f.write(f"sample_id\t{sample_id}\n")
        f.write(f"total_variants\t{total_variants}\n")
        f.write(f"pass_variants\t{pass_variants}\n")
        f.write(f"gq_filtered_variants\t{gq_filtered_variants}\n")
        f.write(f"dp_filtered_variants\t{dp_filtered_variants}\n")
        f.write(f"rare_variants\t{rare_variants}\n")
        f.write(f"output_variants\t{output_variants}\n")
        f.write(f"acmg_yes_variants\t{acmg_yes_variants}\n")
        f.write(f"clinvar_potential_pathogenic\t{clinvar_potential_pathogenic}\n")
        f.write(f"clinvar_pathogenic\t{clinvar_pathogenic}\n")
        f.write(f"intervar_pathogenic\t{intervar_pathogenic}\n")
        f.write(f"conflicting_variants\t{conflicting_variants}\n")

    # 統計情報を返す
    return {
        "sample_id": sample_id,
        "total_variants": total_variants,
        "pass_variants": pass_variants,
        "gq_filtered_variants": gq_filtered_variants,
        "dp_filtered_variants": dp_filtered_variants,
        "rare_variants": rare_variants,
        "output_variants": output_variants,
        "acmg_yes_variants": acmg_yes_variants,
        "clinvar_potential_pathogenic": clinvar_potential_pathogenic,
        "clinvar_pathogenic": clinvar_pathogenic,
        "intervar_pathogenic": intervar_pathogenic,
        "conflicting_variants": conflicting_variants
    }

def main():
    start_time = time.time()
    args = parse_arguments()

    print(f"=== バッチ {args.batch_index} のバリアントフィルタリングプロセスを開始します ===\n")
    print(f"実行日時: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"SF遺伝子ファイル: {args.acmg_gene_file}")
    print(f"ClinVarファイル: {args.clinvar_file}")
    print(f"UUIDリスト: {args.uuid_list}")
    print(f"注意: ToMMo頻度フィルターは削除されました - 全てのバリアントが頻度フィルターを通過します")
    print(f"注意: ClinVarにReviewStatus情報（星評価）が追加されました")

    # 出力ディレクトリが存在しない場合は作成
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)
        print(f"出力ディレクトリを作成しました: {args.output_dir}")

    # サンプル情報を読み込む
    all_samples = load_uuid_list(args.uuid_list)

    # 処理するバッチを決定
    start_idx = args.batch_index * args.batch_size
    end_idx = min(start_idx + args.batch_size, len(all_samples))

    # バッチ内のサンプルを抽出
    batch_samples = all_samples[start_idx:end_idx]

    print(f"バッチ {args.batch_index}: {len(batch_samples)}サンプルを処理します（全{len(all_samples)}サンプル中）")
    print(f"サンプルインデックス範囲: {start_idx}-{end_idx - 1}")

    if not batch_samples:
        print("処理するサンプルがありません。終了します。")
        return

    # ACMG遺伝子リストを読み込む
    acmg_genes = load_acmg_genes(args.acmg_gene_file)

    # ClinVarデータを読み込む
    clinvar_data = load_clinvar_data(args.clinvar_file)

    # ステップ1-6: 各サンプルのフィルタリングを実行
    print("\n=== サンプルごとのフィルタリングステップ ===")

    all_stats = []
    batch_output_file = os.path.join(args.output_dir, f"batch_{args.batch_index}_results.txt")
    batch_stats_file = os.path.join(args.output_dir, f"batch_{args.batch_index}_stats.txt")

    with open(batch_stats_file, 'w') as stats_file:
        stats_file.write("SampleID\t総バリアント数\tPASSバリアント\tGQフィルター\tDPフィルター\t稀なバリアント\t病原性バリアント\tACMG遺伝子\n")

    for idx, sample in enumerate(batch_samples):
        sample_start_time = time.time()
        tumor_sample_name = sample["tumor_sample_name"]  # AAA
        normal_uuid = sample["normal_uuid"]  # BBB

        # 入力ファイルパスを構築
        input_txt = f"/path/to/wgs_results/{tumor_sample_name}/deepvariant/{normal_uuid}/{normal_uuid}.deepvariant.vcf.gz.hg38_multianno.txt"

        # 出力ファイルパスを構築
        sample_output_dir = os.path.join(args.output_dir, tumor_sample_name)
        output_txt = os.path.join(sample_output_dir, f"{normal_uuid}.deepvariant.vcf.gz.hg38_multianno_filtered.txt")

        print(f"\n処理中: サンプル {idx+1}/{len(batch_samples)} - {tumor_sample_name}")
        print(f"入力ファイル: {input_txt}")
        print(f"出力ファイル: {output_txt}")

        # ファイルが存在するか確認
        if not os.path.exists(input_txt):
            print(f"警告: 入力ファイルが見つかりません: {input_txt}")
            print(f"このサンプルはスキップします")
            continue

        # バリアントをフィルタリング
        stats = extract_pathogenic_variants(input_txt, output_txt, acmg_genes, clinvar_data, tumor_sample_name)
        if stats:
            all_stats.append(stats)

            # 統計情報をファイルに追記
            with open(batch_stats_file, 'a') as stats_file:
                stats_file.write(f"{stats['sample_id']}\t{stats['total_variants']}\t{stats['pass_variants']}\t{stats['gq_filtered_variants']}\t{stats['dp_filtered_variants']}\t{stats['rare_variants']}\t{stats['output_variants']}\t{stats['acmg_yes_variants']}\n")

        sample_end_time = time.time()
        print(f"サンプル処理時間: {sample_end_time - sample_start_time:.2f}秒")

    # バッチ情報ファイルを作成（マージスクリプト用）
    batch_info_file = os.path.join(args.output_dir, f"batch_{args.batch_index}_info.txt")
    with open(batch_info_file, 'w') as f:
        for sample in batch_samples:
            tumor_sample_name = sample["tumor_sample_name"]
            normal_uuid = sample["normal_uuid"]
            sample_output_dir = os.path.join(args.output_dir, tumor_sample_name)
            output_txt = os.path.join(sample_output_dir, f"{normal_uuid}.deepvariant.vcf.gz.hg38_multianno_filtered.txt")
            f.write(f"{output_txt}\n")

    end_time = time.time()
    elapsed_time = end_time - start_time

    print(f"\n=== バッチ {args.batch_index} の処理が完了しました ===")
    print(f"バッチ情報: {batch_info_file}")
    print(f"バッチ統計: {batch_stats_file}")
    print(f"合計処理時間: {elapsed_time:.2f}秒 ({elapsed_time/60:.2f}分)")
    print(f"処理完了時刻: {time.strftime('%Y-%m-%d %H:%M:%S')}")

    print("\n処理が正常に完了しました")

if __name__ == "__main__":
    main()
