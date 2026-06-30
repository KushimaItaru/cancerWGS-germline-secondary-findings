#!/bin/bash
# =============================================================================
# run_pharmcat_preprocess_v3.sh
# Pipeline 3 (PharmCAT), step 1: parallel PharmCAT VCF pre-processing.
#   For each sample's normal DeepVariant VCF, runs the PharmCAT VCF preprocessor
#   with --absent-to-ref to produce PharmCAT-ready .preprocessed.vcf.bgz files.
#   Reads sample pairs from a manifest (analysis_status == 1 rows only).
#   This is the exact as-run grid-engine script; an embedded Python driver is
#   written to a temp file and executed. Edit the placeholder paths.
#   No participant data is included (see README.md).
# =============================================================================

#$ -S /bin/bash
#$ -cwd
#$ -l s_vmem=24G,mem_req=24G
#$ -o pharmcat_preprocess.out
#$ -e pharmcat_preprocess.err

# Pythonの環境を設定
source ~/anaconda3/etc/profile.d/conda.sh
conda activate pharmcat

# 必要なモジュールをロード
module use /usr/local/package/modulefiles/
module load java/17.0.9.8.1

echo "バッチジョブからPharmCATプリプロセッサ（並列処理版）を実行します"
echo "開始時刻: $(date)"

# 一時Pythonスクリプトファイルを作成
TMP_PYTHON_SCRIPT=$(mktemp /tmp/pharmcat_preprocess_XXXXXX.py)

# Pythonスクリプトを一時ファイルに書き込む
cat > $TMP_PYTHON_SCRIPT << 'EOF'
#!/usr/bin/env python3
"""
並列処理対応PharmCATプリプロセッサ実行スクリプト
"""
import subprocess
import sys
import os
import time
import csv
import multiprocessing
from functools import partial

def ensure_htslib_tools():
    """htslibツール（bgzip、tabix）が利用可能か確認し、なければインストールする"""
    print("htslibツールを確認しています...")

    # 既存のbgzipを確認
    bgzip_path = os.path.expanduser("~/local/bin/bgzip")
    if os.path.exists(bgzip_path):
        print(f"bgzipが見つかりました: {bgzip_path}")
        # バージョン確認
        try:
            subprocess.run([bgzip_path, "--version"], check=True)
            print("bgzipは正常に動作しています")
            return True
        except:
            print("bgzipは存在しますが、動作に問題があります")
    else:
        print("bgzipが見つかりません")

    # htslibをインストール
    print("htslibをインストールしています...")
    build_dir = os.path.expanduser("~/htslib_build")
    install_dir = os.path.expanduser("~/local")

    # ディレクトリ作成
    os.makedirs(build_dir, exist_ok=True)

    # 現在のディレクトリを保存
    current_dir = os.getcwd()

    try:
        # ビルドディレクトリに移動
        os.chdir(build_dir)

        # GitHubから最新のhtslibソースをダウンロード
        print("htslibソースをダウンロードしています...")
        subprocess.run(
            ["wget", "https://github.com/samtools/htslib/releases/download/1.19/htslib-1.19.tar.bz2"],
            check=True
        )

        # 解凍
        print("ソースを解凍しています...")
        subprocess.run(["tar", "-xjf", "htslib-1.19.tar.bz2"], check=True)

        # ビルドディレクトリに移動
        os.chdir("htslib-1.19")

        # configureを実行
        print("configureを実行しています...")
        subprocess.run(
            ["./configure", f"--prefix={install_dir}"],
            check=True
        )

        # ビルド
        print("makeを実行しています...")
        subprocess.run(["make"], check=True)

        # インストール
        print("make installを実行しています...")
        subprocess.run(["make", "install"], check=True)

        # バージョン確認
        print("インストールされたbgzipのバージョンを確認しています...")
        version_result = subprocess.run(
            [f"{install_dir}/bin/bgzip", "--version"],
            capture_output=True, text=True, check=True
        )
        print(version_result.stdout)

        return True

    except Exception as e:
        print(f"htslibのビルド/インストール中にエラーが発生しました: {e}", file=sys.stderr)
        return False

    finally:
        # 元のディレクトリに戻る
        os.chdir(current_dir)

def read_sample_pairs(uuid_list_file):
    """UUIDリストファイルから腫瘍サンプル名と正常UUIDのペアを読み取る"""
    sample_pairs = []

    try:
        with open(uuid_list_file, 'r') as f:
            # ヘッダー行をスキップ
            for _ in range(2):
                next(f)

            # ヘッダー行を読み込む
            headers = next(f).strip().split()

            # 各行を処理
            for line in f:
                if line.strip() == "########### 共有未実施サンプル ###########":
                    break

                fields = line.strip().split()
                if len(fields) >= 4:  # 少なくとも必要な列があることを確認
                    tumor_sample = fields[0]  # 腫瘍サンプル名
                    normal_uuid = fields[3]   # 正常サンプルのUUID

                    # 1列目（analysis_status）が1であるものだけを処理
                    if len(fields) >= 5 and fields[4] == "1":
                        sample_pairs.append((tumor_sample, normal_uuid))

    except Exception as e:
        print(f"UUIDリストファイルの読み込み中にエラーが発生しました: {e}", file=sys.stderr)

    return sample_pairs

def process_sample(sample_pair, output_dir, pharmcat_script, index=0, total=0):
    """指定されたサンプルペアに対してPharmCATプリプロセッサを実行"""
    tumor_sample, normal_uuid = sample_pair

    # 入力VCFファイルのパス
    vcf_path = f"/path/to/wgs_results/{tumor_sample}/deepvariant/{normal_uuid}/{normal_uuid}.deepvariant.vcf.gz"

    # 出力ファイル名を新しい形式に変更
    output_prefix = f"{tumor_sample}__{normal_uuid}"

    # 入力ファイルの存在確認
    if not os.path.exists(vcf_path):
        print(f"エラー: 入力VCFファイルが見つかりません: {vcf_path}", file=sys.stderr)
        return (tumor_sample, normal_uuid, False)

    # コマンドと引数の設定
    cmd = [
        "python3",
        pharmcat_script,
        "-vcf",
        vcf_path,
        "--absent-to-ref",
        "-o",
        output_dir,
        "-bf",  # base filename
        output_prefix
    ]

    # スレッド数制限のための環境変数設定
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = "2"  # スレッド数をさらに制限（並列処理するため）
    env["OPENBLAS_NUM_THREADS"] = "2"
    env["MKL_NUM_THREADS"] = "2"
    env["NUMEXPR_NUM_THREADS"] = "2"

    # bgzipとbcftoolsのパスを環境変数に設定
    local_bin = os.path.expanduser("~/local/bin")
    env["PATH"] = f"{local_bin}:" + env.get("PATH", "")

    # 実行
    print(f"=== サンプル {index+1}/{total} の処理を開始 ===")
    print(f"サンプル: {tumor_sample} / 正常UUID: {normal_uuid}")
    print(f"プロセスID: {os.getpid()}")
    start_time = time.time()

    try:
        # 標準出力・エラー出力をファイルにリダイレクト（並列処理時の出力を整理するため）
        log_file = os.path.join(output_dir, f"{output_prefix}.preprocess.log")
        with open(log_file, 'w') as log:
            result = subprocess.run(
                cmd,
                env=env,
                text=True,
                stdout=log,
                stderr=subprocess.STDOUT
            )

        end_time = time.time()
        elapsed_time = end_time - start_time

        success = result.returncode == 0
        status = "成功" if success else "失敗"

        print(f"処理が完了しました: {status} (所要時間: {elapsed_time:.2f}秒)")
        return (tumor_sample, normal_uuid, success)

    except Exception as e:
        print(f"実行中に例外が発生しました: {e}", file=sys.stderr)
        return (tumor_sample, normal_uuid, False)

def main():
    print("HPCクラスタでPharmCATプリプロセッサを実行します（並列処理版）")

    # htslibツールの確認とインストール
    if not ensure_htslib_tools():
        print("htslibツールのインストールに失敗しました。")
        sys.exit(1)

    # javaモジュールのロード
    print("\n必要なモジュールをロードしています...")
    try:
        subprocess.run("module use /usr/local/package/modulefiles/", shell=True, check=True)
        subprocess.run("module load java/17.0.9.8.1", shell=True, check=True)
        print("javaモジュールのロードが完了しました")
    except subprocess.CalledProcessError as e:
        print(f"モジュールのロード中にエラーが発生しました: {e}", file=sys.stderr)
        sys.exit(1)

    # 出力ディレクトリの確認と作成
    output_dir = "/path/to/pharmcat/output"
    if not os.path.exists(output_dir):
        try:
            os.makedirs(output_dir)
            print(f"出力ディレクトリを作成しました: {output_dir}")
        except OSError as e:
            print(f"出力ディレクトリの作成に失敗しました: {e}", file=sys.stderr)
            sys.exit(1)
    else:
        # ログディレクトリを作成（既存のディレクトリは保持）
        log_dir = os.path.join(output_dir, "logs")
        os.makedirs(log_dir, exist_ok=True)

    # pharmcat_vcf_preprocessor.pyのパス
    pharmcat_script = "/path/to/pharmcat/pharmcat_vcf_preprocessor.py"

    # UUIDリストファイルからサンプルペアを読み込む
    uuid_list_file = "/path/to/sample_manifest.txt"
    sample_pairs = read_sample_pairs(uuid_list_file)

    total_samples = len(sample_pairs)
    print(f"処理対象サンプル数: {total_samples}")

    # 利用可能なCPUコア数を取得（ジョブに割り当てられたコア数に応じて調整）
    # SGEの環境変数から取得するか、またはデフォルト値を使用
    num_cores = int(os.environ.get('NSLOTS', multiprocessing.cpu_count()))
    # メモリ使用量を考慮して、実際に使用するプロセス数を決定
    # メモリ使用量が高いため、コア数の半分程度に制限
    num_processes = min(10, max(2, num_cores // 2))

    print(f"並列処理するプロセス数: {num_processes} (利用可能コア数: {num_cores})")

    # 並列処理の準備
    process_func = partial(
        process_sample,
        output_dir=output_dir,
        pharmcat_script=pharmcat_script,
        total=total_samples
    )

    # プロセスプールを作成して並列処理を実行
    results = []
    with multiprocessing.Pool(processes=num_processes) as pool:
        # サンプルとインデックスのペアを作成
        indexed_samples = [(sample, i) for i, sample in enumerate(sample_pairs)]

        # 部分関数を使わない場合の書き方
        results = pool.starmap(
            process_sample,
            [(sample, output_dir, pharmcat_script, i, total_samples) for i, sample in enumerate(sample_pairs)]
        )

    # 処理結果を集計
    successful = [r for r in results if r[2]]
    failed = [r for r in results if not r[2]]

    # 処理結果の表示
    print("\n=== 処理結果サマリー ===")
    print(f"合計サンプル数: {total_samples}")
    print(f"成功: {len(successful)}")
    print(f"失敗: {len(failed)}")

    if failed:
        print("\n失敗したサンプル:")
        for tumor_sample, normal_uuid, _ in failed:
            print(f"  - {tumor_sample} ({normal_uuid})")

    # 出力ファイルの確認
    if os.path.exists(output_dir):
        output_files = [f for f in os.listdir(output_dir) if os.path.isfile(os.path.join(output_dir, f)) and f.endswith('.bgz')]
        if output_files:
            print(f"\n主な出力ファイル数: {len(output_files)}")
            print(f"例: {', '.join(output_files[:3])}...")
        else:
            print("\n警告: 出力ディレクトリに.bgzファイルが見つかりません")
    else:
        print(f"\n警告: 出力ディレクトリ {output_dir} が見つかりません")

if __name__ == "__main__":
    main()
EOF

# 一時スクリプトに実行権限を付与
chmod +x $TMP_PYTHON_SCRIPT

# スクリプトを実行
python3 $TMP_PYTHON_SCRIPT

# 一時ファイルを削除
rm -f $TMP_PYTHON_SCRIPT

echo "終了時刻: $(date)"
