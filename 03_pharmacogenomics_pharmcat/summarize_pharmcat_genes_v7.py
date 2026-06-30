#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# summarize_pharmcat_genes_v7.py
#
# 【スクリプト概要】
# ファイル名: summarize_pharmcat_genes_v7.py
# Pipeline 3 (PharmCAT), step 3: aggregate per-sample phenotypes.
# 処理内容:
#   - PharmCATの *.phenotype.json を解析し、DPYD/TPMT/UGT1A1/NUDT15/G6PD の phenotype を抽出
#   - phenotype.json の構造（geneReports -> CPIC -> <GENE> -> recommendationDiplotypes/sourceDiplotypes -> phenotypes）に対応
#   - 1行=1サンプル、各遺伝子=1列のTSV（phenotypeのみ）を output/ に出力
#   - 実行時間（開始・終了・経過秒）と、各遺伝子でNAでない件数をログに記録
#
# De-identification note:
#   入出力は --output-dir / ./output で指定します（ハードコードされた個人パスはありません）。
#   参加者個別データは本リポジトリに含まれません（README.md 参照）。

import argparse
import csv
import glob
import json
import os
import re
import sys
import time
from datetime import datetime
from typing import Any, Dict, List, Optional, Tuple

DEFAULT_GENES = ["DPYD", "TPMT", "UGT1A1", "NUDT15", "G6PD"]

def now_str() -> str:
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")

def norm(s: str) -> str:
    return re.sub(r"\s+", "", s.strip().upper())

def load_json(path: str) -> Optional[Any]:
    try:
        with open(path, "r", encoding="utf-8") as f:
            return json.load(f)
    except Exception as e:
        print(f"[WARN] JSON読み込み失敗: {path} ({e})", file=sys.stderr)
        return None

def parse_sample_prefix(prefix: str) -> Tuple[str, str, str]:
    """
    prefix例) SAMPLE001__<normal_uuid>
    戻り: (sample_prefix, tumor_sample, uuid)
    """
    tumor_sample = ""
    uuid = ""
    if "__" in prefix:
        tumor_sample, uuid = prefix.split("__", 1)
    return (prefix, tumor_sample, uuid)

def has_phenotype_json(dir_path: str) -> bool:
    return len(glob.glob(os.path.join(dir_path, "*.phenotype.json"))) > 0

def detect_output_dir(cli_output_dir: Optional[str]) -> str:
    """
    優先順位:
      1) --output-dir があればそれ
      2) スクリプト配置dirの ./output（*.phenotype.jsonがある）
      3) カレントdirの ./output（*.phenotype.jsonがある）
    """
    if cli_output_dir:
        od = os.path.abspath(cli_output_dir)
        if not os.path.isdir(od):
            raise FileNotFoundError(f"[ERROR] 指定された --output-dir が存在しません: {od}")
        return od

    script_dir = os.path.abspath(os.path.dirname(__file__))
    cand1 = os.path.join(script_dir, "output")
    if os.path.isdir(cand1) and has_phenotype_json(cand1):
        return os.path.abspath(cand1)

    cwd = os.path.abspath(".")
    cand2 = os.path.join(cwd, "output")
    if os.path.isdir(cand2) and has_phenotype_json(cand2):
        return os.path.abspath(cand2)

    msg = (
        "[ERROR] *.phenotype.json が見つかりません。\n"
        f"  - 探索した場所: {cand1}\n"
        f"  - 探索した場所: {cand2}\n"
        "対処:\n"
        "  - `--output-dir /path/to/output` を指定してください。\n"
    )
    raise FileNotFoundError(msg)

def safe_join_phenotypes(phs: Any) -> str:
    """
    phenotypes は通常 list[str]。複数ある場合は ';' で結合。
    """
    if phs is None:
        return "NA"
    if isinstance(phs, list):
        phs2 = [str(x).strip() for x in phs if str(x).strip() != ""]
        return ";".join(phs2) if phs2 else "NA"
    s = str(phs).strip()
    return s if s else "NA"

def extract_gene_phenotype_from_gene_report(gene_report: Dict[str, Any]) -> str:
    """
    gene report から phenotype を抽出する。
    優先順位:
      1) recommendationDiplotypes[0].phenotypes
      2) sourceDiplotypes[0].phenotypes
      3) （保険）gene_report.phenotypes
    """
    if not isinstance(gene_report, dict):
        return "NA"

    # 1) recommendationDiplotypes
    rec = gene_report.get("recommendationDiplotypes")
    if isinstance(rec, list) and len(rec) > 0 and isinstance(rec[0], dict):
        phs = rec[0].get("phenotypes")
        val = safe_join_phenotypes(phs)
        if val != "NA":
            return val

    # 2) sourceDiplotypes
    src = gene_report.get("sourceDiplotypes")
    if isinstance(src, list) and len(src) > 0 and isinstance(src[0], dict):
        phs = src[0].get("phenotypes")
        val = safe_join_phenotypes(phs)
        if val != "NA":
            return val

    # 3) gene_report 直下に phenotypes がある場合
    phs = gene_report.get("phenotypes")
    val = safe_join_phenotypes(phs)
    return val

def extract_gene_phenotype(pheno_obj: Any, gene: str) -> str:
    """
    phenotype.json 全体から gene の phenotype を抽出する。
    対応する構造:
      pheno_obj["geneReports"]["CPIC"][<GENE>] ...
    CPICが無い場合は、geneReports直下の他ソース（最初に見つかったもの）も試す。
    """
    if not isinstance(pheno_obj, dict):
        return "NA"

    gene_norm = norm(gene)

    gene_reports = pheno_obj.get("geneReports")
    if not isinstance(gene_reports, dict):
        return "NA"

    # ソースの優先順（CPIC優先、次に他）
    sources: List[str] = []
    if "CPIC" in gene_reports and isinstance(gene_reports["CPIC"], dict):
        sources.append("CPIC")
    # そのほか存在するソース
    for k, v in gene_reports.items():
        if k == "CPIC":
            continue
        if isinstance(v, dict):
            sources.append(k)

    for src in sources:
        src_block = gene_reports.get(src)
        if not isinstance(src_block, dict):
            continue

        # 1) まずキー一致（"UGT1A1" など）
        if gene in src_block and isinstance(src_block[gene], dict):
            val = extract_gene_phenotype_from_gene_report(src_block[gene])
            if val != "NA":
                return val

        # 2) 大文字小文字差・表記差に備えてキーを正規化して探索
        for k, v in src_block.items():
            if not isinstance(k, str):
                continue
            if norm(k) == gene_norm and isinstance(v, dict):
                val = extract_gene_phenotype_from_gene_report(v)
                if val != "NA":
                    return val

    return "NA"

def write_tsv(path: str, rows: List[Dict[str, str]], header: List[str]) -> None:
    with open(path, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=header, delimiter="\t")
        w.writeheader()
        for r in rows:
            w.writerow(r)

def main() -> None:
    start_ts = time.time()
    print(f"[INFO] Start: {now_str()}")

    ap = argparse.ArgumentParser(description="PharmCAT phenotype.jsonから5遺伝子（DPYD/TPMT/UGT1A1/NUDT15/G6PD）のphenotypeのみをwide形式で集計")
    ap.add_argument("--output-dir", default=None, help="PharmCAT出力ディレクトリ（*.phenotype.jsonがある場所）")
    ap.add_argument("--genes", nargs="+", default=DEFAULT_GENES, help="集計したい遺伝子名（複数指定可）")
    ap.add_argument("--out-file", default="oncology_5genes_phenotype_summary.tsv", help="出力TSVファイル名（output-dir内に作成）")
    ap.add_argument("--log-file", default="summarize_pharmcat_genes_v7.log", help="実行ログファイル名（output-dir内に作成）")
    args = ap.parse_args()

    outdir = detect_output_dir(args.output_dir)
    os.makedirs(outdir, exist_ok=True)

    log_path = os.path.join(outdir, args.log_file)

    def log(msg: str) -> None:
        print(msg)
        try:
            with open(log_path, "a", encoding="utf-8") as lf:
                lf.write(msg + "\n")
        except Exception:
            pass

    genes = [g.strip().upper() for g in args.genes if g.strip()]
    log(f"[INFO] Output dir: {outdir}")
    log(f"[INFO] Genes: {', '.join(genes)}")

    pheno_paths = sorted(glob.glob(os.path.join(outdir, "*.phenotype.json")))
    if not pheno_paths:
        log("[ERROR] *.phenotype.json が見つかりません。output-dirを確認してください。")
        sys.exit(1)

    # prefix一覧
    sample_prefixes: List[str] = []
    for pp in pheno_paths:
        bn = os.path.basename(pp)
        prefix = bn[:-len(".phenotype.json")]
        sample_prefixes.append(prefix)

    log(f"[INFO] Found samples (phenotype.json): {len(sample_prefixes)}")

    header = ["sample_prefix", "tumor_sample", "uuid"] + genes
    rows: List[Dict[str, str]] = []
    filled_counts = {g: 0 for g in genes}

    for i, prefix in enumerate(sample_prefixes, start=1):
        pheno_path = os.path.join(outdir, f"{prefix}.phenotype.json")
        pheno_obj = load_json(pheno_path)

        (sample_prefix, tumor_sample, uuid) = parse_sample_prefix(prefix)

        row: Dict[str, str] = {
            "sample_prefix": sample_prefix,
            "tumor_sample": tumor_sample,
            "uuid": uuid,
        }

        for gene in genes:
            val = extract_gene_phenotype(pheno_obj, gene)
            row[gene] = val
            if val != "NA":
                filled_counts[gene] += 1

        rows.append(row)

        if i % 25 == 0 or i == len(sample_prefixes):
            log(f"[INFO] Processed {i}/{len(sample_prefixes)} samples...")

    out_path = os.path.join(outdir, args.out_file)
    write_tsv(out_path, rows, header)
    log(f"[INFO] Wrote: {out_path} (rows={len(rows)})")

    log("[INFO] Phenotype filled counts (non-NA):")
    for g in genes:
        log(f"[INFO]   {g}: {filled_counts[g]}/{len(rows)}")

    end_ts = time.time()
    elapsed = end_ts - start_ts
    log(f"[INFO] End: {now_str()}")
    log(f"[INFO] Elapsed seconds: {elapsed:.2f}")

if __name__ == "__main__":
    main()
