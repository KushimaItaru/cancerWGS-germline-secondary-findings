#!/bin/bash
# =============================================================================
# install_annotsv_local.sh
# Pipeline 2 (CNV), helper: local (no-root) installation of AnnotSV v3.3.6 into
#   ${HOME}/local/annotsv and creation of a setup_annotsv.sh environment file.
#   Provided for reproducibility of the annotation toolchain. No participant data.
# =============================================================================

echo "=== Installing AnnotSV locally ==="

# 必要なモジュールをロード
module load gcc
module load python/3.9

# インストールディレクトリ
INSTALL_DIR="${HOME}/local/annotsv"
mkdir -p ${INSTALL_DIR}
cd ${INSTALL_DIR}

# AnnotSVのダウンロード
echo "1. Downloading AnnotSV v3.3.6..."
wget https://github.com/lgmgeo/AnnotSV/archive/refs/tags/v3.3.6.tar.gz

if [ ! -f "v3.3.6.tar.gz" ]; then
    echo "ERROR: Failed to download AnnotSV"
    exit 1
fi

# 解凍
echo "2. Extracting..."
tar -xzf v3.3.6.tar.gz
cd AnnotSV-3.3.6

# AnnotSVは実際にはインストール不要（スクリプトベース）
echo "3. Setting up AnnotSV..."

# 必要なPerlモジュールの確認
echo "4. Checking Perl modules..."
perl -e "use YAML::XS;" 2>/dev/null || echo "  Warning: YAML::XS not found"
perl -e "use Sort::Key::Natural;" 2>/dev/null || echo "  Warning: Sort::Key::Natural not found"

# 環境設定ファイルの作成
ANNOTSV_PATH="${INSTALL_DIR}/AnnotSV-3.3.6"
cat > ${HOME}/setup_annotsv.sh << EOF
#!/bin/bash
# AnnotSV environment setup
export ANNOTSV="${ANNOTSV_PATH}"
export PATH="\${ANNOTSV}/bin:\${PATH}"
export PERL5LIB="\${ANNOTSV}/share/perl5:\${PERL5LIB}"

# TCL設定（AnnotSVはTclで書かれている）
module load tcl 2>/dev/null || true

echo "AnnotSV environment loaded"
echo "ANNOTSV=\${ANNOTSV}"
EOF

chmod +x ${HOME}/setup_annotsv.sh

# アノテーションデータのダウンロード設定
echo ""
echo "5. Setting up annotation databases..."
mkdir -p ${ANNOTSV_PATH}/share/AnnotSV/Annotations_Human

# 基本的なテスト
echo ""
echo "6. Testing installation..."
source ${HOME}/setup_annotsv.sh
${ANNOTSV}/bin/AnnotSV -help 2>&1 | head -20 || echo "Note: AnnotSV may need additional setup"

echo ""
echo "Installation completed!"
echo ""
echo "To use AnnotSV:"
echo "1. source ${HOME}/setup_annotsv.sh"
echo "2. AnnotSV -SVinputFile <your.bed> -outputFile <output>"
