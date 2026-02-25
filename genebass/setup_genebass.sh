#!/usr/bin/env bash
# setup_genebass.sh — Install prerequisites and run GeneBass download
# Usage: bash setup_genebass.sh

set -euo pipefail

echo "============================================"
echo " GeneBass pLoF & Missense Download Setup"
echo "============================================"

# ---------- Detect OS ----------
OS="unknown"
if [[ "$OSTYPE" == "darwin"* ]]; then
    OS="mac"
elif [[ "$OSTYPE" == "linux-gnu"* ]]; then
    OS="linux"
fi
echo "Detected OS: $OS"

# ---------- Check / Install Java ----------
echo ""
echo "--- Checking Java ---"
if java -version 2>&1 | grep -q 'version "1[1-7]\.\|version "11\|version "17'; then
    echo "✓ Compatible Java found"
    java -version 2>&1 | head -1
else
    echo "Java 11+ not found. Installing..."
    if [[ "$OS" == "mac" ]]; then
        brew install openjdk@11
        export JAVA_HOME="$(brew --prefix openjdk@11)"
        export PATH="$JAVA_HOME/bin:$PATH"
        echo "Add to your shell profile:"
        echo "  export JAVA_HOME=$(brew --prefix openjdk@11)"
        echo "  export PATH=\$JAVA_HOME/bin:\$PATH"
    elif [[ "$OS" == "linux" ]]; then
        sudo apt-get update && sudo apt-get install -y openjdk-11-jdk
        export JAVA_HOME=/usr/lib/jvm/java-11-openjdk-amd64
    else
        echo "Please install Java 11 manually: https://adoptium.net/"
        exit 1
    fi
fi

# ---------- Check / Install Python ----------
echo ""
echo "--- Checking Python ---"
PYTHON_CMD=""
if command -v python3 &>/dev/null; then
    PYTHON_CMD="python3"
elif command -v python &>/dev/null; then
    PYTHON_CMD="python"
else
    echo "Python 3.8+ not found. Please install Python first."
    exit 1
fi
echo "✓ Using: $($PYTHON_CMD --version)"

# ---------- Install Hail ----------
echo ""
echo "--- Installing Hail ---"
$PYTHON_CMD -m pip install --upgrade pip
$PYTHON_CMD -m pip install hail
echo "✓ Hail installed"

# ---------- Run Download ----------
echo ""
echo "============================================"
echo " Running GeneBass Download"
echo "============================================"
echo ""
echo "This will download pLoF and missense gene burden results."
echo "Data source: gs://ukbb-exome-public/500k/results/results.mt"
echo ""
echo "NOTE: This dataset is large. Expect:"
echo "  - ~20-60 min download time"
echo "  - ~10-50 GB disk space for output CSVs"
echo "  - ~16 GB+ RAM recommended"
echo ""

$PYTHON_CMD download_genebass.py --annotation both --format csv --output-dir ./genebass_output

echo ""
echo "============================================"
echo " Done! Output in ./genebass_output/"
echo "============================================"
