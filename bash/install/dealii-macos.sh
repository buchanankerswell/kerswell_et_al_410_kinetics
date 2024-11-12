#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

# ===========================
# Get args
# ===========================
if [ $# -ne 1 ]; then
  echo "Usage: $0 <DEALII_VERSION>"
  exit 1
fi

DEALII_VERSION=$1

# ===========================
# Cleanup handler
# ===========================
cleanup() {
    local exit_code=$?
    if [ -n "${TMP_CANDI_DIR:-}" ] && [ -d "$TMP_CANDI_DIR" ]; then
        echo "  Cleaning up temporary candi directory:"
        echo "  $TMP_CANDI_DIR"
        rm -rf "$TMP_CANDI_DIR"
    fi

    if [ $exit_code -ne 0 ]; then
        echo "  Installation failed with exit code: $exit_code"
    fi

    echo "  Cleanup completed."
    echo "============================================="
    exit $exit_code
}

trap cleanup EXIT INT TERM

# ===========================
# Environment setup
# ===========================
echo "============================================="
echo "  Setting environment variables ..."

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
NPROC=${NPROC:-$(sysctl -n hw.ncpu)}
DEALII_DIR=$PROJECT_ROOT/dealii

echo "  PROJCT_ROOT:      $PROJECT_ROOT"
echo "  NPROC:            $NPROC"
echo "  DEALII_DIR:       $DEALII_DIR"

# ===========================
# Skip if already built
# ===========================
if [ -d "$DEALII_DIR/deal.II-$DEALII_VERSION" ]; then
    echo "  ==========================================="
    echo "  deal.II-$DEALII_VERSION found at:"
    echo "    $DEALII_DIR/deal.II-$DEALII_VERSION"
    exit 0
fi

# ===========================
# Clone candi
# ===========================
echo "  ==========================================="
TMP_CANDI_DIR=$(mktemp -d)
echo "  Cloning candi into tmp dir: $TMP_CANDI_DIR"
git clone https://github.com/buchanankerswell/candi.git "$TMP_CANDI_DIR" > /dev/null 2>&1
cd "$TMP_CANDI_DIR"
git checkout update-macos-platform
curl -O https://raw.githubusercontent.com/geodynamics/aspect/main/contrib/install/local.cfg
cat << EOF >> local.cfg
DEAL_II_VERSION=$DEALII_VERSION
EOF

# ===========================
# Install dependencies
# ===========================
echo "  ==========================================="
echo "  Installing deal.II and its dependencies ..."
./candi.sh --prefix="$DEALII_DIR" -j "$NPROC" -y

echo "  ==========================================="
echo "  deal.II $DEALII_VERSION build successful!"
echo "  Installed to: $DEALII_DIR/deal.II-$DEALII_VERSION"
