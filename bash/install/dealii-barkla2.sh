#!/usr/bin/env bash
set -eo pipefail
IFS=$'\n\t'

# ===========================
# Get args
# ===========================
if [ $# -ne 6 ]; then
  echo "Usage: $0 <DEALII_VERSION> <TRILINOS_VERSION> <GCC> <OPENMPI> <OPENBLAS> <CMAKE>"
  exit 1
fi

DEALII_VERSION=$1
TRILINOS_VERSION=$2
GCC=$3
OPENMPI=$4
OPENBLAS=$5
CMAKE=$6

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
NPROC=$SLURM_NTASKS
DEALII_DIR=$PROJECT_ROOT/dealii

echo "  PROJCT_ROOT:      $PROJECT_ROOT"
echo "  NPROC:            $NPROC"
echo "  DEALII_VERSION:   $DEALII_VERSION"
echo "  DEALII_DIR:       $DEALII_DIR"
echo "  GCC:              $GCC"
echo "  OPENMPI:          $OPENMPI"
echo "  OPENBLAS:         $OPENBLAS"
echo "  CMAKE:            $CMAKE"

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
# Load modules
# ===========================
echo "  ==========================================="
echo "  Loading openmpi and openblas modules ..."

module load "$GCC" "$OPENMPI" "$OPENBLAS" "$CMAKE"

export CC=mpicc
export CXX=mpicxx
export FC=mpif90
export FF=mpif77

# ===========================
# Clone candi
# ===========================
echo "  ==========================================="
TMP_CANDI_DIR=$(mktemp -d)
echo "  Cloning candi into tmp dir: $TMP_CANDI_DIR"
git clone https://github.com/dealii/candi.git "$TMP_CANDI_DIR" > /dev/null 2>&1
cd "$TMP_CANDI_DIR"
curl -O https://raw.githubusercontent.com/geodynamics/aspect/main/contrib/install/local.cfg
cat << EOF >> local.cfg
GIVEN_PLATFORM=deal.II-toolchain/platforms/supported/linux_cluster.platform
DEAL_II_VERSION=$DEALII_VERSION
TRILINOS_MAJOR_VERSION=AUTO
BLAS_DIR=$OPENBLASDIR
TRILINOS_CONFOPTS="-D TPL_BLAS_LIBRARIES=$OPENBLASLIB/libopenblas.so -D TPL_LAPACK_LIBRARIES=$OPENBLASLIB/libopenblas.so"
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
