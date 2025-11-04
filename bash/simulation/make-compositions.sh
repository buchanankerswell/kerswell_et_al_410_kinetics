#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

# ===========================
# Usage / args
# ===========================
if [ $# -ne 6 ]; then
  echo "Usage: $0 <FIG_DIR> <OUT_DIR> <TIMESTEP> <ID1> <ID2> <ID3>"
  exit 1
fi

FIG_DIR="$1"
OUT_DIR="$2"
TIMESTEP="$3"
ID1="$4"
ID2="$5"
ID3="$6"
Y_CROP=5

mkdir -p "${OUT_DIR}"

ID1_DIR=${ID1//-/_}
ID2_DIR=${ID2//-/_}
ID3_DIR=${ID3//-/_}

echo "    --------------------------------------------------"
echo "    Composing mesh plots"
echo "    --------------------------------------------------"

# ===========================
# Helper function
# ===========================
compose_figures() {
  local TYPE="$1"
  local MODE="$2"
  local SUFFIXES=("${!3}")
  local OUT_NAME="$4"
  local OUT_PATH="${OUT_DIR}/${OUT_NAME}-${TIMESTEP}.png"

  if [[ -f "$OUT_PATH" ]]; then
    echo " -- Found composition: $OUT_PATH"
    return
  fi

  echo " -> Drawing composition: ${OUT_DIR}/${OUT_NAME}-${TIMESTEP}.png"

  if [[ "$MODE" == "single" ]]; then
    local ID="$5"
    local ID_DIR="${ID//-/_}"
    magick \
      \( "${FIG_DIR}/${TYPE}_${ID_DIR}/tiles/${TYPE}-${ID}-${SUFFIXES[0]}-${TIMESTEP}-tagged-abc.png" -gravity south -chop "0x${Y_CROP}%" \) \
      \( "${FIG_DIR}/${TYPE}_${ID_DIR}/tiles/${TYPE}-${ID}-${SUFFIXES[1]}-${TIMESTEP}-tagged-def.png" -gravity south -chop "0x${Y_CROP}%" \) \
      \( "${FIG_DIR}/${TYPE}_${ID_DIR}/tiles/${TYPE}-${ID}-${SUFFIXES[2]}-${TIMESTEP}-tagged-ghi.png" -gravity south -chop "0x${Y_CROP}%" \) \
      "${FIG_DIR}/${TYPE}_${ID_DIR}/tiles/${TYPE}-${ID}-${SUFFIXES[3]}-${TIMESTEP}-tagged-jkl.png" -append "${OUT_PATH}"

  elif [[ "$MODE" == "triple" ]]; then
    magick \
      \( "${FIG_DIR}/${TYPE}_${ID1_DIR}/tiles/${TYPE}-${ID1}-${SUFFIXES[0]}-${TIMESTEP}-tagged-abc.png" -gravity south -chop "0x${Y_CROP}%" \) \
      \( "${FIG_DIR}/${TYPE}_${ID2_DIR}/tiles/${TYPE}-${ID2}-${SUFFIXES[0]}-${TIMESTEP}-tagged-def.png" -gravity south -chop "0x${Y_CROP}%" \) \
      "${FIG_DIR}/${TYPE}_${ID3_DIR}/tiles/${TYPE}-${ID3}-${SUFFIXES[0]}-${TIMESTEP}-tagged-ghi.png" -append "${OUT_PATH}"
  fi
}

FULL_SET=(
  "temperature-nonadiabatic-pressure-nonadiabatic-density-nonadiabatic"
  "arrhenius-thermodynamic-reaction-rate"
  "X-field-vp-vs"
  "viscosity-sigma-ii-strain-rate"
)
SET1=("temperature-nonadiabatic-reaction-rate-X-field")
SET2=("temperature-nonadiabatic-density-nonadiabatic-vp")

compose_figures "slab" "single" FULL_SET[@] "slab-${ID1}-full-set-composition" "$ID1"
compose_figures "slab" "single" FULL_SET[@] "slab-${ID2}-full-set-composition" "$ID2"
compose_figures "slab" "single" FULL_SET[@] "slab-${ID3}-full-set-composition" "$ID3"

compose_figures "slab" "triple" SET1[@] "slab-${ID1}-${ID2}-${ID3}-set1-composition"
compose_figures "slab" "triple" SET2[@] "slab-${ID1}-${ID2}-${ID3}-set2-composition"

compose_figures "plume" "single" FULL_SET[@] "plume-${ID1}-full-set-composition" "$ID1"
compose_figures "plume" "single" FULL_SET[@] "plume-${ID2}-full-set-composition" "$ID2"
compose_figures "plume" "single" FULL_SET[@] "plume-${ID3}-full-set-composition" "$ID3"

compose_figures "plume" "triple" SET1[@] "plume-${ID1}-${ID2}-${ID3}-set1-composition"
compose_figures "plume" "triple" SET2[@] "plume-${ID1}-${ID2}-${ID3}-set2-composition"
