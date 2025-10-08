#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

# ===========================
# Get args
# ===========================
if [ $# -ne 6 ]; then
  echo "Usage: $0 <FIG_DIR> <OUT_DIR> <TIMESTEP> <SLOW_FACTOR> <MID_FACTOR> <FAST_FACTOR>"
  exit 1
fi

FIG_DIR="$1"
OUT_DIR="$2"
TIMESTEP="$3"
SLOW_FACTOR="$4"
MID_FACTOR="$5"
FAST_FACTOR="$6"

Y_CROP=5

mkdir -p "${OUT_DIR}"

echo "    --------------------------------------------------"
echo "    Composing mesh plots"
echo "    --------------------------------------------------"

# Full set slow
echo " -> [slab] full set slow"
magick \
  \( "${FIG_DIR}/slab_${SLOW_FACTOR}/tiles/slab-${SLOW_FACTOR}-temperature-nonadiabatic-pressure-nonadiabatic-density-nonadiabatic-${TIMESTEP}-tagged-abc.png" \
     -gravity south -chop "0x${Y_CROP}%" \) \
  \( "${FIG_DIR}/slab_${SLOW_FACTOR}/tiles/slab-${SLOW_FACTOR}-arrhenius-thermodynamic-reaction-rate-${TIMESTEP}-tagged-def.png" \
     -gravity south -chop "0x${Y_CROP}%" \) \
  "${FIG_DIR}/slab_${SLOW_FACTOR}/tiles/slab-${SLOW_FACTOR}-X-field-vp-vs-${TIMESTEP}-tagged-ghi.png" \
  -append "${OUT_DIR}/slab-${SLOW_FACTOR}-full-set-composition-${TIMESTEP}.png"

# Full set default
echo " -> [slab] full set default"
magick \
  \( "${FIG_DIR}/slab_${MID_FACTOR}/tiles/slab-${MID_FACTOR}-temperature-nonadiabatic-pressure-nonadiabatic-density-nonadiabatic-${TIMESTEP}-tagged-abc.png" \
     -gravity south -chop "0x${Y_CROP}%" \) \
  \( "${FIG_DIR}/slab_${MID_FACTOR}/tiles/slab-${MID_FACTOR}-arrhenius-thermodynamic-reaction-rate-${TIMESTEP}-tagged-def.png" \
     -gravity south -chop "0x${Y_CROP}%" \) \
  "${FIG_DIR}/slab_${MID_FACTOR}/tiles/slab-${MID_FACTOR}-X-field-vp-vs-${TIMESTEP}-tagged-ghi.png" \
  -append "${OUT_DIR}/slab-${MID_FACTOR}-full-set-composition-${TIMESTEP}.png"

# Full set fast
echo " -> [slab] full set fast"
magick \
  \( "${FIG_DIR}/slab_${FAST_FACTOR}/tiles/slab-${FAST_FACTOR}-temperature-nonadiabatic-pressure-nonadiabatic-density-nonadiabatic-${TIMESTEP}-tagged-abc.png" \
     -gravity south -chop "0x${Y_CROP}%" \) \
  \( "${FIG_DIR}/slab_${FAST_FACTOR}/tiles/slab-${FAST_FACTOR}-arrhenius-thermodynamic-reaction-rate-${TIMESTEP}-tagged-def.png" \
     -gravity south -chop "0x${Y_CROP}%" \) \
  "${FIG_DIR}/slab_${FAST_FACTOR}/tiles/slab-${FAST_FACTOR}-X-field-vp-vs-${TIMESTEP}-tagged-ghi.png" \
  -append "${OUT_DIR}/slab-${FAST_FACTOR}-full-set-composition-${TIMESTEP}.png"

# Slow, default, fast set1
echo " -> [slab] slow–fast set1"
magick \
  \( "${FIG_DIR}/slab_${SLOW_FACTOR}/tiles/slab-${SLOW_FACTOR}-temperature-nonadiabatic-reaction-rate-X-field-${TIMESTEP}-tagged-abc.png" \
     -gravity south -chop "0x${Y_CROP}%" \) \
  \( "${FIG_DIR}/slab_${MID_FACTOR}/tiles/slab-${MID_FACTOR}-temperature-nonadiabatic-reaction-rate-X-field-${TIMESTEP}-tagged-def.png" \
     -gravity south -chop "0x${Y_CROP}%" \) \
  "${FIG_DIR}/slab_${FAST_FACTOR}/tiles/slab-${FAST_FACTOR}-temperature-nonadiabatic-reaction-rate-X-field-${TIMESTEP}-tagged-ghi.png" \
  -append "${OUT_DIR}/slab-${SLOW_FACTOR}-${MID_FACTOR}-${FAST_FACTOR}-set1-composition-${TIMESTEP}.png"

# Slow, default, fast set2
echo " -> [slab] slow–fast set2"
magick \
  \( "${FIG_DIR}/slab_${SLOW_FACTOR}/tiles/slab-${SLOW_FACTOR}-temperature-nonadiabatic-density-nonadiabatic-vp-${TIMESTEP}-tagged-abc.png" \
     -gravity south -chop "0x${Y_CROP}%" \) \
  \( "${FIG_DIR}/slab_${MID_FACTOR}/tiles/slab-${MID_FACTOR}-temperature-nonadiabatic-density-nonadiabatic-vp-${TIMESTEP}-tagged-def.png" \
     -gravity south -chop "0x${Y_CROP}%" \) \
  "${FIG_DIR}/slab_${FAST_FACTOR}/tiles/slab-${FAST_FACTOR}-temperature-nonadiabatic-density-nonadiabatic-vp-${TIMESTEP}-tagged-ghi.png" \
  -append "${OUT_DIR}/slab-${SLOW_FACTOR}-${MID_FACTOR}-${FAST_FACTOR}-set2-composition-${TIMESTEP}.png"

# Full set slow
echo " -> [plume] full set slow"
magick \
  \( "${FIG_DIR}/plume_${SLOW_FACTOR}/tiles/plume-${SLOW_FACTOR}-temperature-nonadiabatic-pressure-nonadiabatic-density-nonadiabatic-${TIMESTEP}-tagged-abc.png" \
     -gravity south -chop "0x${Y_CROP}%" \) \
  \( "${FIG_DIR}/plume_${SLOW_FACTOR}/tiles/plume-${SLOW_FACTOR}-arrhenius-thermodynamic-reaction-rate-${TIMESTEP}-tagged-def.png" \
     -gravity south -chop "0x${Y_CROP}%" \) \
  "${FIG_DIR}/plume_${SLOW_FACTOR}/tiles/plume-${SLOW_FACTOR}-X-field-vp-vs-${TIMESTEP}-tagged-ghi.png" \
  -append "${OUT_DIR}/plume-${SLOW_FACTOR}-full-set-composition-${TIMESTEP}.png"

# Full set default
echo " -> [plume] full set default"
magick \
  \( "${FIG_DIR}/plume_${MID_FACTOR}/tiles/plume-${MID_FACTOR}-temperature-nonadiabatic-pressure-nonadiabatic-density-nonadiabatic-${TIMESTEP}-tagged-abc.png" \
     -gravity south -chop "0x${Y_CROP}%" \) \
  \( "${FIG_DIR}/plume_${MID_FACTOR}/tiles/plume-${MID_FACTOR}-arrhenius-thermodynamic-reaction-rate-${TIMESTEP}-tagged-def.png" \
     -gravity south -chop "0x${Y_CROP}%" \) \
  "${FIG_DIR}/plume_${MID_FACTOR}/tiles/plume-${MID_FACTOR}-X-field-vp-vs-${TIMESTEP}-tagged-ghi.png" \
  -append "${OUT_DIR}/plume-${MID_FACTOR}-full-set-composition-${TIMESTEP}.png"

# Full set fast
echo " -> [plume] full set fast"
magick \
  \( "${FIG_DIR}/plume_${FAST_FACTOR}/tiles/plume-${FAST_FACTOR}-temperature-nonadiabatic-pressure-nonadiabatic-density-nonadiabatic-${TIMESTEP}-tagged-abc.png" \
     -gravity south -chop "0x${Y_CROP}%" \) \
  \( "${FIG_DIR}/plume_${FAST_FACTOR}/tiles/plume-${FAST_FACTOR}-arrhenius-thermodynamic-reaction-rate-${TIMESTEP}-tagged-def.png" \
     -gravity south -chop "0x${Y_CROP}%" \) \
  "${FIG_DIR}/plume_${FAST_FACTOR}/tiles/plume-${FAST_FACTOR}-X-field-vp-vs-${TIMESTEP}-tagged-ghi.png" \
  -append "${OUT_DIR}/plume-${FAST_FACTOR}-full-set-composition-${TIMESTEP}.png"

# Slow, default, fast set1
echo " -> [plume] slow–fast set1"
magick \
  \( "${FIG_DIR}/plume_${SLOW_FACTOR}/tiles/plume-${SLOW_FACTOR}-temperature-nonadiabatic-reaction-rate-X-field-${TIMESTEP}-tagged-abc.png" \
     -gravity south -chop "0x${Y_CROP}%" \) \
  \( "${FIG_DIR}/plume_${MID_FACTOR}/tiles/plume-${MID_FACTOR}-temperature-nonadiabatic-reaction-rate-X-field-${TIMESTEP}-tagged-def.png" \
     -gravity south -chop "0x${Y_CROP}%" \) \
  "${FIG_DIR}/plume_${FAST_FACTOR}/tiles/plume-${FAST_FACTOR}-temperature-nonadiabatic-reaction-rate-X-field-${TIMESTEP}-tagged-ghi.png" \
  -append "${OUT_DIR}/plume-${SLOW_FACTOR}-${MID_FACTOR}-${FAST_FACTOR}-set1-composition-${TIMESTEP}.png"

# Slow, default, fast set2
echo " -> [plume] slow–fast set2"
magick \
  \( "${FIG_DIR}/plume_${SLOW_FACTOR}/tiles/plume-${SLOW_FACTOR}-temperature-nonadiabatic-density-nonadiabatic-vp-${TIMESTEP}-tagged-abc.png" \
     -gravity south -chop "0x${Y_CROP}%" \) \
  \( "${FIG_DIR}/plume_${MID_FACTOR}/tiles/plume-${MID_FACTOR}-temperature-nonadiabatic-density-nonadiabatic-vp-${TIMESTEP}-tagged-def.png" \
     -gravity south -chop "0x${Y_CROP}%" \) \
  "${FIG_DIR}/plume_${FAST_FACTOR}/tiles/plume-${FAST_FACTOR}-temperature-nonadiabatic-density-nonadiabatic-vp-${TIMESTEP}-tagged-ghi.png" \
  -append "${OUT_DIR}/plume-${SLOW_FACTOR}-${MID_FACTOR}-${FAST_FACTOR}-set2-composition-${TIMESTEP}.png"
