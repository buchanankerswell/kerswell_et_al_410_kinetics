#!/usr/bin/env bash
set -e

echo "    --------------------------------------------------"
echo "    Drawing composition plots"
echo "    --------------------------------------------------"

# Full set slow
echo " -> [slab] full set slow"
magick ../../figs/simulation/2d_box/meshes/slab_lnk18_Ha274_D1e12_060ppm/tiles/slab-lnk18-Ha274-D1e12-060ppm-temperature-nonadiabatic-pressure-nonadiabatic-density-nonadiabatic-0100-tagged-abc.png ../../figs/simulation/2d_box/meshes/slab_lnk18_Ha274_D1e12_060ppm/tiles/slab-lnk18-Ha274-D1e12-060ppm-thermodynamic-growth-rate-reaction-rate-0100-tagged-def.png ../../figs/simulation/2d_box/meshes/slab_lnk18_Ha274_D1e12_060ppm/tiles/slab-lnk18-Ha274-D1e12-060ppm-X-field-vp-vs-0100-tagged-ghi.png -append ../../figs/simulation/2d_box/compositions/slab-lnk18-Ha274-D1e12-060ppm-full-set-composition-0100.png

# Full set default
echo " -> [slab] full set default"
magick ../../figs/simulation/2d_box/meshes/slab_lnk18_Ha274_d5mm_060ppm/tiles/slab-lnk18-Ha274-d5mm-060ppm-temperature-nonadiabatic-pressure-nonadiabatic-density-nonadiabatic-0100-tagged-abc.png ../../figs/simulation/2d_box/meshes/slab_lnk18_Ha274_d5mm_060ppm/tiles/slab-lnk18-Ha274-d5mm-060ppm-thermodynamic-growth-rate-reaction-rate-0100-tagged-def.png ../../figs/simulation/2d_box/meshes/slab_lnk18_Ha274_d5mm_060ppm/tiles/slab-lnk18-Ha274-d5mm-060ppm-X-field-vp-vs-0100-tagged-ghi.png -append ../../figs/simulation/2d_box/compositions/slab-lnk18-Ha274-d5mm-060ppm-full-set-composition-0100.png

# Full set fast
echo " -> [slab] full set fast"
magick ../../figs/simulation/2d_box/meshes/slab_lnk18_Ha300_d5mm_060ppm/tiles/slab-lnk18-Ha300-d5mm-060ppm-temperature-nonadiabatic-pressure-nonadiabatic-density-nonadiabatic-0100-tagged-abc.png ../../figs/simulation/2d_box/meshes/slab_lnk18_Ha300_d5mm_060ppm/tiles/slab-lnk18-Ha300-d5mm-060ppm-thermodynamic-growth-rate-reaction-rate-0100-tagged-def.png ../../figs/simulation/2d_box/meshes/slab_lnk18_Ha300_d5mm_060ppm/tiles/slab-lnk18-Ha300-d5mm-060ppm-X-field-vp-vs-0100-tagged-ghi.png -append ../../figs/simulation/2d_box/compositions/slab-lnk18-Ha300-d5mm-060ppm-full-set-composition-0100.png

# Slow, default, fast set1
echo " -> [slab] slow–fast set1"
magick ../../figs/simulation/2d_box/meshes/slab_lnk18_Ha300_d5mm_060ppm/tiles/slab-lnk18-Ha300-d5mm-060ppm-temperature-nonadiabatic-reaction-rate-X-field-0100-tagged-abc.png ../../figs/simulation/2d_box/meshes/slab_lnk18_Ha274_d5mm_060ppm/tiles/slab-lnk18-Ha274-d5mm-060ppm-temperature-nonadiabatic-reaction-rate-X-field-0100-tagged-def.png ../../figs/simulation/2d_box/meshes/slab_lnk18_Ha274_D1e12_060ppm/tiles/slab-lnk18-Ha274-D1e12-060ppm-temperature-nonadiabatic-reaction-rate-X-field-0100-tagged-ghi.png -append ../../figs/simulation/2d_box/compositions/slab-slow-default-fast-set1-composition-0100.png

# Slow, default, fast set2
echo " -> [slab] slow–fast set2"
magick ../../figs/simulation/2d_box/meshes/slab_lnk18_Ha300_d5mm_060ppm/tiles/slab-lnk18-Ha300-d5mm-060ppm-temperature-nonadiabatic-density-nonadiabatic-vp-0100-tagged-abc.png ../../figs/simulation/2d_box/meshes/slab_lnk18_Ha274_d5mm_060ppm/tiles/slab-lnk18-Ha274-d5mm-060ppm-temperature-nonadiabatic-density-nonadiabatic-vp-0100-tagged-def.png ../../figs/simulation/2d_box/meshes/slab_lnk18_Ha274_D1e12_060ppm/tiles/slab-lnk18-Ha274-D1e12-060ppm-temperature-nonadiabatic-density-nonadiabatic-vp-0100-tagged-ghi.png -append ../../figs/simulation/2d_box/compositions/slab-slow-default-fast-set2-composition-0100.png


# Full set slow
echo " -> [plume] full set slow"
magick ../../figs/simulation/2d_box/meshes/plume_lnk18_Ha274_D1e12_060ppm/tiles/plume-lnk18-Ha274-D1e12-060ppm-temperature-nonadiabatic-pressure-nonadiabatic-density-nonadiabatic-0100-tagged-abc.png ../../figs/simulation/2d_box/meshes/plume_lnk18_Ha274_D1e12_060ppm/tiles/plume-lnk18-Ha274-D1e12-060ppm-thermodynamic-growth-rate-reaction-rate-0100-tagged-def.png ../../figs/simulation/2d_box/meshes/plume_lnk18_Ha274_D1e12_060ppm/tiles/plume-lnk18-Ha274-D1e12-060ppm-X-field-vp-vs-0100-tagged-ghi.png -append ../../figs/simulation/2d_box/compositions/plume-lnk18-Ha274-D1e12-060ppm-full-set-composition-0100.png

# Full set default
echo " -> [plume] full set default"
magick ../../figs/simulation/2d_box/meshes/plume_lnk18_Ha274_d5mm_060ppm/tiles/plume-lnk18-Ha274-d5mm-060ppm-temperature-nonadiabatic-pressure-nonadiabatic-density-nonadiabatic-0100-tagged-abc.png ../../figs/simulation/2d_box/meshes/plume_lnk18_Ha274_d5mm_060ppm/tiles/plume-lnk18-Ha274-d5mm-060ppm-thermodynamic-growth-rate-reaction-rate-0100-tagged-def.png ../../figs/simulation/2d_box/meshes/plume_lnk18_Ha274_d5mm_060ppm/tiles/plume-lnk18-Ha274-d5mm-060ppm-X-field-vp-vs-0100-tagged-ghi.png -append ../../figs/simulation/2d_box/compositions/plume-lnk18-Ha274-d5mm-060ppm-full-set-composition-0100.png

# Full set fast
echo " -> [plume] full set fast"
magick ../../figs/simulation/2d_box/meshes/plume_lnk18_Ha300_d5mm_060ppm/tiles/plume-lnk18-Ha300-d5mm-060ppm-temperature-nonadiabatic-pressure-nonadiabatic-density-nonadiabatic-0100-tagged-abc.png ../../figs/simulation/2d_box/meshes/plume_lnk18_Ha300_d5mm_060ppm/tiles/plume-lnk18-Ha300-d5mm-060ppm-thermodynamic-growth-rate-reaction-rate-0100-tagged-def.png ../../figs/simulation/2d_box/meshes/plume_lnk18_Ha300_d5mm_060ppm/tiles/plume-lnk18-Ha300-d5mm-060ppm-X-field-vp-vs-0100-tagged-ghi.png -append ../../figs/simulation/2d_box/compositions/plume-lnk18-Ha300-d5mm-060ppm-full-set-composition-0100.png

# Slow, default, fast set
echo " -> [plume] slow–fast set1"
magick ../../figs/simulation/2d_box/meshes/plume_lnk18_Ha300_d5mm_060ppm/tiles/plume-lnk18-Ha300-d5mm-060ppm-temperature-nonadiabatic-reaction-rate-X-field-0100-tagged-abc.png ../../figs/simulation/2d_box/meshes/plume_lnk18_Ha274_d5mm_060ppm/tiles/plume-lnk18-Ha274-d5mm-060ppm-temperature-nonadiabatic-reaction-rate-X-field-0100-tagged-def.png ../../figs/simulation/2d_box/meshes/plume_lnk18_Ha274_D1e12_060ppm/tiles/plume-lnk18-Ha274-D1e12-060ppm-temperature-nonadiabatic-reaction-rate-X-field-0100-tagged-ghi.png -append ../../figs/simulation/2d_box/compositions/plume-slow-default-fast-set1-composition-0100.png

# Slow, default, fast set2
echo " -> [plume] slow–fast set2"
magick ../../figs/simulation/2d_box/meshes/plume_lnk18_Ha300_d5mm_060ppm/tiles/plume-lnk18-Ha300-d5mm-060ppm-temperature-nonadiabatic-density-nonadiabatic-vp-0100-tagged-abc.png ../../figs/simulation/2d_box/meshes/plume_lnk18_Ha274_d5mm_060ppm/tiles/plume-lnk18-Ha274-d5mm-060ppm-temperature-nonadiabatic-density-nonadiabatic-vp-0100-tagged-def.png ../../figs/simulation/2d_box/meshes/plume_lnk18_Ha274_D1e12_060ppm/tiles/plume-lnk18-Ha274-D1e12-060ppm-temperature-nonadiabatic-density-nonadiabatic-vp-0100-tagged-ghi.png -append ../../figs/simulation/2d_box/compositions/plume-slow-default-fast-set2-composition-0100.png

