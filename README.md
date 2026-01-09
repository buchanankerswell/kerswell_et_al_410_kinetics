![](draft/repo-banner.png)

***Figure:*** *Measured 410 displacement and width versus maximum reaction rates $`\dot{X}_{\mathrm{max}}`$ in plume and slab simulations with intermediate strength contrasts ($`B`$ = 4) after 100 Ma. Structure of the 410 near plumes (left column) shows minimal dependence on $`\dot{X}_{\mathrm{max}}`$, with both displacement and width remaining nearly constant across seven orders of magnitude variation in $`\dot{X}_{\mathrm{max}}`$. In contrast, 410 structure near slabs (right column) changes distinctly across three kinetic regimes: (1) quasi-equilibrium at high $`\dot{X}_{\mathrm{max}}`$, where 410 widths are narrow and elevated; (2) an intermediate regime where decreasing reaction rates $`\dot{X}_{\mathrm{max}}`$ progressively widen and deepen the 410; and (3) an ultra-sluggish regime at low $`\dot{X}_{\mathrm{max}}`$, where the 410 narrows while deepening, and slabs completely stall and pond.*

# Kerswell et al. (2026; JGR: Solid Earth)

## Repository

This repository provides all materials for the manuscript *Beyond Equilibrium: Kinetic Thresholds and Rheological Feedbacks Create a Potentially Complex 410 in Slab Regions* (Kerswell et al., 2026; in revision at JGR: Solid Earth), including all datasets required to compile the study and scripts to reproduce all results and figures.

## Prerequisite software

### Python

I recommend installing the [miniforge](https://github.com/conda-forge/miniforge) python distribution. This includes a minimal installation of python (plus some dependencies) and the package manager [conda](https://docs.conda.io/en/latest/), which is required to build the necessary python environment for this study.

### R

R is a programming language used to visualize the results in this study. R can be downloaded from the [R Project homepage](https://www.r-project.org). Follow their instructions to install R on your machine.

### Pandoc

Pandoc is a universal document converter used to build a PDF version of the manuscript, which written in Markdown. Pandoc can be downloaded from the [Pandoc homepage](https://pandoc.org). Follow their instructions to install Pandoc on your machine.

## Reproducing the study

The full set of ASPECT solution files and 410 structural measurement dataset can be downloaded from the [Open Science Framework repository](https://osf.io/9phwc/files). The following steps will download all necessary data, reproduce the results and figures, and build a pdf version of the manuscript:

```bash
# Clone this repository
git clone https://github.com/buchanankerswell/kerswell_et_al_410_kinetics.git

# Change into the directory
cd kerswell_et_al_410_kinetics

# Get data and reproduce figures
make build
```

## Coauthors

- [John Wheeler](https://scholar.google.com/citations?user=jsfp2-8AAAAJ&hl=en&oi=ao) (Department of Earth, Oceans and Ecological Sciences, University of Liverpool)
- [Rene Gassmöller](https://scholar.google.com/citations?user=Vk8SmssAAAAJ&hl=en&oi=ao) (GEOMAR Helmholtz Centre for Ocean Research Kiel)
- [J. Huw Davies](https://scholar.google.com/citations?user=T5ygdwcAAAAJ&hl=en&oi=ao) (School of Earth and Environmental Sciences, Cardiff University)
- Isabel Papanagnou (Bullard Laboratories, Department of Earth Sciences, University of Cambridge)
- [Sanne Cottaar](https://scholar.google.com/citations?user=l5JtmzkAAAAJ&hl=en&oi=ao) (Bullard Laboratories, Department of Earth Sciences, University of Cambridge)

## Acknowledgement

This work was funded by the UKRI NERC Large Grant no. NE/V018477/1. All computations were undertaken on Barkla2, part of the High Performance Computing facilities at the University of Liverpool, who graciously provided expert support. We thank the Computational Infrastructure for Geodynamics ([https://geodynamics.org](https://geodynamics.org)) which is funded by the National Science Foundation under award EAR-0949446 and EAR-1550901 for supporting the development of ASPECT. We also extend our gratitude towards Sujoy Ghosh for his editorial handling, and two anonymous reviewers for their constructive feedback that improved the manuscript.

## Data Availability

All data, code, and relevant information for reproducing this work can be found at [https://github.com/buchanankerswell/kerswell_et_al_410_kinetics](https://github.com/buchanankerswell/kerswell_et_al_410_kinetics), and at [https://doi.org/10.17605/OSF.IO/9PHWC](https://doi.org/10.17605/OSF.IO/9PHWC), the official Open Science Framework data repository. All code within these repositories is MIT Licensed and free for use and distribution (see license details). ASPECT version 3.0.0, (Bangerth et al., 2024a, 2024b; Clevenger & Heister, 2021; Fraters et al., 2019; Fraters, 2020; Gassmöller et al., 2018; Heister et al., 2017; Kronbichler et al., 2012) used in these computations is freely available under the GPL v2.0 or later license through its software landing page [https://geodynamics.org/resources/aspect](https://geodynamics.org/resources/aspect) or [https://aspect.geodynamics.org](https://aspect.geodynamics.org) and is being actively developed on GitHub and can be accessed via [https://github.com/geodynamics/aspect](https://github.com/geodynamics/aspect).

## Abstract

The seismic expression of Earth's 410 km discontinuity varies across tectonic settings, from sharp, high-amplitude interfaces to broad transitions---patterns that cannot be explained by equilibrium thermodynamics without invoking large-scale thermal or compositional heterogeneities. Laboratory experiments show the olivine $`\Leftrightarrow`$ wadsleyite transition responsible for the 410 is rate-limited, yet previous numerical studies have not directly evaluated the sensitivity of 410 structure to kinetic and rheological factors. Here we investigate these relationships by coupling a grain-scale, interface-controlled olivine $`\Leftrightarrow`$ wadsleyite growth model to compressible simulations of mantle plumes and subducting slabs. We vary kinetic parameters across seven orders of magnitude and quantify the resulting 410 displacements and widths. Our results reveal an asymmetry between hot and cold environments. In plumes, high temperatures produce sharp 410s (2--3 km wide) regardless of kinetics. In slabs, kinetics exert first-order control on 410 structure through three regimes: (1) quasi-equilibrium conditions producing narrow, uplifted 410s and continuous slab descent; (2) intermediate reaction rates generating broader, deeper 410s with metastable olivine wedges resisting downward slab motion; and (3) ultra-sluggish reaction rates causing slab stagnation with re-sharpened, deeply displaced 410s ($`\lesssim`$ 100 km). Rheological contrasts modulate these kinetic effects by controlling slab geometry and residence time in the phase transition zone. These findings demonstrate that reaction rates strongly influence 410 structure in subduction zones, establishing the 410 as a potential seismological constraint on upper mantle kinetic processes, particularly in cold environments where disequilibrium effects are amplified.

# License

MIT License

Copyright (c) 2024 Buchanan Kerswell

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
