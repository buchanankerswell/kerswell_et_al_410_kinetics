# Erratum: correction to Kerswell et al. (2026)

**Re:** Correction of the units label on the kinetic prefactor $Z$

## Overview

In our manuscript (doi: 10.1029/2026JB033781), the kinetic prefactor $Z$ in the interface-controlled olivine $\Leftrightarrow$ wadsleyite reaction-rate law (Equation 13) is reported with units of K$^{-1}$ s$^{-1}$:

$$
  \dot{X} = Z\, T\, \exp\!\left(-\frac{H^{\ast} + P V^{\ast}}{R T}\right) \left(1 - \exp\!\left[-\frac{\Delta G}{R T}\right]\right) \left(1 - X\right)
$$

***The units label is wrong.*** Our ASPECT simulations are set to advance time in years, so the prefactor they apply is also in units of years: $Z$ = 3.0e0--7.0e7 **K$^{-1}$ yr$^{-1}$**, not in units of K$^{-1}$ s$^{-1}$ as originally published.

This is a units label correction only, not a correction of results. Since $Z$ scales the rate linearly, and every result in the paper is read directly from the simulation output, all quantitative results, the qualitative regime thresholds, and all figures and conclusions are correct as published. Relabeling $Z$ units does not change any reported values.

We have corrected the units label wherever it appears in the manuscript, Supplementary Information, and revised two passages in the Methods (+1 sentence) and Uncertainties (+3 sentences) sections.

One point of interpretation follows from the correction, and we have stated it clearly in the revised text. Per AGU's correction policy, we will also publish a footnote with the corrected article describing the change.

We apologize for the error and thank the editor for the opportunity to correct it.

Sincerely,

Buchanan Kerswell (on behalf of all coauthors)

\clearpage

## Correction Summary

The consequence of the mislabeled units can be quantified directly from our simulation outputs. If simulations are ran using the laboratory-derived $Z$ with mislabeled units, every phase transition would behave as quasi-equilibrium (transformation rates between 3 yr to 3 hr; Table 1). Correcting the units label avoids such errors and confusion when reproducing our results.

| $Z$ used in simulations (K$^{-1}$ yr$^{-1}$) | Published rate (Ma$^{-1}$) | $Z$ from mislabeled units (K$^{-1}$ yr$^{-1}$) | Rate from mislabeld units (Ma$^{-1}$) | Characteristic transformation time from mislabeled units |
|:------|:------|:-------|:------|:---------------|
| 3.0e0 | 0.012 | 9.5e7  | 7.4e1 | $\sim$ 13.7 yr |
| 2.0e2 | 0.068 | 6.3e9  | 4.9e3 | $\sim$  205 yr |
| 4.7e2 | 0.30  | 1.5e10 | 1.2e4 | $\sim$   87 yr |
| 1.8e5 | 2.2   | 5.7e12 | 4.4e6 | $\sim$  2.7 mo |
| 7.0e7 | 103   | 2.2e15 | 1.7e9 | $\sim$  5.1 hr |

Table: Approximate transformation time scales if simulation use $Z$ with mislabeled units (transformation time is 1/rate).

| Location | Change |
|:---------|:-------|
| Symbol table ($Z$) | K$^{-1}$ s$^{-1}$ $\rightarrow$ K$^{-1}$ yr$^{-1}$ |
| Symbol table (reaction rate) | s$^{-1}$ $\rightarrow$ yr$^{-1}$ |
| Methods, $Z$-range paragraph | Units relabeled; added one sentence |
| Methods, regime boundaries | K$^{-1}$ s$^{-1}$ $\rightarrow$ K$^{-1}$ yr$^{-1}$ |
| Discussion Uncertainties | Added three sentences |
| Figures 5 and 6 captions | K$^{-1}$ s$^{-1}$ $\rightarrow$ K$^{-1}$ yr$^{-1}$ |
| SI captions | K s$^{-1}$ $\rightarrow$ K yr$^{-1}$ |

Table: Complete list of changes

### Methods edit

> *The range of kinetic prefactors $Z$ used in our numerical experiments (3.0e0--7.0e7 **K$^{-1}$ yr$^{-1}$)** was determined by holding $\Gamma$ = $\exp\!\left(-18\right)$ m s$^{-1}$ K$^{-1}$ ppm$_\mathrm{OH}^{-n}$, $H^{\ast}$ = 274 kJ mol$^{-1}$, $V^{\ast}$ = 3.0e-6 m$^3$ mol$^{-1}$, and $n$ = 3.2 constant, while varying water content $C_\mathrm{OH}$ from 50--5000 ppm and grain size $d$ from 1--10 mm. These water contents and grain sizes are consistent with the experimental conditions of @hosoya2005, previous numerical studies of metastable olivine wedges [@rubie1994], and typical grain sizes of upper mantle xenoliths [$\sim$ 3--10 mm, @karato1984; @karato2008]. **Because the rate law of @hosoya2005 is calibrated per second, realizing these conditions at geological time scales requires extrapolating the growth-rate prefactor $\Gamma$ by approximately 7 orders of magnitude.** Thus, our experiments approximate kinetic conditions ranging from slow kinetics in dry rocks with large grain sizes (50 ppm OH; 10 mm; $Z$ = 3.0e0 **K$^{-1}$ yr$^{-1}$**) to fast kinetics in hydrated rocks with small grain sizes (5000 ppm OH; 1 mm; $Z$ = 7.0e7 **K$^{-1}$ yr$^{-1}$**).*

### Uncertainties edit

> *The primary quantitative uncertainty in our analysis stems from the kinetic prefactor $Z$ in the interface-controlled olivine $\Leftrightarrow$ wadsleyite growth model (Equation 13), which spans several orders of magnitude reflecting variable water contents (50--5000 ppm) and grain sizes (1--10 mm). Water content and grain size were intentionally absorbed into the single kinetic prefactor $Z$ in our formulation (Section 2.2.3.2), allowing characterization of first-order kinetic behavior in terms of "fast versus slow" kinetics without presupposing specific water or grain-size distributions. While the OH-dependence of the transformation rate [@hosoya2005] and the thermodynamic effect of water on the phase boundary [@smyth1987; @smyth2002] both motivate explicit treatment of water content, mantle water concentrations and their spatial distribution within subducted slabs remain uncertain [@houser2016; @karato2011; @ishii2021; @cerpa2022]. Laboratory studies also reveal large uncertainties in kinetic parameters $n$, $\Gamma$, $H^{\ast}$, and $V^{\ast}$ that depend strongly on water content, grain size, Mg-Fe composition, and microstructural evolution [@rubie1994; @liu1998; @kubo2004; @hosoya2005; @perrillat2013; @ledoux2023]. **Moreover, using $Z$ values derived directly from the kinetic parameters reported in @hosoya2005 produce conditions where olivine transforms instantaneously on geological time scales, holding the 410 at thermodynamic equilibrium. Testing the hypothesis that kinetics control 410 structure and slab descent requires extrapolating kinetic conditions by approximately 7 orders of magnitude from experimental to geological rates**. Our simulations therefore explore only a limited subset of potential reaction rates in Earth's upper mantle, and we interpret our results within the conceptual framework of "fast versus slow" kinetics rather than "wet versus dry" conditions and/or "fine versus coarse" grain sizes. **These uncertainties motivate new experiments to constrain reaction rates at relevant mantle conditions.***

### Footnote for the corrected version of record

> **Correction:** our ASPECT simulations are set to advance time in years, so the prefactor they apply is also in units of years: $Z$ = 3.0e0--7.0e7 **K$^{-1}$ yr$^{-1}$**, not in units of K$^{-1}$ s$^{-1}$ as originally published. The units label has been corrected in the version of record.
>
> **Implication:** the kinetic regimes we describe, and especially the ultra-sluggish regime, emerge only because the simulated kinetic conditions are $\sim$ 7 orders of magnitude slower than laboratory conditions. Realizing these geological rates requires extrapolating the growth-rate prefactor $\Gamma$ by roughly 7 orders of magnitude from the experimental conditions reported by @hosoya2005. This time-scale extrapolation is a first-order uncertainty in our results and motivates new experiments at relevant mantle conditions.


\clearpage

# References {.unnumbered #sec:references}

::: {#refs}
:::
