# Erratum: correction to Kerswell et al. (2026)

**Re:** Correction of the kinetic prefactor values

## Overview

In our manuscript (doi: 10.1029/2026JB033781), the kinetic prefactor $Z$ in the interface-controlled olivine $\Leftrightarrow$ wadsleyite rate law (Equation 13 in the main text) is reported with units of K$^{-1}$ s$^{-1}$:

$$
  \dot{X} = Z\, T\, \exp\!\left(-\frac{H^{\ast} + P V^{\ast}}{R T}\right) \left(1 - \exp\!\left[-\frac{\Delta G}{R T}\right]\right) \left(1 - X\right)
$$ {#eq:reaction-rate}

Because our ASPECT simulations were set to output time in years, the prefactor range we reported in the manuscript was inadvertently shifted by a factor of 3.15e7 (seconds in a year). Rather than $Z$ = 3.0e0--7.0e7 K$^{-1}$ s$^{-1}$ as originally published, the $Z$ values actually used by ASPECT were 9.5e-8--2.2e0 K$^{-1}$ s$^{-1}$. This shift was an accidental user error and not a bug in the ASPECT code; we discovered the issue when testing subsequent models.

We have corrected the $Z$ values wherever they appear in the manuscript and Supplementary Information. To contextualize the correction in terms of experimental constraints, we also revised one passage in Section 2.2.3.2 (Reaction Kinetics methodology) and another passage in Section 4.1 (Uncertainties and Limitations discussion). Per AGU's correction policy, we will also publish a footnote with the corrected article describing the change.

Since every result in the paper is read directly from the simulation output, all quantitative results, qualitative regime thresholds, figures, and conclusions are correct as published. This erratum only clarifies what $Z$ values were actually used in our simulations.

We apologize for the error and thank the editor for the opportunity to correct it.

Sincerely,

Buchanan Kerswell (on behalf of all coauthors)

\clearpage

## Correction Summary

Assuming that there is plenty of reactant ($X$ = 0) and $\Delta G$ of the reaction is sufficiently large such that the driving force term $\left(1 - \exp\!\left[-\frac{\Delta G}{R T}\right]\right)$ is equal to one, Erratum Eq. -@eq:reaction-rate simplifies to:

$$
  \dot{X} = Z\, T\, \exp\!\left(-\frac{H^{\ast} + P V^{\ast}}{R T}\right)
$$ {#eq:reaction-rate-simple}

The consequence of using the uncorrected $Z$ range can be approximated by applying ambient mantle conditions to the above equation. Below is a table that summarizes the different reaction rates $\dot{X}$ (Ma$^{-1}$) and characteristic transformation times $\tau$ = 1/$\dot{X}$ (Ma) one could roughly expect for the corrected and uncorrected $Z$ range values (K$^{-1}$ s$^{-1}$), evaluated at $P$ = 10 GPa, $T$ = 1706 K, $H^{\ast}$ = 274 kJ mol$^{-1}$, $V^{\ast}$ = 3.0e-6 m$^3$ mol$^{-1}$.

| $Z$ | $\dot{X}$ | $\tau$ | $\widetilde{Z}$ | $\widetilde{\dot{X}}$ | $\widetilde{\tau}$ |
|:------|:-------|:--------|:-------|:------|:-------|
| 3.0e0 | 7.9e7  | 1.3e-8  | 9.5e-8 | 2.5   | 4.0e-1 |
| 2.0e2 | 5.3e9  | 1.9e-10 | 6.3e-6 | 1.7e2 | 6.0e-3 |
| 1.8e5 | 4.8e12 | 2.1e-13 | 5.7e-3 | 1.5e5 | 6.6e-6 |
| 7.0e7 | 1.9e15 | 5.4e-16 | 2.2e0  | 5.8e7 | 1.7e-8 |

Table: Approximate reaction rates and transformation times using the published values (first three columns) versus the corrected values (marked with a tilde). {#tbl:simple-calculation}

Running ASPECT simulations with the published $Z$ values results in characteristic transformation times between 17 milliseconds and 5 days (i.e., instantaneous, quasi-equilibrium transformations at ambient mantle conditions). Using the corrected values reproduces the results in the original publication, which are qualitatively consistent with previous work that showed metastable olivine persisting up to 150 km past its equilibrium transformation depth [@rubie1994; @dassler1996b].

To place the corrected and uncorrected $Z$ values in the context of the experimental data of @hosoya2005, we computed the range of $Z$ values permitted by their experimental uncertainties. Using $Z$ = $\frac{6.67}{d}\,\Gamma\,C_\mathrm{OH}^n$ with a growth-rate prefactor $\ln \Gamma$ ranging from -21.8 to -14.2 (corresponding to the reported best-fit $\ln A$ = -18.0 $\pm$ 3.8), a water-content exponent $n$ = 3.2, water contents of 1--5000 ppm OH and grain sizes $d$ = 1 and 10 mm:

| OH content | Grain Size | $\ln \Gamma$ | $Z$ |
|:-----|:---|:------|:-----------|
| 1    | 1  | -14.2 | 4.5e-3     |
| 1    | 1  | -21.8 | 2.3e-6     |
| 10   | 1  | -14.2 | 7.2e0      |
| 10   | 1  | -21.8 | 3.6e-3     |
| 100  | 1  | -14.2 | 1.1e4      |
| 100  | 1  | -21.8 | 5.7e0      |
| 1000 | 1  | -14.2 | 1.8e7      |
| 1000 | 1  | -21.8 | 9.0e3      |
| 5000 | 1  | -14.2 | **3.1e9**  |
| 5000 | 1  | -21.8 | 1.6e6      |
| 1    | 10 | -14.2 | 4.5e-4     |
| 1    | 10 | -21.8 | **2.3e-7** |
| 10   | 10 | -14.2 | 7.2e-1     |
| 10   | 10 | -21.8 | 3.6e-4     |
| 100  | 10 | -14.2 | 1.1e3      |
| 100  | 10 | -21.8 | 5.7e-1     |
| 1000 | 10 | -14.2 | 1.8e6      |
| 1000 | 10 | -21.8 | 9.0e2      |
| 5000 | 10 | -14.2 | 3.1e8      |
| 5000 | 10 | -21.8 | 1.6e5      |

Table: Range of kinetic prefactors $Z$ permitted by the experimental uncertainty of @hosoya2005. {#tbl:z-range-hosoya}

Given the sparsity of experimental data, the range of possible $Z$ values is quite large, something also pointed out by previous studies [@rubie1994]. This is precisely why our published manuscript framed the simulation results as "phenomenological" in Section 4.1 (Uncertainties and Model Limitations).

Both the corrected and uncorrected $Z$ ranges are consistent with the experimental constraints of @hosoya2005. The uncorrected range (3.0e0--7.0e7 K$^{-1}$ s$^{-1}$) lies fully within the permitted values, while the corrected range (9.5e-8--2.2e0 K$^{-1}$ s$^{-1}$) spans growth-rate prefactors from the slowest end of the experimental uncertainty ($\ln \Gamma$ = $-21.8$), combined with a nominally dry rock (< 1 ppm OH) and coarse grain size (10 mm), to the best-fit prefactor ($\ln \Gamma$ = $-18$), combined with a relatively dry rock (tens of ppm OH) and fine grain size (1 mm); equivalent conditions could also be realized with somewhat higher water contents (< 100 ppm OH) at coarser grain sizes.

The correction therefore shifts the reported kinetic conditions toward slower kinetics in relatively dry rocks, rather than the hydrated rocks implied by the uncorrected values. Nevertheless, the corrected $Z$ values remain within the experimental constraints of @hosoya2005 for dry rocks.

## List of Changes

| Location | Change |
|:---------|:-------|
| Methods (Section 2.2.3.2) | Corrected $Z$ and reframed the corresponding kinetic conditions |
| Results (Section 3.2) | Corrected $Z$ |
| Figure captions | Corrected $Z$ |
| Discussion (Section 4.1) | Contextualized the corrected $Z$ with @hosoya2005 |
| SI Table and captions | Corrected $Z$ |

Table: Complete list of changes. {#tbl:changes}

**Methods (Section 2.2.2).** The sentence describing the determination of the $Z$ range was revised to:

> The range of kinetic prefactors $Z$ used in our numerical experiments (9.5e-8--2.2e0 K$^{-1}$ s$^{-1}$) was determined by holding $H^{\ast}$ = 274 kJ mol$^{-1}$, $V^{\ast}$ = 3.0e-6 m$^3$ mol$^{-1}$, and $n$ = 3.2 constant while varying the growth-rate prefactor $\Gamma$, water content $C_\mathrm{OH}$, and grain size $d$. The kinetic conditions range from nominally dry rocks ($<$ 1 ppm OH) with coarse grain sizes (10 mm) and growth-rate prefactors at the slowest end of the experimental uncertainty of @hosoya2005 ($\ln \Gamma$ = -21.8) to relatively dry rocks (tens of ppm OH) with fine grain sizes (1 mm) at the best-fit growth-rate prefactor ($\ln \Gamma$ = -18); equivalent conditions could also be realized with somewhat higher water contents ($<$ 100 ppm OH) at coarser grain sizes. The chosen grain sizes are consistent with typical grain sizes of upper mantle xenoliths [$\sim$ 3--10 mm, @karato1984; @karato2008], and the $Z$ range is consistent with the experimental conditions of @hosoya2005 and previous numerical studies of metastable olivine wedges [@rubie1994]. Thus, our experiments approximate kinetic conditions ranging from slow kinetics in dry, coarse-grained rocks ($<$ 1 ppm OH; 10 mm; $Z$ = 9.5e-8 K$^{-1}$ s$^{-1}$) to fast kinetics in relatively dry, fine-grained rocks (tens of ppm OH; 1 mm; $Z$ = 2.2e0 K$^{-1}$ s$^{-1}$).

**Uncertainties and Model Limitations (Section 4.1).** The opening passage of the first paragraph was revised to:

> The primary quantitative uncertainty in our analysis stems from the kinetic prefactor $Z$ in the interface-controlled olivine $\Leftrightarrow$ wadsleyite growth model (Equation \ref{eq:reaction-rate}), which spans several orders of magnitude reflecting variable water contents, grain sizes (1--10 mm), and growth-rate prefactors $\Gamma$ within the experimental uncertainty of @hosoya2005. The slowest $Z$ values explored here ($\sim$ 10$^{-7}$ K$^{-1}$ s$^{-1}$) correspond to nominally dry rocks ($<$ 1 ppm OH) with coarse grain sizes (10 mm) and growth rates at the slowest end of the experimental uncertainty ($\ln \Gamma$ = -21.8), so our ultra-sluggish kinetic regime represents the slowest kinetically plausible behavior permitted by the experimental constraints, while the fastest values ($\sim$ 1 K$^{-1}$ s$^{-1}$) approach the best-fit growth rate ($\ln \Gamma$ = -18) in relatively dry (tens of ppm OH), fine-grained (1 mm) rocks.

**Footnote for the corrected version of record**:

> The kinetic prefactor $Z$ in the olivine $\Leftrightarrow$ wadsleyite reaction-rate law (Equation 13) was reported with units of K$^{-1}$ s$^{-1}$, but the published values (3.0e0--7.0e7) were mistakenly passed to ASPECT with time expressed in years (shifting $Z$ values by a factor of 3.15e7). The prefactors actually used were 9.5e-8--2.2e0 K$^{-1}$ s$^{-1}$. The $Z$ range values, figure captions, and the Supplementary Information table have been corrected accordingly. We have edited two passage in Sections 2.2.3.2 (Methods) and 4.1 (Discussion) to contextualize the correction in terms of the experimental constraints of @hosoya2005. While the conclusions of the paper are unchanged, they are the result of slower (dryer) kinetic conditions than implied in the original version.

\clearpage

# References {.unnumbered #sec:references}

::: {#refs}
:::
