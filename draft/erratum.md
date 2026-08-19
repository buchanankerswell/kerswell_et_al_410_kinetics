# Erratum: correction to Kerswell et al. (2026)

**Re:** Correction of the kinetic prefactor values

## Overview

In our manuscript (doi: 10.1029/2026JB033781), the kinetic prefactor $Z$ in the interface-controlled olivine $\Leftrightarrow$ wadsleyite reaction-rate law (Equation 13) is reported with units of K$^{-1}$ s$^{-1}$:

$$
  \dot{X} = Z\, T\, \exp\!\left(-\frac{H^{\ast} + P V^{\ast}}{R T}\right) \left(1 - \exp\!\left[-\frac{\Delta G}{R T}\right]\right) \left(1 - X\right)
$$ {#eq:reaction-rate}

Because our ASPECT simulations were set to advance time in years, the prefactor range we reported in the manuscript was shifted by a factor of 3.15e7 (seconds in a year). Rather than $Z$ = 3.0e0--7.0e7 K$^{-1}$ s$^{-1}$ as originally published, the $Z$ values actually used by ASPECT were 9.5e-8--2.2e0 K$^{-1}$ s$^{-1}$. This shift was an inadvertent mistake and not a bug in the ASPECT code; we discovered the issue when testing subsequent models.

This correction only applies to the range of $Z$ values, not to any of the published simulations or results. Since $Z$ scales the rate linearly, and every result in the paper is read directly from the simulation output with the correct unit conversions applied, all quantitative results, qualitative regime thresholds, figures, and conclusions are correct as published. Correcting the $Z$ range values does not change any other reported values.

We have corrected the $Z$ range values wherever they appear in the manuscript, Supplementary Information, and revised two passages in the Methods and Uncertainties sections.

One point of interpretation follows from the correction, and we have stated it clearly in the revised text. Per AGU's correction policy, we will also publish a footnote with the corrected article describing the change.

We apologize for the error and thank the editor for the opportunity to correct it.

Sincerely,

Buchanan Kerswell (on behalf of all coauthors)

\clearpage

## Correction Summary

Assuming that there is plenty of reactant ($X$ = 0) and $\Delta G$ of the reaction is sufficiently large such that the driving force term $\left(1 - \exp\!\left[-\frac{\Delta G}{R T}\right]\right)$ is equal to one, @eq:reaction-rate simplifies to:

$$
  \dot{X} = Z\, T\, \exp\!\left(-\frac{H^{\ast} + P V^{\ast}}{R T}\right)
$$ {#eq:reaction-rate-simple}

The consequence of using the uncorrected $Z$ range can be quantified by applying characteristic values to @eq:reaction-rate. Using the ambient mantle conditions at the top of the model domain ($P$ = 10 GPa; $T$ = 1706 K), $H^{\ast}$ = 274 kJ mol$^{-1}$, $V^{\ast}$ = 3.0e-6 m$^3$ mol$^{-1}$, and a driving-force term equal to one, we computed the reaction rate $\dot{X}$ and characteristic transformation time $\tau$ = 1/$\dot{X}$ at the five $Z$ values spanning the corrected and uncorrected ranges (Table \ref{tbl:simple-calculation}). The factor of 3.15e7 matters most in cold slabs, where the Arrhenius term $\exp\!\left(-\frac{H^{\ast} + P V^{\ast}}{R T}\right)$ is small and reaction rates are correspondingly slow; in hot plumes, reaction rates remain effectively instantaneous even after the correction.

Below is a table that summarizes the different reaction rates one could roughly expect for the corrected and uncorrected $Z$ range values.

| Published Uncorrected $Z$ (K$^{-1}$ s$^{-1}$) | Published Uncorrected Rate (Ma$^{-1}$) | Uncorrected Characteristic Transformation Time (Ma) | Corrected $Z$ (K$^{-1}$ s$^{-1}$) | Corrected Rate (Ma$^{-1}$) | Corrected Characteristic Transformation Time (Ma) |
|:------|:------|:---------------|:-------|:------|:---------------|
| 3.0e0 | 7.9e7  | 1.3e-8  | 9.5e-8 | 2.5    | 4.0e-1  |
| 2.0e2 | 5.3e9  | 1.9e-10 | 6.3e-6 | 1.7e2  | 6.0e-3  |
| 4.7e2 | 1.2e10 | 8.0e-11 | 1.5e-5 | 4.0e2  | 2.5e-3  |
| 1.8e5 | 4.8e12 | 2.1e-13 | 5.7e-3 | 1.5e5  | 6.6e-6  |
| 7.0e7 | 1.9e15 | 5.4e-16 | 2.2e0  | 5.8e7  | 1.7e-8  |

Table: Approximate reaction rates and characteristic transformation time scales (transformation time is 1/rate) if simulations used the uncorrected versus corrected $Z$ range values, evaluated at $X$ = 0, $P$ = 10 GPa, $T$ = 1706 K, $H^{\ast}$ = 274 kJ mol$^{-1}$, $V^{\ast}$ = 3.0e-6 m$^3$ mol$^{-1}$, and a driving-force term equal to one (Equation \ref{eq:reaction-rate-simple}). {#tbl:simple-calculation}

To place the corrected and uncorrected $Z$ ranges in the context of the experimental data, we computed the range of $Z$ values permitted by the experimental constraints of @hosoya2005 (Table \ref{tbl:z-range-hosoya}). Using $Z$ = (6.67/d)$\Gamma$C$_\mathrm{OH}^n$ with a growth-rate prefactor $\Gamma$ ranging from $\exp(-21.8)$ to $\exp(-14.2)$ (corresponding to the reported best-fit $\ln A$ = $-18.0 \pm 3.8$), a water-content exponent $n$ = 3.2, and grain sizes $d$ = 1 and 10 mm:

| Water content (ppm H$_2$O) | $Z$ range, $d$ = 1 mm (K$^{-1}$ s$^{-1}$) | $Z$ range, $d$ = 10 mm (K$^{-1}$ s$^{-1}$) |
|:------|:------|:------|
| 1    | 2.3e-6--4.5e-3 | 2.3e-7--4.5e-4 |
| 10   | 3.6e-3--7.2e0  | 3.6e-4--7.2e-1 |
| 100  | 5.7e0--1.1e4   | 5.7e-1--1.1e3  |
| 1000 | 9.0e3--1.8e7   | 9.0e2--1.8e6   |
| 5000 | 1.6e6--3.1e9   | 1.6e5--3.1e8   |

Table: Range of kinetic prefactors $Z$ permitted by the experimental uncertainty of @hosoya2005 (growth-rate prefactor $\Gamma$ = $\exp(-21.8)$--$\exp(-14.2)$, water-content exponent $n$ = 3.2) for water contents of 1--5000 ppm H$_2$O and grain sizes of 1 and 10 mm. {#tbl:z-range-hosoya}

Both the corrected and uncorrected $Z$ ranges are consistent with the experimental constraints of @hosoya2005 (@tbl:z-range-hosoya): the uncorrected range (3.0e0--7.0e7 K$^{-1}$ s$^{-1}$) lies fully within the permitted values, while the corrected range (9.5e-8--2.2e0 K$^{-1}$ s$^{-1}$) spans growth-rate prefactors from the slowest end of the experimental uncertainty ($\Gamma$ = $\exp(-21.8)$), combined with a nominally dry rock (< 1 ppm H$_2$O) and coarse grain size (10 mm), to the best-fit prefactor ($\Gamma$ = $\exp(-18)$), combined with a relatively dry rock (tens of ppm H$_2$O) and fine grain size (1 mm); equivalent conditions could also be realized with somewhat higher water contents (< 100 ppm H$_2$O) at coarser grain sizes. The correction therefore shifts the reported kinetic conditions toward slower kinetics in relatively dry rocks, rather than the hydrated rocks implied by the uncorrected values. Given the sparsity of experimental data, the range of possible $Z$ values is quite large, something also pointed out by previous studies [@rubie1994]. This is precisely why our published manuscript framed the simulation results as "phenomenological" in Section 4.1 (Uncertainties and Model Limitations).

The correction shifts the reported kinetic conditions toward relatively slow reaction rates in relatively dry rocks, conditions that nonetheless fall within the experimental constraints of @hosoya2005. These simulations produce results qualitatively consistent with previous work that showed metastable olivine persisting up to 150 km past its equilibrium transformation depth [@rubie1994; @dassler1996b]. Those studies also calibrated their kinetic model on available experimental data, but reached similar outcomes for different reasons: they prescribed much colder slab conditions (823--973 K versus our ~1200--1300 K at 410 km depth), whereas our sluggish kinetics arise from growth-rate prefactors and water contents within the experimental range at more moderate slab temperatures.

| Location | Change |
|:---------|:-------|
| Methods (Section 2.2.2), Equation 13 | Corrected the reported $Z$ range from 3.0e0--7.0e7 to 9.5e-8--2.2e0 K$^{-1}$ s$^{-1}$ and reframed the corresponding kinetic conditions |
| Results (Section 3.2), kinetic-regime thresholds | Corrected the $Z$ thresholds from ($Z$ $\gtrsim$ 1.8e5; 2.0e2 $\lesssim$ $Z$ $\lesssim$ 1.8e5; $Z$ $\lesssim$ 2.0e2) to ($Z$ $\gtrsim$ 5.7e-3; 6.3e-6 $\lesssim$ $Z$ $\lesssim$ 5.7e-3; $Z$ $\lesssim$ 6.3e-6) K$^{-1}$ s$^{-1}$ |
| Figures 3 and 4 captions | Corrected the reported $Z$ values (3.0e0, 4.7e2, 7.0e7) to (9.5e-8, 1.5e-5, 2.2e0) K$^{-1}$ s$^{-1}$ |
| Uncertainties (Section 4.1) | Clarified that the corrected $Z$ range reflects water contents, grain sizes, and growth-rate prefactors within the experimental constraints of @hosoya2005 |
| Supplementary Information (Table S1) | Corrected all kinetic prefactor values from K$^{-1}$ yr$^{-1}$ to K$^{-1}$ s$^{-1}$ |
| Supplementary Information (Figures S1--S8 captions) | Corrected the reported $Z$ values and units (K s$^{-1}$ to K$^{-1}$ s$^{-1}$) |

Table: Complete list of changes. {#tbl:changes}

### Methods edit

> **Methods (Section 2.2.2).** The sentence describing the determination of the $Z$ range was revised to:
>
> > The range of kinetic prefactors $Z$ used in our numerical experiments (9.5e-8--2.2e0 K$^{-1}$ s$^{-1}$) was determined by combining the growth-rate prefactor $\Gamma$, water content $C_\mathrm{OH}$, and grain size $d$ in $Z$ = (6.67/d)$\Gamma$C$_\mathrm{OH}^n$ with $H^{\ast}$ = 274 kJ mol$^{-1}$, $V^{\ast}$ = 3.0e-6 m$^3$ mol$^{-1}$, and $n$ = 3.2. These kinetic conditions range from nominally dry rocks (< 1 ppm OH) with coarse grain sizes (10 mm) and growth-rate prefactors at the slowest end of the experimental uncertainty of @hosoya2005 ($\Gamma$ = $\exp\!\left(-21.8\right)$) to relatively dry rocks (tens of ppm OH) with fine grain sizes (1 mm) at the best-fit growth-rate prefactor ($\Gamma$ = $\exp\!\left(-18\right)$); equivalent conditions could also be realized with somewhat higher water contents (< 100 ppm OH) at coarser grain sizes. The grain sizes are consistent with typical grain sizes of upper mantle xenoliths [$\sim$ 3--10 mm, @karato1984; @karato2008], and both endpoints are consistent with the experimental conditions of @hosoya2005 and previous numerical studies of metastable olivine wedges [@rubie1994]. Thus, our experiments approximate kinetic conditions ranging from slow kinetics in dry, coarse-grained rocks (< 1 ppm OH; 10 mm; $Z$ = 9.5e-8 K$^{-1}$ s$^{-1}$) to fast kinetics in relatively dry, fine-grained rocks (tens of ppm OH; 1 mm; $Z$ = 2.2e0 K$^{-1}$ s$^{-1}$).

### Uncertainties edit

> **Uncertainties and Model Limitations (Section 4.1).** The opening passage of the first paragraph was revised to:
>
> > The primary quantitative uncertainty in our analysis stems from the kinetic prefactor $Z$ in the interface-controlled olivine $\Leftrightarrow$ wadsleyite growth model (Equation \ref{eq:reaction-rate}), which spans several orders of magnitude reflecting variable water contents, grain sizes (1--10 mm), and growth-rate prefactors $\Gamma$ within the experimental uncertainty of @hosoya2005. The slowest $Z$ values explored here ($\sim$ 10$^{-7}$ K$^{-1}$ s$^{-1}$) correspond to nominally dry rocks (< 1 ppm H$_2$O) with coarse grain sizes (10 mm) and growth rates at the slowest end of the experimental uncertainty ($\Gamma$ = $\exp(-21.8)$), so our ultra-sluggish kinetic regime represents the slowest kinetically plausible behavior permitted by the experimental constraints, while the fastest values ($\sim$ 1 K$^{-1}$ s$^{-1}$) approach the best-fit growth rate ($\Gamma$ = $\exp(-18)$) in relatively dry (tens of ppm H$_2$O), fine-grained (1 mm) rocks.

### Footnote for the corrected version of record

> **Correction:** The kinetic prefactor $Z$ in the olivine $\Leftrightarrow$ wadsleyite reaction-rate law (Equation 13) was reported with units of K$^{-1}$ s$^{-1}$, but the published values (3.0e0--7.0e7) were prescribed to ASPECT with time expressed in years. The prefactors actually used were 9.5e-8--2.2e0 K$^{-1}$ s$^{-1}$. The $Z$ range values, kinetic-regime thresholds, figure captions, and the Supplementary Information table have been corrected accordingly.
>
> **Implications:** All quantitative results, figures, and conclusions are unaffected, as every result was read directly from the simulation output. The corrected $Z$ range lies at the slow end of the range permitted by the experimental constraints of @hosoya2005, corresponding to kinetic conditions of relatively dry rocks: from nominally dry (< 1 ppm H$_2$O), coarse-grained (10 mm) rocks with growth rates at the slowest end of the experimental uncertainty ($\Gamma$ = $\exp(-21.8)$) to relatively dry (tens of ppm H$_2$O), fine-grained (1 mm) rocks at the best-fit prefactor ($\Gamma$ = $\exp(-18)$). The correction thus shifts the reported kinetic conditions toward slower kinetics in relatively dry rocks. The qualitative conclusions of the paper are unchanged.


\clearpage

# References {.unnumbered #sec:references}

::: {#refs}
:::
