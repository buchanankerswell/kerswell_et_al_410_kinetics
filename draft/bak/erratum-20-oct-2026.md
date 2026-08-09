# Correction to 10.1029/2026JB033781 {.unnumbered #sec:correction}

The kinetic prefactor $Z$ in the interface-controlled olivine $\Leftrightarrow$ wadsleyite rate law (Equation 13 in 10.1029/2026JB033781) was reported with units of K$^{-1}$ s$^{-1}$:

$$
  \dot{X} = Z\, T\, \exp\!\left(-\frac{H^{\ast} + P V^{\ast}}{R T}\right) \left(1 - \exp\!\left[-\frac{\Delta G}{R T}\right]\right) \left(1 - X\right)
$$ {#eq:reaction-rate}

Because our ASPECT simulations were set to use time in years, $Z$ was inadvertently scaled by a factor of 3.15e7. This was user error in parameter passing, not a bug in ASPECT. The published range $Z$ = 3.0e0--7.0e7 K$^{-1}$ s$^{-1}$ effectively ran as $Z^{\prime}$ = 9.5e-8--2.2e0 K$^{-1}$ s$^{-1}$.

The $Z$ values have been corrected everywhere they appear in 10.1029/2026JB033781 and the Supplementary Information. We have edited two passages in Sections 2.2.3.2 and 4.1 to contextualize the correction in terms of the experimental constraints of @hosoya2005. The accompanying dataset on OSF [@kerswell2026a] is read-only and has not changed. The latest code release on Zenodo [@kerswell2026c] reproduces the corrected manuscript directly from the archived OSF dataset (by dividing $Z$ by 3.15e7). All else remains unchanged.

Every result comes directly from the simulation output, so all figures, regime thresholds, and conclusions stand as published in 10.1029/2026JB033781. However, they reflect slower, drier kinetics than the original text implied.

@hosoya2005 reports $\ln \Gamma$ = -18.0 $\pm$ 3.8 (assumed 1$\sigma$). With $Z = \frac{6.67}{d}\, \Gamma\, C_\mathrm{OH}^n$ and $n$ = 3.2, $C_\mathrm{OH}$ = 50--5000 ppm, $d$ = 1--10 mm, the permitted $Z$ range across $\pm 2\sigma$ runs from 1.4e-3 K$^{-1}$ s$^{-1}$ (10-mm grains, 50 ppm OH at $-2\sigma$) to 1.4e11 K$^{-1}$ s$^{-1}$ (1-mm grains, 5000 ppm OH at $+2\sigma$). Extending to 1 ppm OH widens the full $\pm 2\sigma$ range to 5.1e-09--1.4e11 K$^{-1}$ s$^{-1}$. The published (3.0e0--7.0e7 K$^{-1}$ s$^{-1}$) and effective (9.5e-8--2.2e0 K$^{-1}$ s$^{-1}$) ranges both fall inside the extended dry conditions (Figure \ref{fig:erratum-parameter-ranges}). Our simulations therefore effectively sampled conditions from nominally dry, coarse-grains to relatively dry (hundreds of ppm OH), fine-grains. This is consistent with @ishii2021, who showed that olivine and wadsleyite can remain kinetically dry even in hydrated slabs.

Under ambient mantle conditions at the top of the domain ($P$ = 10 GPa, $T$ = 1706 K, $H^{\ast}$ = 274 kJ mol$^{-1}$, $V^{\ast}$ = 3.0e-6 m$^3$ mol$^{-1}$), published $Z$ values give transformation times of milliseconds to days, and effective $Z^{\prime}$ values give 0.02 yr to 0.4 Ma. Our simulations therefore retain metastable olivine at warmer slab temperatures (1350--1400 K at 15 GPa, or 500 K colder than the ambient mantle) due to sluggish kinetics. This contrasts with earlier work that used different kinetics and attributed metastable olivine to much colder slab temperatures [773--923 K at 15 GPa, @rubie1994; @dassler1996b; @devaux1997; @hosoya2005].

\floatplacement{figure}{H}

![Published $Z$ range (red) and effective $Z^{\prime}$ range passed to ASPECT (blue). The thin solid box marks the experimental conditions of @hosoya2005; the dashed box marks extended dry conditions (down to 1 ppm OH). Dashed contours are calculated $Z$ values. Thick solid contours mark the permitted minimum (left panel, 5.1e-09 K$^{-1}$ s$^{-1}$) and maximum (right panel, 1.4e11 K$^{-1}$ s$^{-1}$) $Z$ values for the extended dry conditions. The blue and red ranges overlap at the 95% confidence interval ($\pm 2\sigma \ln \Gamma$). The arrow (middle panel) marks the accidental 3.15e7 scaling in $Z$.](../figs/erratum-kinetics.png){#fig:erratum-parameter-ranges width=100%}

# References {.unnumbered #sec:references}

::: {#refs}
:::
