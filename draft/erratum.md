# Correction to 10.1029/2026JB033781 {.unnumbered #sec:correction}

The kinetic prefactor $Z$ in the interface-controlled olivine $\Leftrightarrow$ wadsleyite rate law (Equation 13 in 10.1029/2026JB033781) was reported with units of K$^{-1}$ s$^{-1}$:

$$
  \dot{X} = Z\, T\, \exp\!\left(-\frac{H^{\ast} + P V^{\ast}}{R T}\right) \left(1 - \exp\!\left[-\frac{\Delta G}{R T}\right]\right) \left(1 - X\right)
$$ {#eq:reaction-rate}

Because our ASPECT simulations were set to use time in years, $Z$ was inadvertently scaled by a factor of 3.15e7. This was user error in parameter passing, not a bug in ASPECT. The published range $Z$ = 3.0e0--7.0e7 K$^{-1}$ s$^{-1}$ effectively ran as $Z^{\prime}$ = 9.5e-8--2.2e0 K$^{-1}$ s$^{-1}$.

Every result derives directly from simulation output, so all figures, regime thresholds, and conclusions stand as published. However, they reflect slower, dryer kinetics than the original text implied. Our simulations effectively sampled conditions from nominally dry, coarse-grained behavior to relatively dry (hundreds of ppm OH), fine-grained behavior.

@hosoya2005 reports $\ln \Gamma$ = -18.0 $\pm$ 3.8 (assumed 1$\sigma$). With $Z = \frac{6.67}{d}\, \Gamma\, C_\mathrm{OH}^n$ and $n$ = 3.2, $C_\mathrm{OH}$ = 1--5000 ppm, $d$ = 1--10 mm, the permitted $Z$ value range is 5.1e-09 to 1.4e11 K$^{-1}$ s$^{-1}$. The effective $Z^{\prime}$ value range falls inside this envelope near its slow end, from roughly 2.5 ppm OH in 10-mm grains at $\ln \Gamma$ = -25.6 to roughly 50 ppm OH in 1-mm grains at best fit ($\ln \Gamma$ = -18.0) or 500 ppm OH in 10-mm grains at $\ln \Gamma$ = -10.4 (Figure \ref{fig:erratum-parameter-ranges}).

![Published (red) and effective (blue) $Z$ ranges at best-fit $\Gamma$ (solid), with $\pm2\sigma$ ghosts (dashed) inside the permitted envelope (grey).](../figs/erratum-kinetics.png){#fig:erratum-parameter-ranges width=60%}

Under ambient mantle conditions at the top of the domain ($P$ = 10 GPa, $T$ = 1706 K, $H^{\ast}$ = 274 kJ mol$^{-1}$, $V^{\ast}$ = 3.0e-6 m$^3$ mol$^{-1}$), published $Z$ values give transformation times of milliseconds to days, and effective $Z^{\prime}$ values give 0.02 yr to 0.4 Ma. Our simulations therefore retain metastable olivine at warmer slab temperatures (1350–1400 K at 15 GPa, or 500 K colder than the ambient mantle) due to sluggish kinetics. This contrasts with earlier work that used different kinetics and attributed metastable olivine to much colder slab temperatures [773--923 K at 15 GPa, @rubie1994; @dassler1996b; @devaux1997; @hosoya2005].

Readers should divide $Z$ values in the manuscript by a factor of 3.15e7 to reproduce the effective $Z^{\prime}$ used in our simulations. This applies to Section 2.2.3.2, Section 3.2, the figure captions, the Figure 8 x-axis, and the Supplementary Information table. All else remains unchanged.

# References {.unnumbered #sec:references}

::: {#refs}
:::
