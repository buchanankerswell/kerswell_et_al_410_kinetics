---
#title: "Evaluating Mantle Weather"
#title: "Evaluating the Magnitudes of Nonadiabatic pressures, Temperatures, and Stresses in Earth's Mantle"
#title: "The Effects of Dynamic Stress on Phase Changes in Earth's Mantle"
#title: "Evaluating the Effects of Nonadiabatic Temperature on Phase Changes in Earth's Mantle"
title: "Evaluating the Effects of Dynamic Pressure, Temperature, and Stress on Phase Changes Earth's Mantle"
subtitle: "How Plumes and Slabs Displace Seismic Discontinuities"
author: ["Kerswell B.", "Wheeler J."]
date: "07 June 2025"
subject: "Mantle Convection"
keywords: [mantle convection, phase changes, geodynamics, numerical modeling]
lang: "en"
titlepage: true
titlepage-color: "2E7A40"
titlepage-text-color: "FFFFFF"
titlepage-rule-color: "FFFFFF"
titlepage-rule-height: 2
listings: true
tables: true
graphics: true
caption-justification: "justified"
numbersections: true
toc: true
toc-own-page: true
colorlinks: true
link-citations: true
---

# Introduction

## Research Questions

  1. Can displacements in seismic discontinuities at 410 km and 660 km be entirely explained by thermal effects? [e.g., @cottaar2016]
  2. How does the bulk composition affect thermal displacements of seismic discontinuities? [e.g., @cottaar2016]
  3. How does compressibility affect thermal displacements of seismic discontinuities? [e.g., @gassmoller2020]
  4. What thermomechanical feedbacks exist between phase changes and nonadiabatic pressures (deviatoric stresses) that might partially, or wholly, explain observed displacements in seismic discontinuities at 410 km and 660 km? [e.g., @powell2018; @wheeler2014; @tajcmanova2015; @wheeler2018; @wheeler2020]
  5. Can we describe such thermomechanical feedbacks mathematically and implement them into numerical geodynamic simulations of mantle convection? [e.g., @lu2024]
  6. Can we determine the timescales of such feedbacks?
  7. What are the implications of such feedbacks for mantle flow and observed seismic structures?

## Hypothesis Statements

  1. Nonadiabatic pressure tends to [counteract ?] [amplify ?] thermal effects at the [410 km] [660 km] seismic discontinuity
  2. The timescales of nonadiabatic pressure effects on phase transitions is [shorter ?] [longer ?] than the timescales of mantle flow
  3. Mantle flow responds to [increases] [decreases] of nonadiabatic pressure by [?]

\cleardoublepage

|Parameter name|Symbol|Unit|
|:--------------|:------|:----|
|Activation volume|$V_a$|m$^3$ mol$^{-1}$|
|Activation energy|$E_a$|J mol$^{-1}$|
|Adiabatic reference thermal expansivity|$\bar{\alpha}$|Pa$^{-1}$|
|Adiabatic reference isentropic compressibility|$\bar{\beta_S}$|Pa$^{-1}$|
|Adiabatic reference density|$\bar{\rho}$|kg m$^{-3}$|
|Adiabatic reference pressure|$\bar{P}$|K|
|Adiabatic reference specific heat|$\bar{C_p}$|J kg$^{-1}$ K$^{-1}$|
|Adiabatic reference temperature|$\bar{T}$|K|
|Adiabatic reference viscosity|$\bar{\eta}$|Pa s|
|Cohesion|$C$|Pa|
|Deviatoric stess tensor|$\sigma^{\prime}$|Pa|
|Density|$\rho$|kg m$^{-3}$|
|Effective deviatoric strain rate|$\dot{\epsilon}^{\prime}_{\text{II}}$|s$^{-1}$|
|Effective deviatoric stress|$\sigma_{\text{II}}$|Pa|
|Gas constant|$R$|J K$^{-1}$ mol$^{-1}$|
|Grain size|$d$|1 m|
|Grain size exponent|$m$|-|
|Gravitational constant|$G$|m$^3$ kg$^{-1}$ s$^{-2}$|
|Gravity vector|$g$|m s$^{-2}$|
|Initial viscosity|$\eta_0$|Pa s|
|Internal angle of friction|$\phi$|$^\circ$|
|Latent heat|$Q_L$|J kg$^{-1}$|
|Lateral viscosity variation factor|$\bar{V}$|K|
|Nonadiabatic density|$\hat{\rho}$|kg m$^{-3}$|
|Nonadiabatic pressure|$\hat{P}$|Pa|
|Nonadiabatic temperature|$\hat{T}$|K|
|Pre-exponential factor|$A$|Pa$^{-n}$ m$^{m}$ s$^{-1}$|
|Radial depth|$r$|m|
|Radial mass distribution of a spherical shell|$M(r)$|m|
|Specific heat|$C_p$|J kg$^{-1}$ K$^{-1}$|
|Deviatoric strain rate tensor|$\dot{\epsilon}^{\prime}$|s$^{-1}$|
|Stress exponent|$n$|-|
|Temperature|$T$|K|
|Thermal conductivity|$k$|W m$^{-1}$ K$^{-1}$|
|Thermal expansivity|$\alpha$|K$^{-1}$|
|Thermal viscosity exponent|$D$|-|
|Time|$t$|s|
|Total pressure|$P$|Pa|
|Velocity vector|$\vec{u}$|m s$^{-1}$|
|Viscosity|$\eta$|Pa s|
|Viscosity scaling factor|$B$|Pa$^{-n}$ m$^{m}$ s$^{-1}$|
|Yield strength|$\sigma_{\text{yield}}$|Pa|

Table: Definition of symbols. {#tbl:symbol-definitions}

\cleardoublepage

# Methods

## Simulating Compressible Mantle Flow

Mantle flow is simulated by using the finite-element geodynamic code ASPECT [v3.0.0, @aspect-doi-v3.0.0; @aspectmanual; @heister2017; @kronbichler2012; @gassmoller2018; @clevenger2021; @fraters2019; @fraters2020] to find the velocity $\vec{u}$, pressure $P$, and temperature $T$ fields that satisfy the following equations:

\begin{equation}
  \nabla P - \nabla \cdot \sigma^{\prime} = \rho g
  \label{eq:navier-stokes-no-inertia}
\end{equation}

\begin{equation}
  \frac{\partial \rho}{\partial t} + \nabla \cdot (\rho \vec{u}) = 0
  \label{eq:continuity-compressible}
\end{equation}

\begin{equation}
  \rho C_p \left( \frac{\partial T}{\partial t} + \vec{u} \cdot \nabla T \right) - \nabla \cdot \left( k \nabla T \right) = \sigma^{\prime} : \dot{\epsilon}^{\prime} + \alpha T \left( \vec{u} \cdot \nabla P \right) + Q_L
  \label{eq:energy}
\end{equation}

where $\sigma^{\prime}$ is the deviatoric stress tensor, $\rho$ is density, $g$ is gravitational acceleration, $t$ is time, $C_p$ is the specific heat capacity, $k$ is the thermal conductivity, $\alpha$ is the coefficient of thermal expansivity, and $Q_L$ is the latent heat released or absorbed during phase transformations. Equations \ref{eq:navier-stokes-no-inertia} and \ref{eq:continuity-compressible} together describe the buoyancy-driven flow of a nearly-incompressible isotropic fluid with negligible inertia and Equation \ref{eq:energy} describes the conduction, advection, and production (or consumption) of thermal energy [@schubert2001]. Note that the total pressure $P$ in Equation \ref{eq:navier-stokes-no-inertia} is equal to the mean normal stress and is positive under compression: $P = - \frac{\sigma_{xx} + \sigma_{yy} + \sigma_{zz}}{3}$ (see Appendix \ref{sec:momentum-derivation}).

The *isentropic compression approximation* [ICA, @gassmoller2020] is applied to the continuity equation (Equation \ref{eq:continuity-compressible}), which becomes:

\begin{equation}
  \nabla\cdot \vec{u} = - \beta_S \rho g \cdot \vec{u}
  \label{eq:continuity-ica}
\end{equation}

where $\beta_S$ and $\rho$ are the compressibility and density of the fluid, respectively. This formulation of the continuity equation accounts for static compression due to the combined local effects of pressure, temperature, and composition (see Appendix \ref{sec:formulations-appendix}). Our rationale for using the ICA formulation of Equation \ref{eq:continuity-compressible} (the default in ASPECT) is that it is more accurate than an incompressible fluid in cases where materials are undergoing phase transitions and/or being strongly heated or cooled [@gassmoller2020].

## Adiabatic Reference Condition {#sec:adiabatic-reference-condition}

ASPECT requires an initial guess for the material properties in Earth's mantle in order for the nonlinear solver to effectively converge on a solution [@aspectmanual; @kronbichler2012; @heister2017]. For this purpose, we computed material properties along a 1-dimensional isentropic adiabat that is characteristic of a steady-state non-convecting isochemical mantle with a potential temperature of 1573 K (Figure \ref{fig:pyr-profiles}). This adiabatic reference condition is also used as a basis for computing "nonadiabatic" quantities such as the nonadiabatic temperature $\hat{T} = T - \bar{T}$, nonadiabatic pressure $\hat{P} = P - \bar{P}$, and nonadiabatic density $\hat{\rho} = \rho - \bar{\rho}$, that quantify how much the solution deviates from a non-convecting ambient mantle.

![Adiabatic reference condition used in numerical experiments. Profiles were computed using the BurnMan software [@cottaar2014; @myhill2023] and are based on the equations of state and thermodynamic data of @stixrude2022 for a dry pyrolitic mantle.](../figs/perplex/PYR/PYR-adiabatic-profile.png){#fig:pyr-profiles width=90%}

To compute the profiles in Figure \ref{fig:pyr-profiles}, we began by evaluating entropy changes in Earth's mantle---assuming that the mantle is non-convecting (hydrostatic) and in thermodynamic equilibrium. Using the Gibbs Free Energy minimization software Perple_X [v.7.0.9, @connolly2009], we computed entropies and other material properties for a dry pyrolitic mantle [after @green1979] over a PT range of 273–3773 K and 0.001–136 GPa (Figure \ref{fig:pyr-tables}). Phase equilibria were assessed in the Na$_2$O‐CaO‐FeO‐MgO‐Al$_2$O$_3$‐SiO$_2$ (NCFMAS) chemical system with thermodynamic data and solution models of @stixrude2022. Equations of state were included for solid solution phases: olivine, plagioclase, spinel, clinopyroxene, wadsleyite, ringwoodite, perovskite, ferropericlase, high‐pressure C2/c pyroxene, orthopyroxene, akimotoite, post‐perovskite, Ca‐ferrite, garnet, and Na‐Al phase.

![PT lookup tables used in numerical experiments. Panels show density $\rho$ (a), thermal expansivity $\alpha$ (b), isentropic compressibility $\beta_S$ (c), and specific heat capacity $C_p$ (d) for a dry mantle. Material properties were computed with Perple_X using the equations of state and thermodynamic data of @stixrude2022. The white line indicates the adiabatic reference condition used in all numerical experiments. Gray regions indicate that no stable phase assemblages exists.](../figs/perplex/PYR/PYR-material-table.png){#fig:pyr-tables width=90%}

We determined the isentropic adiabat $\left( \frac{\partial P}{\partial T} \right)_S$ by applying the Newton–Raphson root finding algorithm to an array of pressures from Earth's surface to the core-mantle boundary (CMB). The method finds a corresponding temperature for each pressure such that the entropy difference between the local PT and the surface PT is zero:

\begin{equation}
  \Delta S = S(P, T) - S(P_0, T_0) = 0
\end{equation}

where $S(P, T)$ is the entropy at the local PT, and $S(P_0, T_0)$ is the entropy at the surface PT (1e5 Pa, 1573 K). Material properties were then evaluated at each PT point along the isentrope to construct the adiabatic reference profiles shown in Figure \ref{fig:pyr-profiles}.

The reference gravity $\bar{g}$ was refined by numerical integration in an iterative two-step process to ensure its self-consistency with the reference pressure $\bar{P}$ and density $\bar{\rho}$. First, depth was evaluated by integrating gravity and density over the reference pressure profile $\bar{P}$, initially assuming that gravity is a constant 9.81 m/s$^2$:

\begin{equation}
  \begin{aligned}
    \frac{d\bar{P}}{dr} &= \bar{\rho}(\bar{P}) \bar{g}(\bar{P}) \\
    r &= \int \frac{1}{\bar{\rho}(\bar{P}) \bar{g}(\bar{P})} d\bar{P}
  \end{aligned}
  \label{eq:solve-depth}
\end{equation}

where $r$ is radial depth. Once the depth profile was obtained, the reference gravity was refined by integrating the mass distribution over a spherical shell:

\begin{equation}
  \begin{aligned}
    \bar{g}(r) &= \frac{GM(r)}{r^2} \\
    M(r) &= \int 4 \pi r^2 \bar{\rho}(r) dr \\
    \bar{g}(r) &= \frac{G}{r^2} \int 4 \pi r^2 \bar{\rho}(r) dr
  \end{aligned}
  \label{eq:solve-gravity}
\end{equation}

where $M(r)$ is radial mass distribution and $G$ is the gravitational constant. Equations \ref{eq:solve-depth}–\ref{eq:solve-gravity} are repeated five times, refining both depth and gravity profiles at each iteration until gravity becomes self-consistent with the reference pressure and density along the isentrope.

## Material Models {#sec:material-models}

### Adiabatic Reference Profile (2d_shell_bm)

Experiment 2d_shell_bm allows density to vary in response to local deviations in pressure and temperature by applying a first-order Taylor expansion to $\rho(P, T)$ [@jarvis1980; @gassmoller2020]:

\begin{equation}
  \rho(P, T) \approx \bar{\rho} + \left( \frac{\partial \bar{\rho}}{\partial P} \right)_T \Delta P + \left( \frac{\partial \bar{\rho}}{\partial T} \right)_P \Delta T
  \label{eq:density-ala-expansion}
\end{equation}

and then rewriting Equation \ref{eq:density-ala-expansion} using standard thermodynamic relations $\beta = \frac{1}{\rho} \left( \frac{\partial \rho}{\partial P} \right)_T$ and $\alpha = -\frac{1}{\rho} \left( \frac{\partial \rho}{\partial T} \right)_P$ to get the expression:

\begin{equation}
  \rho(P, T) = \bar{\rho} \left( 1 + \bar{\beta_S} \hat{P} - \bar{\alpha} \hat{T} \right)
  \label{eq:density-ala}
\end{equation}

where $\bar{\rho}$, $\bar{\beta_S}$, $\bar{\alpha}$, are the adiabatic reference density, compressibility, and thermal expansivity, respectively, and $\Delta P = \hat{P} = P - \bar{P}$ and $\Delta T = \hat{T} = T - \bar{T}$ are the nonadiabatic pressure and temperature. Besides density, the other material properties used to solve Equations \ref{eq:navier-stokes-no-inertia}–\ref{eq:continuity-ica} are taken directly from the adiabatic reference condition (Figure \ref{fig:pyr-profiles}). Thus, the material model used for experiment 2d_shell_bm effectively eliminates the influence of local PT perturbations on phase transitions by assuming that deviations in material properties from the adiabatic reference condition are negligible. Under these assumptions, the governing equations become:

\begin{equation}
  \nabla P - \nabla \cdot \sigma^{\prime} = \rho(P, T) \bar{g}
  \label{eq:navier-stokes-no-inertia-bm}
\end{equation}

\begin{equation}
  \nabla\cdot \vec{u} = - \bar{\beta_S} \rho(P, T) \bar{g} \cdot \vec{u}
  \label{eq:continuity-ica-bm}
\end{equation}

\begin{equation}
  \rho(P, T) \bar{C_p} \left( \frac{\partial T}{\partial t} + \vec{u} \cdot \nabla T \right) - \nabla \cdot \left( k \nabla T \right) = \sigma^{\prime} : \dot{\epsilon}^{\prime} + \bar{\alpha} T \left( \vec{u} \cdot \nabla P \right) + Q_L
  \label{eq:energy-bm}
\end{equation}

where $\rho(P, T)$ is substituted from Equation \ref{eq:density-ala}, the coefficients $\bar{g}$, $\bar{\beta}$, $\bar{C_p}$, and $\bar{\alpha}$ are taken from the reference adiabatic condition (Figure \ref{fig:pyr-profiles}), and the thermal conductivity $k$ = 4.0 Wm$^{-1}$K$^{-1}$ is constant.

### Thermodynamic Lookup Table (2d_shell_stb and 2d_shell_vp) {#sec:lookup-table-material-model}

In contrast to experiment 2d_shell_bm, experiments 2d_shell_stb and 2d_shell_vp use the thermodynamic calculations shown in Figure \ref{fig:pyr-tables} as PT lookup tables, where material properties are found by interpolation at any PT point within the table. This method ensures that material properties accurately adjust to local changes in PT. For example, major phase transitions found at 410 km and 660 km depths at ambient mantle conditions will be displaced by cold slabs and warm plumes.

The current implementation of this material model in ASPECT [v3.0.0, @aspect-doi-v3.0.0; @aspectmanual] uses the full pressure $P$ to retrieve material properties from the PT lookup table, rather than the adiabatic reference pressure $\bar{P}$. This is not strictly consistent with the thermodynamic lookup tables, which are calculated assuming pressure is under the adiabatic reference condition $\bar{P}$ (i.e., the lookup tables do not account for the nonadiabatic component of pressure $\hat{P}$). To a first-order, the material model used for experiments 2d_shell_stb and 2d_shell_vp effectively account for local thermal and pressure perturbations on phase transitions, despite the internal inconsistency between full pressure $P$ and the adiabatic reference pressure $\bar{P}$. Under these assumptions, the governing equations become:

\begin{equation}
  \nabla P - \nabla \cdot \sigma^{\prime} = \rho(P, T) \bar{g}
  \label{eq:navier-stokes-no-inertia-stb-vp}
\end{equation}

\begin{equation}
  \nabla\cdot \vec{u} = - \beta_S(P, T) \rho(P, T) \bar{g} \cdot \vec{u}
  \label{eq:continuity-ica-stb-vp}
\end{equation}

\begin{equation}
  \rho(P, T) C_p(P, T) \left( \frac{\partial T}{\partial t} + \vec{u} \cdot \nabla T \right) - \nabla \cdot \left( k \nabla T \right) = \sigma^{\prime} : \dot{\epsilon}^{\prime} + \alpha(P, T) T \left( \vec{u} \cdot \nabla P \right) + Q_L
  \label{eq:energy-stb-vp}
\end{equation}

where only the gravity $\bar{g}$ is taken from the reference adiabatic condition, and the other coefficients $\rho(P, T)$, $\beta(P, T)$, $C_p(P, T)$, and $\alpha(P, T)$ are found by interpolation within the PT lookup table (Figure \ref{fig:pyr-tables}) at the full pressure $P$ and temperature $T$. Note that the thermal conductivity $k$ = 4.0 Wm$^{-1}$K$^{-1}$ is constant.

## Rheological Models {#sec:rheological-models}

Our numerical experiments use different rheological models to simulate deformation of mantle rocks. The models generally increase in complexity from a simple layered rheological model (2d_shell_bm), to a rheological model that has been optimized with geophysical observations [2d_shell_stb, @steinberger2006], and finally to a viscoplastic rheological model that takes into account different creep mechanisms and plastic yielding [2d_shell_vp, @glerum2018]. Our rationale for using different rheological models was to determine the effects of rheology on the development of nonadiabatic pressure $\hat{P}$ and deviatoric stress $\sigma^{\prime}$ at important phase transitions.

### Simple Temperature-Dependent Rheological Model (2d_shell_bm)

In the most simple case, the mantle is assumed to have a nominal viscosity of 1e21 Pa s, which is modified by a piecewise function that depends on radial depth $r$:

\begin{equation}
  \bar{\eta}(r) =
  \begin{cases}
    1, & r \leq 150 \text{ km} \\
    0.1, & 150 \text{ km} < r < 410 \text{ km} \\
    1, & 410 \text{ km} \leq r < 660 \text{ km} \\
    10, & r \geq 660 \text{ km}
  \end{cases}
  \label{eq:piecewise-viscosity-function}
\end{equation}

This function simulates a strong lithosphere underlain by a relatively weak asthenosphere, followed by a mantle transition zone (MTZ) with a nominal viscosity of 1e21 Pa s, which sharply transitions to a stronger lower mantle. A thermal dependency is also implemented through an exponential term:

\begin{equation}
  \eta(r, T) = \bar{\eta}(r) \eta_0 \exp \left( -D \frac{\hat{T}}{\bar{T}} \right)
  \label{eq:rheology-model-bm}
\end{equation}

where $\eta_0$ = 1e21 is the nominal background viscosity, $D$ = 10 is the thermal viscosity exponent factor, and $\hat{T} = T - \bar{T}$ is the difference between the full temperature $T$ and the adiabatic reference temperature $\bar{T}$. The resulting initial viscosity $\eta(r, T)$ used in the 2d_shell_bm models is shown in Figure \ref{fig:viscosity-profiles}a.

### Steinberger Rheological Model (2d_shell_stb)

@steinberger2006 described an approach to constrain a layered mantle rheology through an optimization procedure that accounts for Earth's geoid, glacial rebound data, and a radial heat flux model based on mineral physics data. Their optimized rheological model is implemented by modifying a reference viscosity $\bar{\eta}(r)$ with a coefficient $\bar{V}(r)$ that approximates variations in viscosity due to temperature and different creep mechanisms:

\begin{equation}
  \eta(r, T) = \bar{\eta}(r) \exp \left( -\frac{\bar{V}(r) \hat{T}}{T \bar{T}} \right)
  \label{eq:rheology-model-stb}
\end{equation}

where $\hat{T}$ and $\bar{T}$ are the same as in Equation \ref{eq:rheology-model-bm}. Mantle viscosity in this case is also arbitrarily constrained to stay within 1e20 and 5e24 Pa s. The resulting initial viscosity $\eta(r, T)$ is shown in Figure \ref{fig:viscosity-profiles}a.

### Viscoplastic Rheological Model (2d_shell_vp) {#sec:vp-model}

The effective creep viscosity for both Newtonian (linear) diffusion creep and non-Newtonian (non-linear) dislocation creep are defined by the same equation [after @karato2008; @karato1993]:

\begin{equation}
  \eta^{visc} = \frac{1}{2} B A^{-1/n} d^{m/n} \dot{\epsilon}_{\text{II}}^{(1-n)/n} \exp \left( \frac{E_a + V_a P}{nRT} \right)
  \label{eq:viscous-creep}
\end{equation}

where $\eta^{visc}$ is the effective creep viscosity, $d$ is the grain size, $m$ is the grain size exponent, $n$ is the stress exponent, $A$ is the pre-exponential factor, $\dot{\epsilon}^{\prime}_{\text{II}} = \sqrt{\frac{1}{2} \Sigma \dot{\epsilon}^{\prime 2}_{ij}}$ is the second invariant of the deviatoric strain rate tensor (see Appendix \ref{sec:momentum-derivation}), $E_a$ is the activation energy, $P$ is the full pressure, $V_a$ is the activation volume, $R$ is the gas constant, $T$ is the full temperature, and $B$ is a scaling factor that is used to arbitrarily tune the effective viscosity [@glerum2018].

![Evolution of mantle viscosity in numerical experiments. The initial viscosity (a) evolves to its end state (b) at $t$ = 500 Ma. In general, model 2d_shell_bm has the weakest upper and lower mantle, 2d_shell_stb has a relatively strong upper mantle, but moderately weak MTZ and lower mantle, while 2d_shell_vp has the highest overall strength, especially in the lower mantle. The evolution of viscosities with time depends on the rheological models (c.f., Equations \ref{eq:rheology-model-bm}, \ref{eq:rheology-model-stb}, \ref{eq:rheology-model-vp}), which control the vigor of mantle convection.](../figs/misc/combined-viscosity-profiles.png){#fig:viscosity-profiles width=90%}

Both dislocation creep and diffusion creep mechanisms are active to some degree in Earth's mantle, but laboratory experiments suggest that dislocation creep should be the predominant deformation mechanism in Earth's upper mantle while diffusion creep is more active in the lower mantle [@karato2008; @ranalli1995; @schubert2001; @vanhunen2005]. Thus, our numerical experiments assume that diffusion creep is mostly inactive above 660 km ($n$ = 3, $m$ = 1), and dislocation creep is completely inactive below 660 km ($n$ = 1, $m$ = 3), which is common practice in numerical geodynamic models [e.g.,  @billen2007; @thieulot2011; @gerya2019; @glerum2018; @vanhunen2005; @vanderwiel2024]. To achieve this, we set the coefficients $A$, $m$, $n$, $E_a$, and $V_a$ in Equation \ref{eq:viscous-creep} independently in the upper mantle, MTZ, and lower mantle (Table \ref{tbl:rheological-parameters-vanderwiel20204}). The grain size $d$ is assumed to be 1 m for the entire mantle and the scaling factor $B$ is set to 1.

|Depth       |Creep Mechanism|     $A$|$m$|$n$| $E_a$|   $V_a$|$\phi$|$C$|
|:-----------|:--------------|-------:|--:|--:|-----:|-------:|-----:|--:|
|Upper mantle|Diffusion      |6.00e-17|  1|  1|1.50e5| 6.34e-7|1.4337|4e6|
|Upper mantle|Dislocation    |6.51e-16|  1|  3|5.00e5| 1.30e-5|1.4337|4e6|
|MTZ         |Diffusion      |9.00e-17|  1|  1|1.55e5|14.34e-7|1.4337|4e6|
|MTZ         |Dislocation    |8.51e-16|  1|  3|5.00e5| 1.30e-5|1.4337|4e6|
|Lower mantle|Diffusion      |1.00e-18|  3|  1|1.50e5|10.34e-7|1.4337|4e6|
|Lower mantle|Dislocation    |6.51e-16|  1|  1|5.30e5| 1.30e-5|1.4337|4e6|

Table: Material parameters after @vanderwiel2024 for a viscous flow law (Equation \ref{eq:viscous-creep}). Grain size $d$ is set to 1 m and the scale factor $B$ to 1. Units are $A$: $\text{Pa}^{-n} \text{m}^{m} \text{s}^{-1}$, $m$: none, $n$: none, $V_a$: $\text{m}^3 \text{mol}^{-1}$, $E_a$ $\text{J} \text{mol}^{-1}$, $\phi$: $^\circ$, $C$: Pa. {#tbl:rheological-parameters-vanderwiel20204}

A composite viscous creep rheology is calculated by averaging the effective diffusion and dislocation viscosities harmonically [@gerya2019]:

\begin{equation}
  \eta^{\text{creep}} = \left( \frac{1}{\eta^{\text{diff}}} + \frac{1}{\eta^{\text{disl}}} \right)^{-1}
  \label{eq:combined-creep}
\end{equation}

Finally, plastic yielding is implemented by applying the Druker-Prager yield criterion [@davis2005]:

\begin{equation}
  \sigma_{\text{yield}} = C cos(\phi) + sin(\phi) P
  \label{eq:druker-prager}
\end{equation}

where $\sigma_{\text{yield}}$ is the stress at which the material begins to yield in a plastic manner, $C$ is the material cohesion, $\phi$ is the internal friction angle, and $P$ is pressure. An effective yield viscosity is defined as [@kachanov2004]:

\begin{equation}
  \eta^{\text{pl}} = \frac{\sigma_{\text{yield}}}{2 \dot{\epsilon}_{\text{II}}}
  \label{eq:plastic-viscosity}
\end{equation}

which is used to limit the effective creep viscosity if plastic yielding is occurring [@karato2008]:

\begin{equation}
  \eta^{\text{vp}} = \text{min}(\eta^{\text{creep}}\text{, } \eta^{\text{pl}})
  \label{eq:rheology-model-vp}
\end{equation}

## Numerical Setup

Our numerical experiments simulate mantle convection in a 2d annulus with global parameters defined in Table \ref{tbl:global-parameters}. The simulations begin with a cooled outer boundary layer that has an initial thermal age of 50 Ma. A constant temperature of 273 K is maintained at the outer boundary with a imposed steady-state velocity representing present-day plate motions around Earth's equator. The inner boundary is free-slip with a constant temperature of 3773 K. The finite element mesh is initially refined 4 times globally. Adaptive mesh refinement is then applied during each timestep to 1) maintain a minimum refinement level of 9 between 350–720 km radial depth, and 2) adjust refinement levels according to local nonadiabatic temperature $\hat{T}$. The simulations are ran from this initial state for 500 Ma.

|Parameter|Value|Units|
|:--|:--|:--|
|Inner radius|3480|km|
|Outer radius|6370|km|
|Stokes solver tolerance|1e-6|-|
|CFL Number|0.7|-|
|Material properties|Figures \ref{fig:pyr-profiles}-\ref{fig:pyr-tables}|-|
|Viscosity|Figure \ref{fig:viscosity-profiles}|-|
|Mantle potential temperature|1573|K|
|Outer boundary temperature|273|K|
|Inner boundary temperature|3773|K|
|Age of outer boundary layer|50|Ma|
|Min. mesh refinement levels|4|-|
|Max. mesh refinement levels|9|-|
|Formulation|Isentropic compression|-|
|ASPECT version|3.0.0|-|
|Experimental parameters|Appendix \ref{sec:prm-configs}|-|
|Cluster|BARKLA2 (Univ. of Liverpool)|-|
|Nodes/CPUs|3/120|-|
|Cluster configuration|Appendix \ref{sec:barkla2-config}|-|

Table: Global parameters for mantle convection models. {#tbl:global-parameters}

Note that although the mesh geometry and refinement functions are formulated with respect to radial coordinates, the governing equations are solved in a xyz Cartesian coordinate reference frame. For example, the stress and strain rate tensors are defined using indices ij that correspond to the Cartesian directions xyz (see Appendix \ref{sec:momentum-derivation}). The origin of the Cartesian reference frame is at the center of the annulus and the gravity vector always points towards the origin.

\cleardoublepage

# Results

## Experiment 2d_shell_bm

Mantle convection in experiment 2d_shell_bm is driven by cold slabs that sink from the outer boundary where the surface velocity is converging. Two-sided subduction and whole mantle convection develop relatively rapidly (< 100 Ma) due to the relatively weak mantle rheology (Figures \ref{fig:viscosity-profiles} and \ref{fig:bm-mantle-flow}, see also Equation \ref{eq:rheology-model-bm}), especially in the lower mantle where 2d_shell_bm is about 10x weaker than experiments 2d_shell_stb and 2d_shell_vp (Figure \ref{fig:viscosity-profiles}). The cold dense slabs penetrate the MTZ and sink vertically as coherent columns until they deform in the deep lower mantle and begin to pile up against the CMB. Meanwhile, warm buoyant plumes generated at the CMB rise quickly towards the surface where they disperse laterally under the base of the thin lithosphere (Figure \ref{fig:bm-mantle-flow}).

Effective deviatoric stress $\sigma_{\text{II}}$ and nonadiabatic pressure $\hat{P}$ are most significant at: 1) the inner boundary due to piling slab material near the CMB, and 2) the outer boundary due to the imposed surface velocity condition (Figure \ref{fig:bm-binned-depth-zoomed-view}). Besides these boundary effects, the highest deviatoric stress (100–200 MPa) and nonadiabatic pressure (50–100 MPa) in experiment 2d_shell_bm occur within slabs in the lower mantle as they bend (Figure \ref{fig:bm-binned-depth-zoomed-view}). Elevated deviatoric stress also occurs within slabs at the base of the lithosphere (125 km) and at the 410 km and 660 km phase transitions. While some of these regions of high deviatoric stress within slabs correlate with high nonadiabatic pressures (e.g., bending slabs in the lower mantle), the correlation between elevated deviatoric stress and dynamic pressure does not hold everywhere (e.g., the 660 km phase transition).

In contrast to cold slabs, warm plumes do not generally generate deviatoric stress or nonadiabatic pressure greater than $\pm$ 30–50 MPa. There are two exceptions, however, where dynamic effects are observed in plumes: 1) a thin layer of elevated deviatoric stress at 125 km (up to 180 MPa) where plumes are dispersing laterally along the base of the lithosphere, and 2) elevated nonadiabatic pressure in plume heads and tails (up to $\pm$ 100 MPa, Figure \ref{fig:bm-binned-depth-zoomed-view}). The nonadiabatic effects observed in plumes are relatively short-lived due to their relatively high temperatures and high velocities.

\cleardoublepage

![](../figs/simulation/full_mantle_with_piecewise_viscosity_profile/tiles/temperature-nonadiabatic-density-nonadiabatic-viscosity-0020-full-view-1.png)

![](../figs/simulation/full_mantle_with_piecewise_viscosity_profile/tiles/temperature-nonadiabatic-density-nonadiabatic-viscosity-0033-full-view-1.png)

![Experiment 2d_shell_bm after 100 Ma (top), 165 Ma (middle), and 235 Ma (bottom) showing the global distribution of nonadiabatic temperature $\hat{T}$ (left), nonadiabatic density $\hat{\rho}$ (middle), and log viscosity $\eta$ (right). Models were calculated with ASPECT using the ICA formulation (Equation \ref{eq:continuity-ica}) and a simple temperature dependent rheological model (Equation \ref{eq:rheology-model-bm}). Initial and final viscosities are shown in Figure \ref{fig:viscosity-profiles}. Colorbar limits are the same for all numerical experiments. Gray regions indicate where the mantle has not deviated from the adiabatic reference condition (Figure \ref{fig:pyr-profiles}). Values outside of the colorbar limits are filled white.](../figs/simulation/full_mantle_with_piecewise_viscosity_profile/tiles/temperature-nonadiabatic-density-nonadiabatic-viscosity-0047-full-view-1.png){#fig:bm-mantle-flow}

![Experiment 2d_shell_bm after 100 Ma of evolution. A representative view of the model domain (top) showing nonadiabatic density $\hat{\rho}$, nonadiabatic pressure $\hat{P}$, and effective deviatoric stress $\sigma_{\text{II}}$. Colorbar limits are the same for all numerical experiments. Gray regions indicate where the mantle has not deviated from the adiabatic reference condition (Figure \ref{fig:pyr-profiles}). Values outside of the colorbar limits are filled white. Bar charts show global averages, binned approximately every 25 km, for slabs (middle) and plumes (bottom). Thin black lines mark the base of the lithosphere at 125 km and major phase transformations at 410 and 660 km.](../figs/simulation/full_mantle_with_piecewise_viscosity_profile/tiles_comps/density-nonadiabatic-pressure-nonadiabatic-sigma-ii-0020-zoomed-view.png){#fig:bm-binned-depth-zoomed-view}

![Experiment 2d_shell_bm after 100 Ma of evolution. A representative view of the model domain (top) showing the normal (xx and yy) and shear (xy = yx) components of the deviatoric stress tensor $\sigma^{\prime}$ (top). Colorbar limits are the same for all numerical experiments. Gray regions indicate where the mantle has not deviated from the adiabatic reference condition (Figure \ref{fig:pyr-profiles}). Values outside of the colorbar limits are filled white. Bar charts show global averages, binned approximately every 25 km, for slabs (middle) and plumes (bottom). Thin black lines mark the base of the lithosphere at 125 km and major phase transformations at 410 and 660 km.](../figs/simulation/full_mantle_with_piecewise_viscosity_profile/tiles_comps/sigma-deviatoric-xx-sigma-deviatoric-yy-sigma-deviatoric-xy-0020-zoomed-view.png){#fig:bm-binned-depth-zoomed-view2}

\cleardoublepage

## Experiment 2d_shell_stb

Like experiment 2d_shell_bm, mantle convection in experiment 2d_shell_stb is driven by cold downwellings that sink from convergent plate boundaries at the surface. However, the slightly stronger rheological model implemented in 2d_shell_stb (Figures \ref{fig:viscosity-profiles} and \ref{fig:stb-mantle-flow}, see also Equation \ref{eq:rheology-model-stb}) results in slower mantle convection rates and stronger slabs. This is evident by slabs that sink subvertically and bend broadly in the lower mantle while remaining rigid and intact (Figure \ref{fig:stb-mantle-flow}), rather than sink vertically and pile atop the CMB in tight folds (Figure \ref{fig:bm-mantle-flow}). Moreover, it takes about twice as long to initiate whole mantle convection in experiment 2d_shell_stb (> 200 Ma) compared to experiment 2d_shell_bm (> 100 Ma).

As observed in experiment 2d_shell_bm, effective deviatoric stress $\sigma_{\text{II}}$ and nonadiabatic pressure $\hat{P}$ are about an order of magnitude higher in cold slabs than in warm plumes. However, the overall magnitudes of these dynamic effects in experiment 2d_shell_stb are approximately 2–3x higher than in experiment 2d_shell_bm (c.f. Figures \ref{fig:bm-binned-depth-zoomed-view} and \ref{fig:bm-binned-depth-zoomed-view2} with Figures \ref{fig:stb-binned-depth-zoomed-view} and \ref{fig:stb-binned-depth-zoomed-view2}). Additionally, abrupt changes in material properties at major phase transitions are observed in experiment 2d_shell_stb, but not in experiment 2d_shell_bm. For example, nonadiabatic densities of - 0.02 $\leq$ $\hat{\rho}$ $\leq$ 0.15 g/cm$^3$ are observed at the boundaries of the MTZ and correlate with with high nonadiabatic pressures of -500 $\leq$ $\hat{P}$ $\leq$ 0 MPa and effective deviatoric stresses of $\sigma^{\prime}$ $\leq$ 400 MPa (Figure \ref{fig:stb-binned-depth-zoomed-view}). Stronger dynamic effects occur in experiment 2d_shell_stb for two reasons: 1) the stronger rheological model can support higher effective deviatoric stresses, and 2) the material model (thermodynamic lookup table) allows nonadiabatic temperature $\hat{T}$ and pressure $\hat{P}$ in slabs and plumes to drive displacements in phase transition depths (see Section \ref{sec:lookup-table-material-model}).

\cleardoublepage

![](../figs/simulation/full_mantle_with_steinberger_viscosity_profile/tiles/temperature-nonadiabatic-density-nonadiabatic-viscosity-0047-full-view-1.png)

![](../figs/simulation/full_mantle_with_steinberger_viscosity_profile/tiles/temperature-nonadiabatic-density-nonadiabatic-viscosity-0073-full-view-1.png)

![Experiment 2d_shell_stb after 235 Ma (top), 365 Ma (middle), and 500 Ma (bottom) showing the global distribution of nonadiabatic temperature $\hat{T}$ (left), nonadiabatic density $\hat{\rho}$ (middle), and log viscosity $\eta$ (right). Models were calculated with ASPECT using the ICA formulation (Equation \ref{eq:continuity-ica}) and a  rheological model after @steinberger2006 (Equation \ref{eq:rheology-model-stb}). Initial and final viscosities are shown in Figure \ref{fig:viscosity-profiles}. Colorbar limits are the same for all numerical experiments. Gray regions indicate where the mantle has not deviated from the adiabatic reference condition (Figure \ref{fig:pyr-profiles}). Values outside of the colorbar limits are filled white.](../figs/simulation/full_mantle_with_steinberger_viscosity_profile/tiles/temperature-nonadiabatic-density-nonadiabatic-viscosity-0100-full-view-1.png){#fig:stb-mantle-flow}

![Experiment 2d_shell_stb after 235 Ma of evolution. A representative view of the model domain (top) showing nonadiabatic density $\hat{\rho}$, nonadiabatic pressure $\hat{P}$, and effective deviatoric stress $\sigma_{\text{II}}$. Colorbar limits are the same for all numerical experiments. Gray regions indicate where the mantle has not deviated from the adiabatic reference condition (Figure \ref{fig:pyr-profiles}). Values outside of the colorbar limits are filled white. Bar charts show global averages, binned approximately every 25 km, for slabs (middle) and plumes (bottom). Thin black lines mark the base of the lithosphere at 125 km and major phase transformations at 410 and 660 km.](../figs/simulation/full_mantle_with_steinberger_viscosity_profile/tiles_comps/density-nonadiabatic-pressure-nonadiabatic-sigma-ii-0047-zoomed-view.png){#fig:stb-binned-depth-zoomed-view}

![Experiment 2d_shell_stb after 235 Ma of evolution. A representative view of the model domain (top) showing the normal (xx and yy) and shear (xy = yx) components of the deviatoric stress tensor $\sigma^{\prime}$ (top). Colorbar limits are the same for all numerical experiments. Gray regions indicate where the mantle has not deviated from the adiabatic reference condition (Figure \ref{fig:pyr-profiles}). Values outside of the colorbar limits are filled white. Bar charts show global averages, binned approximately every 25 km, for slabs (middle) and plumes (bottom). Thin black lines mark the base of the lithosphere at 125 km and major phase transformations at 410 and 660 km.](../figs/simulation/full_mantle_with_steinberger_viscosity_profile/tiles_comps/sigma-deviatoric-xx-sigma-deviatoric-yy-sigma-deviatoric-xy-0047-zoomed-view.png){#fig:stb-binned-depth-zoomed-view2}

\cleardoublepage

## Experiment 2d_shell_vp

Mantle convection in experiment 2d_shell_vp is the most sluggish overall, which is a result of having the strongest rheological model, especially in the lower mantle (see Equation \ref{eq:rheology-model-vp} and Figure \ref{fig:viscosity-profiles}). For example, warm plumes only begin to rise from the CMB after approximately 235 Ma of evolution (Figure \ref{fig:vp-mantle-flow}), whereas whole mantle convection has already completed in experiment 2d_shell_bm by 100 Ma, and in 2d_shell_stb by 235 Ma (c.f. Figures \ref{fig:bm-mantle-flow}, \ref{fig:stb-mantle-flow}, and \ref{fig:vp-mantle-flow}). Despite having a strong rheological model, however, sinking slabs in experiment 2d_shell_vp become completely dismembered or tightly folded due to high deviatoric stresses (up to 200–500 MPa) that exceed their plastic yield strength (Figure \ref{fig:vp-mantle-flow}).

Unlike 2d_shell_bm and 2d_shell_stb, the magnitudes of nonadiabatic pressures in experiment 2d_shell_vp can be roughly equivalent in cold slabs and warm plumes despite large differences in their deviatoric stress. For example, the average nonadiabatic pressure in both slabs and plumes in the lower mantle is 200–250 MPa, but the average deviatoric stress are approximately 200 MPa in slabs vs. 50 MPa in plumes (Figure \ref{fig:vp-binned-depth-zoomed-view}). Moreover, deviatoric stress in slabs can be extremely high, exceeding 500 MPa where the slabs are bending in the lower mantle (Figure \ref{fig:vp-binned-depth-zoomed-view}). As observed in experiment 2d_shell_stb, these regions of high deviatoric stress generally correlate with high nonadiabatic pressures and nonadiabatic densities (e.g., at phase transitions and deforming slabs), but not everywhere (e.g., near the CMB, Figure \ref{fig:vp-binned-depth-zoomed-view}).

\cleardoublepage

![](../figs/simulation/full_mantle_with_viscoplastic_viscosity_profile/tiles/temperature-nonadiabatic-density-nonadiabatic-viscosity-0020-full-view-1.png)

![](../figs/simulation/full_mantle_with_viscoplastic_viscosity_profile/tiles/temperature-nonadiabatic-density-nonadiabatic-viscosity-0033-full-view-1.png)

![Experiment 2d_shell_vp after 100 Ma (top), 165 Ma (middle), and 235 Ma (bottom) showing the global distribution of nonadiabatic temperature $\hat{T}$ (left), nonadiabatic density $\hat{\rho}$ (middle), and log viscosity $\eta$ (right). Models were calculated with ASPECT using the ICA formulation (Equation \ref{eq:continuity-ica}) and a viscoplastic rheological model (Equation \ref{eq:rheology-model-vp}). Initial and final viscosities are shown in Figure \ref{fig:viscosity-profiles}. Colorbar limits are the same for all numerical experiments. Gray regions indicate where the mantle has not deviated from the adiabatic reference condition (Figure \ref{fig:pyr-profiles}). Values outside of the colorbar limits are filled white.](../figs/simulation/full_mantle_with_viscoplastic_viscosity_profile/tiles/temperature-nonadiabatic-density-nonadiabatic-viscosity-0047-full-view-1.png){#fig:vp-mantle-flow}

![Experiment 2d_shell_vp after 100 Ma of evolution. A representative view of the model domain (top) showing nonadiabatic density $\hat{\rho}$, nonadiabatic pressure $\hat{P}$, and effective deviatoric stress $\sigma_{\text{II}}$. Colorbar limits are the same for all numerical experiments. Gray regions indicate where the mantle has not deviated from the adiabatic reference condition (Figure \ref{fig:pyr-profiles}). Values outside of the colorbar limits are filled white. Bar charts show global averages, binned approximately every 25 km, for slabs (middle) and plumes (bottom). Thin black lines mark the base of the lithosphere at 125 km and major phase transformations at 410 and 660 km.](../figs/simulation/full_mantle_with_viscoplastic_viscosity_profile/tiles_comps/density-nonadiabatic-pressure-nonadiabatic-sigma-ii-0047-zoomed-view.png){#fig:vp-binned-depth-zoomed-view}

![Experiment 2d_shell_vp after 100 Ma of evolution. A representative view of the model domain (top) showing the normal (xx and yy) and shear (xy = yx) components of the deviatoric stress tensor $\sigma^{\prime}$ (top). Colorbar limits are the same for all numerical experiments. Gray regions indicate where the mantle has not deviated from the adiabatic reference condition (Figure \ref{fig:pyr-profiles}). Values outside of the colorbar limits are filled white. Bar charts show global averages, binned approximately every 25 km, for slabs (middle) and plumes (bottom). Thin black lines mark the base of the lithosphere at 125 km and major phase transformations at 410 and 660 km.](../figs/simulation/full_mantle_with_viscoplastic_viscosity_profile/tiles_comps/sigma-deviatoric-xx-sigma-deviatoric-yy-sigma-deviatoric-xy-0047-zoomed-view.png){#fig:vp-binned-depth-zoomed-view2}

\cleardoublepage

# Discussion

  1. What is the relationship between $\sigma_\text{II}$ and $\hat{P}$?
  2. How big is the thermal displacement of the 410 and 660?

# Acknowledgements

This work was funded by the UKRI NERC Large Grant no. NE/V018477/1 awarded to John Wheeler at the University of Liverpool. All computations were undertaken on Barkla2, part of the High Performance Computing facilities at the University of Liverpool, who graciously provided expert support. We thank the Computational Infrastructure for Geodynamics ([geodynamics.org](geodynamics.org)) which is funded by the National Science Foundation under award EAR-0949446 and EAR-1550901 for supporting the development of ASPECT.

# Data Availability

The code modifications, parameter, data, and log files used for the models in the study are available at ...

ASPECT version 3.0.0, [@aspect-doi-v3.0.0; @aspectmanual; @heister2017; @kronbichler2012; @gassmoller2018; @clevenger2021; @fraters2019; @fraters2020] used in these computations is freely available under the GPL v2.0 or later license through its software landing page [https://geodynamics.org/resources/aspect](https://geodynamics.org/resources/aspect) or [https://aspect.geodynamics.org](https://aspect.geodynamics.org) and is being actively developed on GitHub and can be accessed via [https://github.com/geodynamics/aspect](https://github.com/geodynamics/aspect).

# References
::: {#refs}
:::

\cleardoublepage

# Appendix {.unnumbered}

## Governing Equations {.unnumbered}

### The Momentum Equation {.unnumbered #sec:momentum-derivation}

The momentum equation in the following form is referred to as the Navier-Stokes equation and describes the flow of a compressible viscous fluid primarily by buoyancy forces:

\begin{equation}
  \rho \left( \frac{\partial \vec{u}}{\partial t} \right) = \nabla \cdot \sigma^{\prime} - \nabla P + \rho g
  \label{eq:navier-stokes-compressible}
\end{equation}

where $\rho$ is density, $\vec{u}$ is velocity, $P$ is pressure, $\sigma^{\prime}$ is the deviatoric stress tensor, and $g$ is gravitational acceleration. In terms of classical mechanics, the left hand side is analogous to mass times acceleration $ma$, and the right-hand side are the forces $F$ that are acting on the fluid. Hence, the equation describes a balance between force and momentum $ma = F$ [@gerya2019].

The forces in Equation \ref{eq:navier-stokes-compressible} include the pressure gradient $\nabla P$ which acts to drive the fluid from high pressure to low pressure, the viscous forces $\nabla \cdot \sigma^{\prime}$ that dissipate energy by resisting flow, and the buoyancy forces $\rho g$ that drive convection (denser fluids sink while lighter fluids rise). Because the flow of Earth's mantle occurs at such slow rates, however, the inertial term $\left( \frac{\partial \vec{u}}{\partial t} \right)$ on the left hand side of Equation \ref{eq:navier-stokes-compressible} can be ignored and the momentum equation simplifies to:

\begin{equation}
  \nabla P - \nabla \cdot \sigma^{\prime} = \rho g
  \label{eq:navier-stokes-no-inertia-appendix}
\end{equation}

Equation \ref{eq:navier-stokes-no-inertia-appendix} describes a balance between the buoyancy force and the pressure gradient minus the energy dissipation due to deformation. Since the deviatoric stress tensor $\sigma^{\prime}$ can be described in terms of velocity (see Equation \ref{eq:stress-deviatoric-component}), the primary unknowns in Equation \ref{eq:navier-stokes-compressible} are pressure and velocity.

The complete stress tensor can be written as:

\begin{equation}
  \sigma_{ij} = \sigma^{\prime}_{ij} - P \delta_{ij}
  \label{eq:stress-complete}
\end{equation}

where $\sigma_{ij}$ is the total stress, $\sigma^{\prime}_{ij}$ is the deviatoric (non-hydrostatic) component of stress, $P = - \frac{\sigma_{xx} + \sigma_{yy} + \sigma_{zz}}{3}$ is the hydrostatic component of stress, and $\delta_{ij}$ is the Kroneker delta:

\begin{equation}
  \delta_{ij} =
  \begin{cases}
    1, & \text{if } i = j \\
    0, & \text{if } i \neq j
  \end{cases}
\end{equation}

The hydrostatic component of stress acts equally in all directions and therefore affects the fluid's volume (density) but does not change its shape or cause it to flow. Note that the negative sign in Equation \ref{eq:stress-complete} implies that pressure is positive under compression (negative normal stress). This is a convention used in geodynamics that differs from material sciences and other fields.

The deviatoric part of the stress tensor is responsible for deformation and flow of the fluid and is equal to the total stress without the hydrostatic stress component, $\sigma^{\prime}_{ij} = \sigma_{ij} + P \delta_{ij}$, or in full matrix form:

\begin{equation}
  \begin{pmatrix}
  \sigma^{\prime}_{xx} & \sigma^{\prime}_{xy} & \sigma^{\prime}_{xz} \\
  \sigma^{\prime}_{yx} & \sigma^{\prime}_{yy} & \sigma^{\prime}_{yz} \\
  \sigma^{\prime}_{zx} & \sigma^{\prime}_{zy} & \sigma^{\prime}_{zz}
  \end{pmatrix} =
  \begin{pmatrix}
  \sigma_{xx} & \sigma_{xy} & \sigma_{xz} \\
  \sigma_{yx} & \sigma_{yy} & \sigma_{yz} \\
  \sigma_{zx} & \sigma_{zy} & \sigma_{zz}
  \end{pmatrix} +
  \begin{pmatrix}
  -\frac{\sigma_{xx} + \sigma_{yy} + \sigma_{zz}}{3} & 0 & 0 \\
  0 & -\frac{\sigma_{xx} + \sigma_{yy} + \sigma_{zz}}{3} & 0 \\
  0 & 0 & -\frac{\sigma_{xx} + \sigma_{yy} + \sigma_{zz}}{3}
  \end{pmatrix}
  \label{eq:stress-deviatoric}
\end{equation}

In practice, the deviatoric stress tensor $\sigma^{\prime}$ is computed by applying a constitutive relationship between stress and strain to express $\sigma^{\prime}$ in terms of of velocity. In the present case, we apply a generalized linear model that combines shear deformation without rotation and volumetric deformation (dilation):

\begin{equation}
  \sigma^{\prime} = \eta \left( \nabla \vec{u} + \left( \nabla \vec{u} \right)^\intercal \right) - \left( \frac{2}{3} \eta - \zeta \right) \left( \nabla \cdot \vec{u} \right) I
\end{equation}

where $\vec{u}$ is velocity, $\eta$ is the shear viscosity, and $\zeta$ is the bulk viscosity. For nearly-incompressible fluids (very small $\zeta$), the expression for deviatoric stress simplifies to:

\begin{equation}
  \sigma^{\prime} = 2 \eta \dot{\epsilon}^{\prime}
\end{equation}

where $\dot{\epsilon}^{\prime} = \frac{1}{2} \left( \nabla \vec{u} + \left( \nabla \vec{u} \right)^\intercal \right) - \frac{1}{3} \left( \nabla \cdot \vec{u} \right) I$ is the deviatoric strain rate tensor. In full component form the deviatoric stress tensor is:

\begin{equation}
  \sigma^{\prime}_{ij} = \eta \left( \frac{\partial u_i}{\partial x_j} + \frac{\partial u_j}{\partial x_i} \right) - \frac{2}{3} \eta \left( \frac{\partial u_x}{\partial x} + \frac{\partial u_y}{\partial y} + \frac{\partial u_z}{\partial z} \right) \delta_{ij}
  \label{eq:stress-deviatoric-component}
\end{equation}

Note that the deviatoric stress and strain rate tensors are symmetric such that $\sigma^{\prime}_{ij} = \sigma^{\prime}_{ji}$ and $\dot{\epsilon}^{\prime}_{ij} = \dot{\epsilon}^{\prime}_{ji}$, which implies that there is zero rigid-body rotation in the fluid flow. Because of this symmetry, the full matrix form the deviatoric stress tensor can be written as:

\begin{equation}
  \sigma^{\prime} = \begin{pmatrix}
  2 \eta \frac{\partial u_x}{\partial x} - \frac{2}{3} \eta \left( \frac{\partial u_x}{\partial x} + \frac{\partial u_y}{\partial y} + \frac{\partial u_z}{\partial z} \right) & \eta \left( \frac{\partial u_x}{\partial y} + \frac{\partial u_y}{\partial x} \right) & \eta \left( \frac{\partial u_x}{\partial z} + \frac{\partial u_z}{\partial x} \right) \\
  \eta \left( \frac{\partial u_y}{\partial x} + \frac{\partial u_x}{\partial y} \right) & 2 \eta \frac{\partial u_y}{\partial y} - \frac{2}{3} \eta \left( \frac{\partial u_x}{\partial x} + \frac{\partial u_y}{\partial y} + \frac{\partial u_z}{\partial z} \right) & \eta \left( \frac{\partial u_y}{\partial z} + \frac{\partial u_z}{\partial y} \right) \\
  \eta \left( \frac{\partial u_z}{\partial x} + \frac{\partial u_x}{\partial z} \right) & \eta \left( \frac{\partial u_z}{\partial y} + \frac{\partial u_y}{\partial z} \right) & 2 \eta \frac{\partial u_z}{\partial z} - \frac{2}{3} \eta \left( \frac{\partial u_x}{\partial x} + \frac{\partial u_y}{\partial y} + \frac{\partial u_z}{\partial z} \right)
  \end{pmatrix}
\end{equation}

It is often useful to visualize the *second invariant* of the deviatoric stress tensor, which is independent of the coordinate reference frame and quantifies the local deviation of stress from a hydrostatic (non-convecting) state:

\begin{equation}
  \sigma_{\text{II}} = \sqrt{ \frac{1}{2} \left( \text{tr}(\sigma^{\prime 2}) - \text{tr}(\sigma^{\prime})^2 \right) }
  \label{eq:second-invariant-definition}
\end{equation}

where $\text{tr}(\sigma^{\prime 2}) = \Sigma \sigma^{\prime 2}_{ij}$ and $\text{tr}(\sigma^{\prime})^2 = (\sigma^{\prime}_{xx} + \sigma^{\prime}_{yy} + \sigma^{\prime}_{zz})^2$. Note that Equation \ref{eq:second-invariant-definition} uses the convention that compressive stress is positive. It follows from Equation \ref{eq:stress-complete} that the normal deviatoric stresses are:

\begin{equation}
  \begin{aligned}
    \sigma^{\prime}_{xx} &= \sigma_{xx} + P \\
    \sigma^{\prime}_{yy} &= \sigma_{yy} + P \\
    \sigma^{\prime}_{zz} &= \sigma_{zz} + P
  \end{aligned}
\end{equation}

and thus by definition $\text{tr}(\sigma^{\prime}) = \text{tr}(\sigma) + 3P = 0$, since $\text{tr}(\sigma) = -3P$. By this definition, Equation \ref{eq:second-invariant-definition} can be written as:

\begin{equation}
  \begin{aligned}
    \sigma_{\text{II}} &= \sqrt{ \frac{1}{2} \left( \text{tr}(\sigma^{\prime 2}) - 0 \right) } = \sqrt{ \frac{1}{2} \text{tr}(\sigma^{\prime 2}) } = \sqrt{ \frac{1}{2} \sum_{i,j}\sigma^{\prime 2}_{ij} } \\
    \sigma_{\text{II}} &= \sqrt{ \frac{1}{2} (\sigma^{\prime 2}_{xx} + \sigma^{\prime 2}_{yy} + \sigma^{\prime 2}_{zz} + \sigma^{\prime 2}_{xy} + \sigma^{\prime 2}_{yx} + \sigma^{\prime 2}_{xz} + \sigma^{\prime 2}_{zx} + \sigma^{\prime 2}_{yz} + \sigma^{\prime 2}_{zy}) } \\
    \sigma_{\text{II}} &= \sqrt{ \frac{1}{2} (\sigma^{\prime 2}_{xx} + \sigma^{\prime 2}_{yy} + \sigma^{\prime 2}_{zz}) + \sigma^{\prime 2}_{xy} + \sigma^{\prime 2}_{xz} + \sigma^{\prime 2}_{yz} }
  \end{aligned}
\end{equation}

Note also that many engineering applications use the von Mises stress:

\begin{equation}
  \sigma_{\text{vm}} = \sqrt{ \frac{3}{2} \sum_{i,j}\sigma^{\prime 2}_{ij} }
\end{equation}

which is proportional to the second invariant of the deviatoric stress tensor by a factor of $\sqrt{3}$:

\begin{equation}
  \sigma_{\text{vm}} = \sqrt{3} \, \sigma_{\text{II}}
\end{equation}

### The Continuity Equation {.unnumbered}

The continuity equation describes the conservation of mass (per volume) of the system through time:

\begin{equation}
  \frac{\partial \rho}{\partial t} + \nabla \cdot (\rho \vec{u}) = 0
  \label{eq:continuity-compressible-appendix}
\end{equation}

where $\rho$ is density, $t$ is time, and $\vec{u}$ is velocity. The form of this equation is particularly important for geodynamic simulations because density changes cause nonadiabatic pressure changes, which in turn cause changes in volume, which cause changes in density, which cause changes in nonadiabatic pressure. This feedback can lead to oscillations in the solution of mantle flow if the time derivative of density $\frac{\partial \rho}{\partial t}$ is kept in the equation and the timestep of the simulation is less than the viscous relaxation timescale [@curbelo2019].

The paper by @gassmoller2020 discuss the advantages and disadvantages of various treatments of the continuity equation that address numerical instability issues. @gassmoller2020 also propose new ways to solve the compressible momentum and continuity equations that can account for abrupt density changes due to strong heating/cooling and/or phase transformations.

### The Energy Equation {.unnumbered}

The energy equation accounts for temperature change of the fluid due to conductive and advective heat transfer, plus other source heating (or cooling) terms:

\begin{equation}
  \rho C_p \left( \frac{\partial T}{\partial t} + \vec{u} \cdot \nabla T \right) - \nabla \cdot \left( k \nabla T \right) = \sigma^{\prime} : \dot{\epsilon}^{\prime} + \alpha T \left( \vec{u} \cdot \nabla P \right) + Q_L
\end{equation}

where $\rho$ is density, $C_p$ is the specific heat capacity, $T$ is temperature, $t$ is time, $\vec{u}$ is velocity, $k$ is the thermal conductivity, $\alpha$ is the coefficient of thermal expansivity, $P$ is pressure, and $Q_L$ is the latent heat released or absorbed during phase transformations.

\cleardoublepage

## Formulations of Compressible Mantle Flow {.unnumbered #sec:formulations-appendix}

### The (Truncated) Anelastic Liquid Approximation {.unnumbered}

The anelastic liquid approximation (ALA) is based on the assumption that the density of Earth's mantle changes very little in the lateral direction (i.e., it is entirely depth-dependent). Density is assumed not to deviate very far from the adiabatic reference density $\bar{\rho}$ that is predefined in radial coordinates. Given an adiabatic reference condition, density is calculated via a first-order Taylor expansion with respect to pressure and temperature [@jarvis1980; @gassmoller2020]:

\begin{equation}
  \rho(P, T) \approx \bar{\rho} + \left( \frac{\partial \bar{\rho}}{\partial P} \right)_T \Delta P + \left( \frac{\partial \bar{\rho}}{\partial T} \right)_P \Delta T
  \label{eq:density-ala-expansion-appendix}
\end{equation}

and then rewriting Equation \ref{eq:density-ala-expansion-appendix} using standard thermodynamic relations $\beta = \frac{1}{\rho} \left( \frac{\partial \rho}{\partial P} \right)_T$ and $\alpha = -\frac{1}{\rho} \left( \frac{\partial \rho}{\partial T} \right)_P$ to get the expression:

\begin{equation}
  \rho(P, T) = \bar{\rho} \left( 1 + \bar{\beta_S} \hat{P} - \bar{\alpha} \hat{T} \right)
  \label{eq:density-ala-appendix}
\end{equation}

where $\bar{\rho}$, $\bar{\beta_S}$, $\bar{\alpha}$, are the adiabatic reference density, compressibility, and thermal expansivity, respectively, and $\Delta P = \hat{P} = P - \bar{P}$ and $\Delta T = \hat{T} = T - \bar{T}$ are the nonadiabatic pressure and temperature.

The only difference between the ALA and truncated ALA (TALA) is the assumption that deviations from the adiabatic reference density are only temperature-dependent:

\begin{equation}
  \rho(P, T) = \bar{\rho} \left(1 - \bar{\alpha} \hat{T} \right)
  \label{eq:density-tala}
\end{equation}

The other assumption is that deviations from the adiabatic reference density only apply to the gravitational body force term on the right-hand side of the momentum equation (Equation \ref{eq:navier-stokes-no-inertia-appendix}) and do not apply to the continuity equation or energy equation. Thus, the momentum and continuity equations simplify to:

\begin{equation}
  \nabla P - \nabla \cdot \sigma^{\prime} = \rho(P, T) g
\end{equation}

\begin{equation}
  \nabla \cdot (\bar{\rho} \vec{u}) = 0
  \label{eq:continuity-ala}
\end{equation}

### The Boussinesq Approximation {.unnumbered}

Similar to the ALA, the Boussinesq Approximation (BA) assumes that deviations from the adibatic reference density only apply to the right-hand side of the momentum equation (Equation \ref{eq:navier-stokes-no-inertia-appendix}). The BA simplifies the problem further, however, by assuming that the fluid is completely incompressible. In this case, the momentum and continuity equations become:

\begin{equation}
  \nabla P - \nabla \cdot \sigma^{\prime} = \rho(P, T) g
\end{equation}

\begin{equation}
  \nabla \cdot \vec{u} = 0
  \label{eq:continuity-incompressible}
\end{equation}

### The Isentropic Compression Approximation {.unnumbered}

The isentropic compression approximation (ICA) starts with the full continuity equation:

\begin{equation}
  \frac{\partial \rho}{\partial t} + \nabla \cdot (\rho \vec{u}) = 0
\end{equation}

which is expanded by applying the product rule to $\nabla \cdot (\rho \vec{u})$ and multiplying both sides by $\frac{1}{\rho}$ to get to the following expression:

\begin{equation}
  \frac{1}{\rho} \left(  \frac{\partial \rho}{\partial t} + \nabla \cdot (\rho \vec{u}) \right) = \frac{1}{\rho} \frac{\partial \rho}{\partial t} + \nabla \cdot \vec{u} + \left( \frac{1}{\rho} \nabla \rho \right) \cdot \vec{u} = 0
  \label{eq:continuity-expanded-appendix}
\end{equation}

The ICA is similar to the ALA, but replaces the adiabatic reference density $\bar{\rho}$ in the continuity equation (Equation \ref{eq:continuity-ala}) with an approximation of the local effects of compression [@gassmoller2020]. This is realized by substituting:

\begin{equation}
  \frac{1}{\rho} \nabla \rho \approx \frac{1}{\rho} \frac{\partial \rho}{\partial P} \nabla P \approx \beta_S \rho g
\end{equation}

where $\beta_S$, $\rho$, and $P$ are the compressibility, density, and pressure of the local fluid, respectively. The ICA can use an adiabatic reference condition in the continuity equation by setting $\beta_S = \bar{\beta_S}$ and $\rho = \bar{\rho}$, but this need not be the case. For example, using thermodynamic lookup tables as a material model allows $\beta_S$ and $\rho$ to vary as functions of pressure, temperature, and composition. Thus, the ICA is theoretically more accurate than the ALA---especially where materials are layered, undergoing phase transitions, or being strongly heated or cooled [@gassmoller2020].

The ICA is still not completely compressible, however, because it assumes that compaction is instantaneous by neglecting the viscous relaxation timescale introduced by the time derivative $\frac{\partial \rho}{\partial t}$ in the continuity equation (Equation \ref{eq:continuity-compressible-appendix}). Under these assumptions, the continuity equation becomes:

\begin{equation}
  \nabla\cdot \vec{u} = - \bar{\beta_S} \rho g \cdot \vec{u}
  \label{eq:continuity-ica-appendix}
\end{equation}

### The Projected Density Approximation {.unnumbered}

The PDA takes a different approach than the previously mentioned approximations because: 1) it includes the time derivative of density $\partial \rho / \partial t$, and 2) it does not define density as merely a function of PT (and PT derivatives $\partial \rho / \partial T$ and $\partial \rho / \partial P$), but rather uses density as defined in Equation \ref{eq:continuity-expanded-appendix}. In principle, this is a more accurate approach because density is free to evolve dynamically and be influenced by thermal and/or compositional changes. However, the challenge with using Equation \ref{eq:continuity-expanded-appendix} is that density is not a primary variable in the numerical solution of the governing equations. Practically speaking, only the time and spatial derivatives of PT and velocity are defined everywhere across the finite-element (FE) mesh, so the right-hand side of Equation \ref{eq:continuity-expanded-appendix} cannot be evaluated numerically [@gassmoller2020].

The solution proposed by @gassmoller2020 is to project the density field onto the FE mesh where it is discretized in the same manner as the PT and velocity fields, and as a result, Equation \ref{eq:continuity-expanded-appendix} is well-defined everywhere. The main disadvantage of this approach is that it introduces extra non-linearity by including the time derivative of density. This requires more non-linear iterations to converge on a solution—significantly slowing down the simulation in some cases [@gassmoller2020].

\cleardoublepage

## Viscoplastic Flow Laws {.unnumbered}

![Deformation maps of the Earth's upper mantle showing the relationships between temperature $T$, stress $\sigma$, strain rate $\dot{\epsilon}^{\prime}$, and active creep mechanism. The dry mantle is expected be be deforming by dislocation creep at typical values of $\dot{\epsilon}^{\prime}$ = 1e-15 s$^{-1}$and $T$ = 1600 K (a), while a hydrated mantle (or mantle with a very fine grain size) is expected to be deforming by diffusion creep at lower stress (b). At the same conditions, a dry mantle is expected to have a bulk viscosity on the order of 1e20 Pa s (c), while a wet mantle is expected bo be about an order of magnitude weaker (d). Models are based on the flow law from @karato1993: $\dot{\epsilon}^{\prime}_{\text{II}} = A \left( \frac{b}{d} \right)^m \left( \frac{\sigma_\text{II}}{G} \right)^n \exp \left( \frac{-E_a + V_a P}{R T} \right)$ assuming $P$ = 0 and a grain size $d$ of 1 mm. Other relevant parameters are: $G$ = 8e10 Pa is the bulk modulus, $b$ = 5e-10 m is the length of the Burgers vector, and $R$ = 8.314 J/K/mol is the gas constant.](../figs/deformation-map.png){#fig:deformation-map}

|Ref|Rock Type       |Creep Mechanism|     $A$|$m$|$n$| $E_a$| $V_a$|
|:--|:---------------|:--------------|-------:|--:|--:|-----:|-----:|
|1  |Dry mantle      |Diffusion      |1.50e-15|  3|1.0|3.75e5|  6e-6|
|1  |Wet mantle      |Diffusion      |2.50e-17|  3|1.0|3.75e5|  1e-5|
|1  |Dry mantle      |Dislusion      |1.10e-16|  0|3.5|5.30e5|6.5e-6|
|1  |Wet mantle      |Dislusion      |1.60e-18|  0|3.5|5.20e5|2.2e-5|
|2  |Dry olivine     |Dislocation    |2.50e-17|  0|3.5|5.32e5|     0|
|2  |Wet olivine     |Dislocation    |2.00e-21|  0|4.0|4.71e5|     0|
|2  |Rock salt       |Dislocation    |9.98e-32|  0|5.3|1.02e5|     0|
|2  |Quartz          |Dislocation    |1.00e-15|  0|2.0|1.67e5|     0|
|2  |Plag An75       |Dislocation    |2.08e-23|  0|3.2|2.38e5|     0|
|2  |Orthopyroxene   |Dislocation    |1.27e-15|  0|2.4|2.93e5|     0|
|2  |Clinopyroxene   |Dislocation    |3.94e-15|  0|2.6|3.35e5|     0|
|2  |Dry granite     |Dislocation    |1.14e-28|  0|3.2|1.23e5|     0|
|2  |Wet granite     |Dislocation    |7.96e-16|  0|1.9|1.37e5|     0|
|2  |Dry quartzite   |Dislocation    |2.67e-20|  0|2.4|1.56e5|     0|
|2  |Wet quartzite   |Dislocation    |5.07e-18|  0|2.3|1.54e5|     0|
|2  |Quartz Diorite  |Dislocation    |5.18e-18|  0|2.4|2.19e5|     0|
|2  |Diabase         |Dislocation    |7.96e-25|  0|3.4|2.60e5|     0|
|2  |Anorthosite     |Dislocation    |2.02e-23|  0|3.2|2.38e5|     0|
|2  |Felsic granulite|Dislocation    |2.01e-21|  0|3.1|2.43e5|     0|
|2  |Mafic granulite |Dislocation    |8.83e-22|  0|4.2|4.45e5|     0|

Table: Material parameters for a viscous flow law (Equation \ref{eq:viscous-creep}). Units are $A$: $\text{Pa}^{-n} \text{m}^{m} \text{s}^{-1}$, $m$: none, $n$: none, $V_a$: $\text{m}^3 \text{mol}^{-1}$, $E_a$ $\text{J} \text{mol}^{-1}$. References are 1: @hirth2003, 2: @ranalli1995. {#tbl:rheological-parameters}

\cleardoublepage

## Installing ASPECT {.unnumbered}

Native ASPECT installation is preferred to virtual machines or docker containers because it is faster. However, native installation requires configuring and installing a number of dependencies. This process will vary from machine to machine. Thankfully the ASPECT group has written [`candi`](https://github.com/dealii/candi), which tries to Compile AND Install (i.e., "candi") ASPECT's dependencies.

The base packages and libraries needed to run `candi` include:

  - `zlib`
  - `bzip2`
  - `git`
  - `cmake`
  - `boost`
  - `numdiff`
  - `openblas`
  - `scalapack`

These packages are normally pre-installed on your machine, for example via homebrew for MacOS, or `apt`, `yum`, `dnf`, for Linux. If you are using a HPC cluster, some of these packages and libraries might be loaded via `module load`. Once the base packages and libraries are installed and available in your environment, you can use `candi` to install the ASPECT dependencies. See the `candi` [github](https://github.com/dealii/candi) and ASPECT [wiki](https://github.com/geodynamics/aspect/wiki) pages for tips.

### Installing ASPECT on MacOS ARM {.unnumbered #sec:macos-config}

Follow the instructions at [hombrew's homepage](https://brew.sh/) to download and install Homebrew on your machine. Then use homebrew to install the base packages required by `candi`.

~~~{.bash}
brew install zlib bzip2 git cmake boost numdiff openblas scalapack
~~~

**NOTE:** If you already have homebrew, you will need to check if hdf5 is installed via homebrew. The version of hdf5 conflicts with the mpi-enabled hdf5 version installed by `candi`. If you have a homebrew version of hdf5, you will need to uninstall it (and the packages that depend on it) before using `candi`.

~~~{.bash}
brew uninstall hdf5
~~~

The following script will check for local or homebrew installations of the base packages, clone `candi`, and then try to install ASPECT and its dependencies. The script was written based off tips from the ASPECT [documentation](https://aspect-documentation.readthedocs.io/en/latest/user/install/local-installation/using-candi.html), ASPECT [wiki](https://github.com/geodynamics/aspect/wiki/Installation-on-ARM-OSX), and a lot of trial and error. The full ASPECT build takes about 1.5 hours on my Macbook air M1 (2020) with 8 processors.

~~~{#lst:install-dealii-macos .bash caption="Install script for deal.II on MacOS (tested on Macbook Air M1 2020, MacOS 15.3.1)."}
@include ../bash/install/dealii-macos.sh
~~~

~~~{#lst:install-aspect-macos .bash caption="Install script for ASPECT on MacOS (tested on Macbook Air M1 2020, MacOS 15.3.1)."}
@include ../bash/install/aspect-macos.sh
~~~

### Installing ASPECT on BARKLA2 {.unnumbered #sec:barkla2-config}

The following script will load the required modules, clone `candi`, and then try to install ASPECT and its dependencies. The script was written based off a lot of trial and error. The full ASPECT build takes 1.5 hours on BARKLA2 with 16 processors.

~~~{#lst:install-dealii-barkla2 .bash caption="Install script for deal.II on BARKLA2 (University of Liverpool)."}
@include ../bash/install/dealii-barkla2.sh
~~~

~~~{#lst:install-aspect-barkla2 .bash caption="Install script for ASPECT on BARKLA2 (University of Liverpool)."}
@include ../bash/install/aspect-barkla2.sh
~~~

\cleardoublepage

## Experimental Parameters {.unnumbered #sec:prm-configs}

### Experiment 2d_shell_bm {.unnumbered}

Code Block \ref{lst:base-ica-bm-simple-config} shows an example of an ASPECT configuration that uses a simple linear rheological model with material properties from an adiabatic reference condition. The configuration has been modified from the ASPECT cookbook:

  - [2D compressible convection with a reference profile and material properties from BurnMan](https://aspect-documentation.readthedocs.io/en/latest/user/cookbooks/cookbooks/burnman/doc/burnman.html)

~~~{#lst:base-ica-bm-simple-config .bash caption="Example of an ASPECT configuration with a simple rheological model."}
@include ../simulation/2d_shell/configs/full-mantle-with-piecewise-viscosity-profile.prm
~~~

\cleardoublepage

![](../figs/simulation/full_mantle_with_piecewise_viscosity_profile/tiles/temperature-nonadiabatic-density-nonadiabatic-viscosity-0020-full-view.png)

![](../figs/simulation/full_mantle_with_piecewise_viscosity_profile/tiles/strain-rate-pressure-nonadiabatic-sigma-ii-0020-full-view.png)

![Mantle convection model 2d_shell_bm after 100 Ma of evolution. Models were calculated with ASPECT using the ICA formulations of the continuity equation (Equation \ref{eq:continuity-ica}) and a simple temperature dependent rheological model (Equation \ref{eq:rheology-model-bm}). The adiabatic reference condition is shown in Figure \ref{fig:pyr-profiles}. Initial and final viscosities are shown in Figure \ref{fig:viscosity-profiles}.](../figs/simulation/full_mantle_with_piecewise_viscosity_profile/tiles/sigma-deviatoric-xx-sigma-deviatoric-yy-sigma-deviatoric-xy-0020-full-view.png){#fig:bm-mantle-flow-full-view}

\cleardoublepage

![](../figs/simulation/full_mantle_with_piecewise_viscosity_profile/tiles/temperature-nonadiabatic-density-nonadiabatic-viscosity-0020-zoomed-view.png)

![](../figs/simulation/full_mantle_with_piecewise_viscosity_profile/tiles/strain-rate-pressure-nonadiabatic-sigma-ii-0020-zoomed-view.png)

![Mantle convection model 2d_shell_bm after 100 Ma of evolution. Models were calculated with ASPECT using the ICA formulations of the continuity equation (Equation \ref{eq:continuity-ica}) and a simple temperature dependent rheological model (Equation \ref{eq:rheology-model-bm}). The adiabatic reference condition is shown in Figure \ref{fig:pyr-profiles}. Initial and final viscosities are shown in Figure \ref{fig:viscosity-profiles}.](../figs/simulation/full_mantle_with_piecewise_viscosity_profile/tiles/sigma-deviatoric-xx-sigma-deviatoric-yy-sigma-deviatoric-xy-0020-zoomed-view-1.png){#fig:bm-mantle-flow-zoomed-view}

\cleardoublepage

### Experiment 2d_shell_stb {.unnumbered}

Code Block \ref{lst:base-ica-stb-simple-config} shows an example of an ASPECT configuration that uses a more complex rheology after @steinberger2006 with material properties from a Perple_X lookup table. The configuration has been modified from the ASPECT cookbook:

  - [Convection using a pressure–temperature look-up table and the rheology of Steinberger and Calderwood (2006)](https://aspect-documentation.readthedocs.io/en/latest/user/cookbooks/cookbooks/steinberger/doc/steinberger.html)

~~~{#lst:base-ica-stb-simple-config .bash caption="Example of an ASPECT configuration with a steinberger rheological model."}
@include ../simulation/2d_shell/configs/full-mantle-with-steinberger-viscosity-profile.prm
~~~

\cleardoublepage

![](../figs/simulation/full_mantle_with_steinberger_viscosity_profile/tiles/temperature-nonadiabatic-density-nonadiabatic-viscosity-0047-full-view.png)

![](../figs/simulation/full_mantle_with_steinberger_viscosity_profile/tiles/strain-rate-pressure-nonadiabatic-sigma-ii-0047-full-view.png)

![Evolution of mantle convection model 2d_shell_stb. Models were calculated with ASPECT using the ICA formulations of the continuity equation (Equation \ref{eq:continuity-ica}) and a rheological model from @steinberger2006 (Equation \ref{eq:rheology-model-stb}). Modified from the ASPECT cookbook: [Convection using a pressure–temperature look-up table and the rheology of Steinberger and Calderwood (2006)](https://aspect-documentation.readthedocs.io/en/latest/user/cookbooks/cookbooks/steinberger/doc/steinberger.html).](../figs/simulation/full_mantle_with_steinberger_viscosity_profile/tiles/sigma-deviatoric-xx-sigma-deviatoric-yy-sigma-deviatoric-xy-0047-full-view.png){#fig:stb-mantle-flow-full-view}

\cleardoublepage

![](../figs/simulation/full_mantle_with_steinberger_viscosity_profile/tiles/temperature-nonadiabatic-density-nonadiabatic-viscosity-0047-zoomed-view.png)

![](../figs/simulation/full_mantle_with_steinberger_viscosity_profile/tiles/strain-rate-pressure-nonadiabatic-sigma-ii-0047-zoomed-view.png)

![Evolution of mantle convection model 2d_shell_stb. Models were calculated with ASPECT using the ICA formulations of the continuity equation (Equation \ref{eq:continuity-ica}) and a rheological model from @steinberger2006 (Equation \ref{eq:rheology-model-stb}). Modified from the ASPECT cookbook: [Convection using a pressure–temperature look-up table and the rheology of Steinberger and Calderwood (2006)](https://aspect-documentation.readthedocs.io/en/latest/user/cookbooks/cookbooks/steinberger/doc/steinberger.html).](../figs/simulation/full_mantle_with_steinberger_viscosity_profile/tiles/sigma-deviatoric-xx-sigma-deviatoric-yy-sigma-deviatoric-xy-0047-zoomed-view-1.png){#fig:stb-mantle-flow-zoomed-view}

\cleardoublepage

### Experiment 2d_shell_vp {.unnumbered}

Code Block \ref{lst:base-ica-vp-simple-config} shows an identical ASPECT configuration as Code Block \ref{lst:base-ica-vp-simple-config}, but with a composite viscoplastic rheology.

~~~{#lst:base-ica-vp-simple-config .bash caption="Example of an ASPECT configuration with a viscoplastic rheological model."}
@include ../simulation/2d_shell/configs/full-mantle-with-viscoplastic-viscosity-profile.prm
~~~

\cleardoublepage

![](../figs/simulation/full_mantle_with_viscoplastic_viscosity_profile/tiles/temperature-nonadiabatic-density-nonadiabatic-viscosity-0047-full-view.png)

![](../figs/simulation/full_mantle_with_viscoplastic_viscosity_profile/tiles/strain-rate-pressure-nonadiabatic-sigma-ii-0047-full-view.png)

![Evolution of mantle convection model 2d_shell_bm. Models were calculated with ASPECT using the ICA formulations of the continuity equation (Equation \ref{eq:continuity-ica}) and a viscoplastic rheological model (Equation \ref{eq:rheology-model-vp}) with material parameters after @vanderwiel2024 (Table \ref{tbl:rheological-parameters-vanderwiel20204}). Modified from the ASPECT cookbook: [Mantle convection with continents in an annulus](https://aspect-documentation.readthedocs.io/en/latest/user/cookbooks/cookbooks/mantle_convection_with_continents_in_annulus/doc/mantle_convection_in_annulus.html).](../figs/simulation/full_mantle_with_viscoplastic_viscosity_profile/tiles/sigma-deviatoric-xx-sigma-deviatoric-yy-sigma-deviatoric-xy-0047-full-view.png){#fig:vp-mantle-flow-full-view}

\cleardoublepage

![](../figs/simulation/full_mantle_with_viscoplastic_viscosity_profile/tiles/temperature-nonadiabatic-density-nonadiabatic-viscosity-0047-zoomed-view.png)

![](../figs/simulation/full_mantle_with_viscoplastic_viscosity_profile/tiles/strain-rate-pressure-nonadiabatic-sigma-ii-0047-zoomed-view.png)

![Evolution of mantle convection model 2d_shell_bm. Models were calculated with ASPECT using the ICA formulations of the continuity equation (Equation \ref{eq:continuity-ica}) and a viscoplastic rheological model (Equation \ref{eq:rheology-model-vp}) with material parameters after @vanderwiel2024 (Table \ref{tbl:rheological-parameters-vanderwiel20204}). Modified from the ASPECT cookbook: [Mantle convection with continents in an annulus](https://aspect-documentation.readthedocs.io/en/latest/user/cookbooks/cookbooks/mantle_convection_with_continents_in_annulus/doc/mantle_convection_in_annulus.html).](../figs/simulation/full_mantle_with_viscoplastic_viscosity_profile/tiles/sigma-deviatoric-xx-sigma-deviatoric-yy-sigma-deviatoric-xy-0047-zoomed-view-1.png){#fig:vp-mantle-flow-zoomed-view}

\cleardoublepage
