---
title: "Displaced and Faded"
subtitle: "How thermodynamics and reaction kinetics complicate seismic discontinuities"
# title: "Deflection and Sharpness of the 410km Seismic Discontinuity Largely Determined by Olivine -> Wadsleyite Reaction Kinetics"
# subtitle: "How Plumes and Slabs Displace Seismic Discontinuities"
author: ["Kerswell B.", "Wheeler J.", "Gassmöller R."]
date: "08 September 2025"
subject: "Mantle Convection"
keywords: [mantle convection, phase changes, geodynamics, numerical modeling]
lang: "en-US"
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

  1. Can displacements in the 410 km discontinuity be entirely explained by thermal effects? [e.g., @cottaar2016]
  2. How does compressibility affect displacement of the 410 km discontinuity? [e.g., @gassmoller2020]
  3. Do feedbacks exist between reaction kinetics and mantle flow that might explain observed displacements in the 410 km discontinuity? [e.g., @powell2018; @wheeler2014; @tajcmanova2015; @wheeler2018; @wheeler2020]
  4. Can we implement such feedbacks into geodynamic simulations of mantle convection? [e.g., @lu2024]
  5. Can we determine the timescales of such feedbacks?
  6. Can we determine the lengthscales of such feedbacks?
  7. Can we constrain reaction kinetic parameters with empirical seismic data?
  8. What are the implications of such feedbacks for mantle flow and observed seismic structures?

## Hypothesis Statements

  1. Deflection and sharpness of the 410 km discontinuity are inversely correlated with reaction rate
  2. Clearly-defined discontinuities require a reaction rate : strain rate ratio of > (?)
  3. Slab ponding at the 410 discontinuity requires a reaction rate : strain rate ratio of < (?)

\cleardoublepage

|Parameter name|Symbol|Unit|
|:--------------|:------|:----|
|Activation enthalpy|$H^{\ast}$|J mol$^{-1}$|
|Activation volume|$V^{\ast}$|m$^3$ mol$^{-1}$|
|Density|$\rho$|kg m$^{-3}$|
|Deviatoric stess tensor|$\sigma^{\prime}$|Pa|
|Deviatoric strain rate tensor|$\dot{\epsilon}^{\prime}$|s$^{-1}$|
|Dislocation density|$D$|m$^{-2}$|
|Effective deviatoric strain rate|$\dot{\epsilon}^{\prime}_{\text{II}}$|s$^{-1}$|
|Effective deviatoric stress|$\sigma_{\text{II}}$|Pa|
|Gas constant|$R$|J K$^{-1}$ mol$^{-1}$|
|Grain size|$d$|1 m|
|Gravitational acceleration|$g$|m s$^{-2}$|
|Interface geometry factor|$S$|m$^{-1}$|
|Interface growth kinetic prefactor|$B$|m s$^{-1} K^{-1}$|
|Interface growth rate|$\dot{x}$|m s$^{-1}$|
|Latent heat|$Q_L$|J kg$^{-1}$|
|Molar entropy|$\bar{S}$|J mol$^{-1}$ K$^{-1}$|
|Molar Gibbs free-energy|$\bar{G}$|J mol$^{-1}$|
|Molar volume|$\bar{V}$|m$^{3}$ mol$^{-1}$|
|Nominal background viscosity|$\eta_0$|Pa s|
|Nonadiabatic density|$\hat{\rho}$|kg m$^{-3}$|
|Nonadiabatic pressure|$\hat{P}$|Pa|
|Nonadiabatic temperature|$\hat{T}$|K|
|Phase transformation rate|$\dot{X}$|s$^{-1}$|
|Piecewise viscosity function|$\bar{\eta}$|Pa s|
|Pressure|$P$|Pa|
|Reference density|$\bar{\rho}$|kg m$^{-3}$|
|Reference compressibility|$\bar{\beta}$|Pa$^{-1}$|
|Reference pressure|$\bar{P}$|K|
|Reference specific heat capacity|$\bar{C}_p$|J kg$^{-1}$ K$^{-1}$|
|Reference temperature|$\bar{T}$|K|
|Reference thermal conductivity|$\bar{k}$|W m$^{-1}$ K$^{-1}$|
|Reference thermal expansivity|$\bar{\alpha}$|Pa$^{-1}$|
|Temperature|$T$|K|
|Thermal viscosity exponent|$A$|-|
|Time|$t$|s|
|Velocity|$\vec{u}$|m s$^{-1}$|
|Viscosity|$\eta$|Pa s|
|Volume fraction|$X$|-|
|Water content|$C_{OH}$|ppm|
|Water content exponent|$n$|-|

Table: Definition of symbols. {#tbl:symbol-definitions}

\cleardoublepage

# Methods {#sec:methods}

## Governing Equations for Compressible Mantle Flow {#sec:governing-equations}

Mantle flow is simulated by using the finite-element geodynamic code ASPECT [v3.0.0, @aspect-doi-v3.0.0; @aspectmanual; @heister2017; @kronbichler2012; @gassmoller2018; @clevenger2021; @fraters2019; @fraters2020] to find the velocity $\vec{u}$, pressure $P$, and temperature $T$ fields that satisfy the following equations:

\begin{equation}
  \nabla P - \nabla \cdot \sigma^{\prime} = \rho\, g
  \label{eq:navier-stokes-no-inertia}
\end{equation}

\begin{equation}
  \frac{\partial \rho}{\partial t} + \nabla \cdot (\rho\, \vec{u}) = 0
  \label{eq:continuity-compressible}
\end{equation}

\begin{equation}
  \rho\, \bar{C}_p \left(\frac{\partial T}{\partial t} + \vec{u} \cdot \nabla T \right) - \nabla \cdot \left(\bar{k}\, \nabla T \right) = \sigma^{\prime} : \dot{\epsilon}^{\prime} + \bar{\alpha}\, T \left(\vec{u} \cdot \nabla P \right) + Q_L
  \label{eq:energy}
\end{equation}

where $\sigma^{\prime}$ is the deviatoric stress tensor, $\rho$ is density, $g$ is gravitational acceleration, $t$ is time, $\bar{C}_p$, $\bar{k}$, $\bar{\alpha}$ are the adiabatic reference specific heat capacity, thermal conductivity, and thermal expansivity, respectively, and $Q_L$ is the latent heat released or absorbed during phase transformations. Equations \ref{eq:navier-stokes-no-inertia} and \ref{eq:continuity-compressible} together describe the buoyancy-driven flow of a nearly-incompressible isotropic fluid with negligible inertia and Equation \ref{eq:energy} describes the conduction, advection, and production (or consumption) of thermal energy [@schubert2001]. Note that the pressure $P$ in this context is equal to the mean normal stress of the solution and is positive under compression: $P = - \frac{\sigma_{xx} + \sigma_{yy}}{2}$ (see Appendix \ref{sec:momentum-derivation}).

The continuity equation is reformulated using the *projected density approximation* [PDA, @gassmoller2020] by applying the product rule to $\nabla \cdot (\rho\, \vec{u})$ and multiplying both sides of Equation \ref{eq:continuity-compressible} by $\frac{1}{\rho}$ to get to the following expression:

\begin{equation}
  \frac{1}{\rho} \left(\frac{\partial \rho}{\partial t} + \nabla \cdot (\rho\, \vec{u}) \right) = \frac{1}{\rho} \frac{\partial \rho}{\partial t} + \nabla \cdot \vec{u} + \left(\frac{1}{\rho} \nabla \rho \right) \cdot \vec{u} = 0
  \label{eq:continuity-expanded}
\end{equation}

This formulation of Equation \ref{eq:continuity-compressible} accounts for dynamic compression due to the combined local effects of pressure, temperature, and composition (see Appendix \ref{sec:formulations-appendix}). Our rationale for using the PDA formulation is that it is more accurate than an incompressible fluid in cases where materials are undergoing phase transitions and/or being strongly heated or cooled [@gassmoller2020].

## Numerical Setup {#sec:numerical-setup}

### Adiabatic Reference Conditions {#sec:adiabatic-reference-conditions}

In order to effectively converge on a solution, we needed to initialize our ASPECT simulations with reasonable guesses for the PT fields and material properties in Earth's upper mantle [@aspectmanual; @kronbichler2012; @heister2017]. For this purpose, we began by evaluating entropy changes over a PT range of 1573–1823 K and 0.001–20 GPa (Figure \ref{fig:isentrope}) using the Gibbs free-energy minimization software Perple_X [v.7.0.9, @connolly2009]. We assumed a dry pyrolitic bulk composition after @green1979 and phase equilibria were evaluated in the Na$_2$O‐CaO‐FeO‐MgO‐Al$_2$O$_3$‐SiO$_2$ (NCFMAS) chemical system with thermodynamic data and solution models of @stixrude2022. Equations of state were included for solid solution phases: olivine, plagioclase, spinel, clinopyroxene, wadsleyite, ringwoodite, perovskite, ferropericlase, high‐pressure C2/c pyroxene, orthopyroxene, akimotoite, post‐perovskite, Ca‐ferrite, garnet, and Na‐Al phase.

![Entropy (a) and density (b) changes in Earth's upper mantle under thermodynamic equilibrium. Material properties were computed with Perple_X using the equations of state and thermodynamic data of @stixrude2022. The black box indicates the approximate PT range of our ASPECT simulations, while the white line indicates the isentropic adiabat used to calculate material properties.](../figs/PYR-material-table.png){#fig:isentrope width=80%}

We then determined the mantle adiabat by applying the Newton–Raphson algorithm to find corresponding temperatures for each pressure such that entropy remains constant (white line in Figure \ref{fig:isentrope}). Material properties were evaluated at each PT point along the isentrope to construct the adiabatic reference conditions shown in Figure \ref{fig:material-property-profiles}. These reference conditions serve three main purposes: 1) initializing the PT fields and material properties in our ASPECT simulations (see Section \ref{sec:initialization-and-boundary-conditions}), 2) updating the material model during the ASPECT simulations (see Section \ref{sec:material-model}), and 3) serving as a basis for computing "nonadiabatic" quantities, such as the nonadiabatic temperature $\hat{T} = T - \bar{T}$ and nonadiabatic pressure $\hat{P} = P - \bar{P}$, and nonadiabatic density $\hat{\rho} = \rho - \bar{\rho}$, that quantify how much the solution deviates from a non-convecting ambient mantle.

![Reference material properties used in our ASPECT simulations. Profiles were computed using the BurnMan software [@cottaar2014; @myhill2023] and are based on the equations of state and thermodynamic data of @stixrude2022 for pure Mg$_{0.9}$Fe$_{0.01}$ olivine (ol) and wadsleyite (wad).](../figs/material-property-profiles.png){#fig:material-property-profiles}

### Initialization and Boundary Conditions {#sec:initialization-and-boundary-conditions}

Our ASPECT simulations were initialized with pure Mg$_{0.9}$Fe$_{0.1}$ olivine and wadsleyite within a 396 $\times$ 264 km rectangular model domain (Figure \ref{fig:initial-setup}). PT conditions of 10 GPa and 1706 K were applied at the top boundary such that the "surface" of the model is at approximately 280 km depth and the olivine $\Leftrightarrow$ wadsleyite transition occurs approximately half-way down the model. The initial PT fields were then computed by numerically integrating the following equations:

\begin{equation}
  \frac{d\bar{T}}{dy} = \frac{\bar{\alpha}\, \bar{T}\, g}{\bar{C}_p}
  \label{eq:adiabatic-temperature}
\end{equation}

\begin{equation}
  \frac{d\bar{P}}{dy} = \bar{\rho}\, g
  \label{eq:adiabatic-pressure}
\end{equation}

where $d\bar{P}/dy$ and $d\bar{T}/dy$ are vertical profiles of pressure and temperature applied uniformly across the model domain, $g$ is gravitational acceleration, and the depth-dependent adiabatic reference density $\bar{\rho}$, thermal expansivity $\bar{\alpha}$, and specific heat capacity $\bar{C}_p$ are determined from the material model (Figure \ref{fig:material-property-profiles}).

![Initial setup and boundary conditions for our ASPECT simulations. For slab models (top), a negative thermal anomaly and prescribed inflow of $\vec{u}_x$ = 0.5 and $\vec{u}_z$ -1.5 cm/yr are applied at the top boundary. For plume models (bottom), a positive thermal anomaly and prescribed inflow of $\vec{u}_y$ = 0.0 and $\vec{u}_y$ = 1.5 cm/yr are applied at the bottom boundary. All outflow boundaries have a constant vertical stress equal to the initial hydrostatic stress determined from Equation \ref{eq:adiabatic-pressure}. The top boundary has a constant PT of 10 GPa and 1706 K.](../figs/initial-setup.png){#fig:initial-setup width=65%}

For slab models, an initial negative thermal anomaly of 500 K is applied at the top boundary using a normal Guassian distribution, accompanied by a constant inflow velocity of $\vec{u}_x$ = 0.5 and $\vec{u}_y$ = -1.5 cm/yr (Figure \ref{fig:initial-setup}). All other boundaries (right, bottom, left) are prescribed a constant horizontal velocity of $\vec{u}_x$ = 0.5 cm/yr and constant vertical stress equal to the initial hydrostatic stress determined by numerical integration of Equation \ref{eq:adiabatic-pressure}. Boundary conditions for plume models essentially mirror those of slab models with two key changes: 1) plume models have zero horizontal velocity, and 2) a positive thermal anomaly of 500 K is applied at the bottom (inflow) boundary. These boundary conditions ensures that outflow of material from the model must be driven by dynamic (nonadiabatic) pressures generated by convection and/or volume changes due to the olivine $\Leftrightarrow$ wadsleyite phase transition.

### Material Model {#sec:material-model}

#### Material Properties {#sec:material-properties}

Material properties were updated during our ASPECT simulations by referencing the pressure-dependent adiabatic profiles shown in Figure \ref{fig:material-property-profiles}. Besides density, no PT corrections were applied to material properties---effectively assuming that deviations in material properties from a non-convecting ambient mantle were negligible. For density, however, we applied a nonadiabatic PT correction through a first-order Taylor expansion [@jarvis1980; @gassmoller2020]:

\begin{equation}
  \rho \approx \bar{\rho} + \left(\frac{\partial \bar{\rho}}{\partial P} \right)_T \Delta P + \left(\frac{\partial \bar{\rho}}{\partial T} \right)_P \Delta T
  \label{eq:density-ala-expansion}
\end{equation}

Equation \ref{eq:density-ala-expansion} is rewritten using standard thermodynamic relations $\beta = \frac{1}{\rho} \left(\frac{\partial \rho}{\partial P}\right)_T$ and $\alpha = -\frac{1}{\rho} \left(\frac{\partial \rho}{\partial T}\right)_P$ to get the expression:

\begin{equation}
  \rho = \bar{\rho} \left(1 + \bar{\beta}\, \hat{P} - \bar{\alpha}\, \hat{T} \right)
  \label{eq:density-ala}
\end{equation}

where $\bar{\rho}$, $\, \bar{\beta}$, $\, \bar{\alpha}$, are the adiabatic reference density, compressibility, and thermal expansivity, respectively, and $\Delta P = \hat{P} = P - \bar{P}$ and $\Delta T = \hat{T} = T - \bar{T}$ are the nonadiabatic PT. Note that the thermal conductivity $k$ = 4.0 Wm$^{-1}$K$^{-1}$ is constant in all our numerical experiments.

#### Phase Transition Kinetics {#sec:phase-transition-kinetics}

In our material model, the kinetics of the olivine $\Leftrightarrow$ wadsleyite phase transition are governed entirely by interface growth, as nucleation is assumed to saturate rapidly and does not limit the phase transformation [@cahn1956]. Following @faccenda2017, the transformed volume fraction is given by:

\begin{equation}
  X = 1 - \exp\left(-S\, \dot{x}\, t \right)
  \label{eq:volume-fraction}
\end{equation}

where $X$ is the volume fraction of the product phase (olivine or wadsleyite), $S$ is a geometric factor that accounts for nucleation sites, $\dot{x}$ is the interface growth rate, and $t$ is the elapsed time after site saturation. For inter-crystalline grain-boundary controlled growth, $S = 6.67/d$, where $d$ is grain size. For intra-crystalline dislocation-controlled growth, $S = 2\sqrt{D}$, where $D$ is dislocation density.

Since we assume growth-controlled kinetics, the expression for the interface growth rate $\dot{x}$ determines the overall transformation rate:

\begin{equation}
  \dot{x} = B\, T\, C_{OH}^n \exp\left(-\frac{H^{\ast} + P V^{\ast}}{R T}\right) \left(1 - \exp\left[-\frac{\Delta G}{R T}\right] \right)
  \label{eq:growth-rate}
\end{equation}

where $B$ is the interface growth kinetic prefactor, $C_{OH}$ is the concentration of water in the reactant phase, $n$ is the water content exponent, $H^{\ast}$ is activation enthalphy, $V^{\ast}$ is activation volume, $P$ is pressure, $T$ is temperature, $R$ is the gas constant, and $\Delta G$ is the excess Gibbs free-energy difference between olivine and wadsleyite. The excess Gibbs free-energy is approximated by:

\begin{equation}
  \Delta G \approx \Delta \bar{G} + \hat{P}\, \Delta \bar{V} - \hat{T}\, \Delta \bar{S}
  \label{eq:excess-gibbs}
\end{equation}

where $\Delta \bar{G}$, $\Delta \bar{V}$, and $\Delta \bar{S}$ are the Gibbs free-energy, volume, and entropy differences between olivine and wadsleyite along the adiabatic reference profile (Figure \ref{fig:phase-transition-kinetics-profile}), and $\hat{P}$ and $\hat{T}$ are the nonadiabatic PT.

![Reference thermodynamic properties used in our ASPECT simulations. Profiles were computed using the BurnMan software [@cottaar2014; @myhill2023] and are based on the equations of state and thermodynamic data of @stixrude2022 for pure Mg$_{0.9}$Fe$_{0.01}$ olivine (ol) and wadsleyite (wad).](../figs/phase-transition-kinetics-profile.png){#fig:phase-transition-kinetics-profile width=70%}

In this formulation, the time evolution of the olivine $\Leftrightarrow$ wadsleyite phase transformation is fully described by the interplay of pressure, temperature, and growth parameters (Table \ref{tbl:growth-parameters}), without explicit consideration of nucleation kinetics. Our material model therefore computes the macro-scale phase transformation rate by taking the time derivative of Equation \ref{eq:volume-fraction}:

\begin{equation}
  \frac{dX}{dt} = \dot{X} = S\, \dot{x}\, \left(1 - X \right)
  \label{eq:phase-transformation-rate}
\end{equation}

|Parameter|Value Range|Units|Reference|
|:--|--:|--:|--:|
|$B$|exp(-18)–exp(-21)|m s$^{-1}$ K$^{-1}$|@hosoya2005|
|$C_{OH}$|60–120|ppm|@hosoya2005|
|$n$|3.2||@hosoya2005|
|$H^{\ast}$|274e$^3$–300e$^3$|J mol$^{-1}$|@hosoya2005|
|$V^{\ast}$|3.3e$^{-6}$|m$^3$ mol$^{-1}$|@hosoya2005|
|$\Delta G$|Figure \ref{fig:phase-transition-kinetics-profile}|J mol$^{-1}$|@hosoya2005|

Table: Parameters used for the phase transition kinetics material model. {#tbl:growth-parameters}

Since the phase transformation rate $\dot{X}$ is slower than the advection timescale in our ASPECT simulations, we employ a first-order operator splitting scheme to decouple advection from phase-change kinetics. In this approach, the phase fraction is updated in two sequential steps within each overall time step $\Delta t$:

  1. **Advection step:** Solve the transport of material without phase changes over the time interval $\Delta t$ to yield an intermediate composition $X^\ast$:

     \begin{equation}
       \frac{\partial X}{\partial t} + \vec{u} \cdot \nabla X = 0
     \end{equation}

  2. **Reaction step:** integrate the growth-controlled kinetics over the same time interval (using a smaller sub-step $\delta t \le \Delta t$) to account for the transformation:

     \begin{equation}
       \frac{dX}{dt} = S\, \dot{x}\, \left(1 - X \right)
     \end{equation}

starting from $X^\ast$ to obtain the updated composition $X^{n+1}$. This operator splitting scheme ensures numerical stability while accurately capturing slow kinetics without restricting the convective timestep.

### Rheological Model {#sec:rheological-model}

Our ASPECT simulations use a simple rheological model where the mantle is assumed to have a nominal viscosity of 1e21 Pa s, which is modified by a depth-dependent piecewise function:

\begin{equation}
  \bar{\eta} =
  \begin{cases}
    1, & y \leq 132 \text{km} \\
    1, & y < 132 \text{km}
  \end{cases}
  \label{eq:piecewise-viscosity-function}
\end{equation}

This function produces a mantle with a uniform background viscosity across the olivine $\Leftrightarrow$ wadsleyite phase transition. A thermal dependency is then implemented through an exponential term:

\begin{equation}
  \eta = \bar{\eta}\, \eta_0 \exp \left(-A\, \frac{\hat{T}}{\bar{T}} \right)
  \label{eq:rheology-model-bm}
\end{equation}

where $\eta_0$ = 1e21 is the nominal background viscosity, $A$ = 1 is the thermal viscosity exponent factor, and $\hat{T}$ is the nonadiabatic temperature.

\cleardoublepage

# Results

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
  \rho \left(\frac{\partial \vec{u}}{\partial t} \right) = \nabla \cdot \sigma^{\prime} - \nabla P + \rho g
  \label{eq:navier-stokes-compressible}
\end{equation}

where $\rho$ is density, $\vec{u}$ is velocity, $P$ is pressure, $\sigma^{\prime}$ is the deviatoric stress tensor, and $g$ is gravitational acceleration. In terms of classical mechanics, the left-hand side is analogous to mass times acceleration $ma$, and the right-hand side are the forces $F$ that are acting on the fluid. Hence, the equation describes a balance between force and momentum $ma = F$ [@gerya2019].

The forces in Equation \ref{eq:navier-stokes-compressible} include the pressure gradient $\nabla P$ which acts to drive the fluid from high pressure to low pressure, the viscous forces $\nabla \cdot \sigma^{\prime}$ that dissipate energy by resisting flow, and the buoyancy forces $\rho g$ that drive convection (denser fluids sink while lighter fluids rise). Because the flow of Earth's mantle occurs at such slow rates, however, the inertial term $\left(\frac{\partial \vec{u}}{\partial t} \right)$ on the left-hand side of Equation \ref{eq:navier-stokes-compressible} can be ignored, and the momentum equation simplifies to:

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

where $\sigma_{ij}$ is total stress, $\sigma^{\prime}_{ij}$ is the deviatoric (non-hydrostatic) component of stress, $P = - \frac{\sigma_{xx} + \sigma_{yy} + \sigma_{zz}}{3}$ is the hydrostatic component of stress, and $\delta_{ij}$ is the Kroneker delta:

\begin{equation}
  \delta_{ij} =
  \begin{cases}
    1, & \text{if} i = j \\
    0, & \text{if} i \neq j
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

In practice, the deviatoric stress tensor $\sigma^{\prime}$ is computed by applying a constitutive relationship between stress and strain to express $\sigma^{\prime}$ in terms of velocity. In the present case, we apply a generalized linear model that combines shear deformation without rotation and volumetric deformation (dilation):

\begin{equation}
  \sigma^{\prime} = \eta \left(\nabla \vec{u} + \left(\nabla \vec{u} \right)^\intercal \right) - \left(\frac{2}{3} \eta - \zeta \right) \left(\nabla \cdot \vec{u} \right) I
\end{equation}

where $\vec{u}$ is velocity, $\eta$ is shear viscosity, and $\zeta$ is bulk viscosity. For nearly-incompressible fluids (very small $\zeta$), the expression for deviatoric stress simplifies to:

\begin{equation}
  \sigma^{\prime} = 2 \eta \dot{\epsilon}^{\prime}
\end{equation}

where $\dot{\epsilon}^{\prime} = \frac{1}{2} \left( \nabla \vec{u} + \left( \nabla \vec{u} \right)^\intercal \right) - \frac{1}{3} \left( \nabla \cdot \vec{u} \right) I$ is the deviatoric strain rate tensor. In full component form the deviatoric stress tensor is:

\begin{equation}
  \sigma^{\prime}_{ij} = \eta \left(\frac{\partial u_i}{\partial x_j} + \frac{\partial u_j}{\partial x_i} \right) - \frac{2}{3} \eta \left(\frac{\partial u_x}{\partial x} + \frac{\partial u_y}{\partial y} + \frac{\partial u_z}{\partial z} \right) \delta_{ij}
  \label{eq:stress-deviatoric-component}
\end{equation}

Note that the deviatoric stress and strain rate tensors are symmetric such that $\sigma^{\prime}_{ij} = \sigma^{\prime}_{ji}$ and $\dot{\epsilon}^{\prime}_{ij} = \dot{\epsilon}^{\prime}_{ji}$, which implies that there is zero rigid-body rotation in the fluid flow. Because of this symmetry, the full matrix form the deviatoric stress tensor can be written as:

\begin{equation}
  \sigma^{\prime} = \begin{pmatrix}
  2 \eta \frac{\partial u_x}{\partial x} - \frac{2}{3} \eta \left(\frac{\partial u_x}{\partial x} + \frac{\partial u_y}{\partial y} + \frac{\partial u_z}{\partial z} \right) & \eta \left(\frac{\partial u_x}{\partial y} + \frac{\partial u_y}{\partial x} \right) & \eta \left(\frac{\partial u_x}{\partial z} + \frac{\partial u_z}{\partial x} \right) \\
  \eta \left(\frac{\partial u_y}{\partial x} + \frac{\partial u_x}{\partial y} \right) & 2 \eta \frac{\partial u_y}{\partial y} - \frac{2}{3} \eta \left(\frac{\partial u_x}{\partial x} + \frac{\partial u_y}{\partial y} + \frac{\partial u_z}{\partial z} \right) & \eta \left(\frac{\partial u_y}{\partial z} + \frac{\partial u_z}{\partial y} \right) \\
  \eta \left(\frac{\partial u_z}{\partial x} + \frac{\partial u_x}{\partial z} \right) & \eta \left(\frac{\partial u_z}{\partial y} + \frac{\partial u_y}{\partial z} \right) & 2 \eta \frac{\partial u_z}{\partial z} - \frac{2}{3} \eta \left(\frac{\partial u_x}{\partial x} + \frac{\partial u_y}{\partial y} + \frac{\partial u_z}{\partial z} \right)
  \end{pmatrix}
\end{equation}

It is often useful to visualize the *second invariant* of the deviatoric stress tensor, which is independent of the coordinate reference frame and quantifies the local deviation of stress from a hydrostatic (non-convecting) state:

\begin{equation}
  \sigma_{\text{II}} = \sqrt{\frac{1}{2} \left(\text{tr}(\sigma^{\prime 2}) - \text{tr}(\sigma^{\prime})^2 \right)}
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
    \sigma_{\text{II}} &= \sqrt{\frac{1}{2} \left(\text{tr}(\sigma^{\prime 2}) - 0 \right)} = \sqrt{\frac{1}{2} \text{tr}(\sigma^{\prime 2})} = \sqrt{\frac{1}{2} \sum_{i, j}\sigma^{\prime 2}_{ij}} \\
    \sigma_{\text{II}} &= \sqrt{\frac{1}{2} (\sigma^{\prime 2}_{xx} + \sigma^{\prime 2}_{yy} + \sigma^{\prime 2}_{zz} + \sigma^{\prime 2}_{xy} + \sigma^{\prime 2}_{yx} + \sigma^{\prime 2}_{xz} + \sigma^{\prime 2}_{zx} + \sigma^{\prime 2}_{yz} + \sigma^{\prime 2}_{zy})} \\
    \sigma_{\text{II}} &= \sqrt{\frac{1}{2} (\sigma^{\prime 2}_{xx} + \sigma^{\prime 2}_{yy} + \sigma^{\prime 2}_{zz}) + \sigma^{\prime 2}_{xy} + \sigma^{\prime 2}_{xz} + \sigma^{\prime 2}_{yz}}
  \end{aligned}
\end{equation}

Note also that many engineering applications use the von Mises stress:

\begin{equation}
  \sigma_{\text{vm}} = \sqrt{\frac{3}{2} \sum_{i, j}\sigma^{\prime 2}_{ij}}
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
  \rho \bar{C}_p \left(\frac{\partial T}{\partial t} + \vec{u} \cdot \nabla T \right) - \nabla \cdot \left(k \nabla T \right) = \sigma^{\prime} : \dot{\epsilon}^{\prime} + \alpha T \left(\vec{u} \cdot \nabla P \right) + Q_L
\end{equation}

where $\rho$ is density, $\bar{C}_p$ is specific heat capacity, $T$ is temperature, $t$ is time, $\vec{u}$ is velocity, $k$ is thermal conductivity, $\alpha$ is thermal expansivity, $P$ is pressure, and $Q_L$ is the latent heat released or absorbed during phase transformations.

\cleardoublepage

## Formulations of Compressible Mantle Flow {.unnumbered #sec:formulations-appendix}

### The (Truncated) Anelastic Liquid Approximation {.unnumbered}

The anelastic liquid approximation (ALA) is based on the assumption that the density of Earth's mantle changes very little in the lateral direction (i.e., it is entirely depth-dependent). Density is assumed not to deviate very far from the adiabatic reference density $\bar{\rho}$ that is predefined in radial coordinates. Given an adiabatic reference condition, density is calculated via a first-order Taylor expansion with respect to pressure and temperature [@jarvis1980; @gassmoller2020]:

\begin{equation}
  \rho(P, T) \approx \bar{\rho} + \left(\frac{\partial \bar{\rho}}{\partial P} \right)_T \Delta P + \left(\frac{\partial \bar{\rho}}{\partial T} \right)_P \Delta T
  \label{eq:density-ala-expansion-appendix}
\end{equation}

and then rewriting Equation \ref{eq:density-ala-expansion-appendix} using standard thermodynamic relations $\beta = \frac{1}{\rho} \left( \frac{\partial \rho}{\partial P} \right)_T$ and $\alpha = -\frac{1}{\rho} \left( \frac{\partial \rho}{\partial T} \right)_P$ to get the expression:

\begin{equation}
  \rho(P, T) = \bar{\rho} \left(1 + \bar{\beta_S} \hat{P} - \bar{\alpha} \hat{T} \right)
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
  \frac{1}{\rho} \left(\frac{\partial \rho}{\partial t} + \nabla \cdot (\rho \vec{u}) \right) = \frac{1}{\rho} \frac{\partial \rho}{\partial t} + \nabla \cdot \vec{u} + \left(\frac{1}{\rho} \nabla \rho \right) \cdot \vec{u} = 0
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

These packages are normally pre-installed on your machine, for example via homebrew for macOS, or `apt`, `yum`, `dnf`, for Linux. If you are using an HPC cluster, some of these packages and libraries might be loaded via `module load`. Once the base packages and libraries are installed and available in your environment, you can use `candi` to install the ASPECT dependencies. See the `candi` [GitHub](https://github.com/dealii/candi) and ASPECT [wiki](https://github.com/geodynamics/aspect/wiki) pages for tips.

### Installing ASPECT on macOS ARM {.unnumbered #sec:macos-config}

Follow the instructions at [hombrew's homepage](https://brew.sh/) to download and install Homebrew on your machine. Then use homebrew to install the base packages required by `candi`.

~~~{.bash}
brew install zlib bzip2 git cmake boost numdiff openblas scalapack
~~~

**NOTE:** If you already have homebrew, you will need to check if hdf5 is installed via homebrew. The version of hdf5 conflicts with the mpi-enabled hdf5 version installed by `candi`. If you have a homebrew version of hdf5, you will need to uninstall it (and the packages that depend on it) before using `candi`.

~~~{.bash}
brew uninstall hdf5
~~~

The following script will check for local or homebrew installations of the base packages, clone `candi`, and then try to install ASPECT and its dependencies. The script was written based off tips from the ASPECT [documentation](https://aspect-documentation.readthedocs.io/en/latest/user/install/local-installation/using-candi.html), ASPECT [wiki](https://github.com/geodynamics/aspect/wiki/Installation-on-ARM-OSX), and a lot of trial and error. The full ASPECT build takes about 1.5 hours on my Macbook Air M1 (2020) with 8 processors.

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
