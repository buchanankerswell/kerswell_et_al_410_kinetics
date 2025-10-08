---
title: "Displaced and Faded"
subtitle: "How thermodynamics and kinetics collude to complicate seismic structures in Earth's mantle"
# title: "Displacement and Width of the 410km Seismic Discontinuity Largely Determined by Olivine -> Wadsleyite Reaction Kinetics"
# subtitle: "How Plumes and Slabs Displace Seismic Discontinuities"
author: ["Kerswell B.", "Wheeler J.", "Gassmöller R."]
date: "08 October 2025"
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

# Abstract {.unnumbered #sec:abstract}

The seismic expression of Earth’s 410 km discontinuity varies from sharp, high-amplitude interfaces to broad weak signals---patterns not fully explained by equilibrium thermodynamics. Laboratory studies show that the olivine $\Leftrightarrow$ wadsleyite transition is strongly rate-limited, implying that kinetics may strongly influence seismic structures, but uncertainties in growth and nucleation span several orders of magnitude. To quantify the effects of kinetics on seismic structures, we embed a growth-controlled kinetic model into compressible geodynamic simulations of mantle plumes and slabs. We ask how thermodynamic, kinetic, and dynamic feedbacks shape the 410, what timescales and lengthscales characterize these processes, and how seismic observables can be used to constrain kinetic parameters. Our results show that plumes follow smooth, monotonic trends: sluggish kinetics broaden and uplift the 410, whereas faster kinetics sharpen it to seismically resolvable widths. Slabs instead display thresholded behavior: moderately-slow kinetics generate metastable wedges, broad phase transition zones, and stagnation at the 410, while both ultra-slow ($\dot{X}$ < 0.02 Ma$^{-1}$) and ultra-fast ($\dot{X}$ > 1 Ma$^{-1}$) kinetics promote high-contrast seismic discontinuities. We identify kinetic thresholds of approximately $\dot{X}$ < 1 Ma$^{-1}$ that separate metastable, slab ponding regimes from equilibrated, slab penetrating ones, providing a direct link between kinetic rates and seismic signatures in subduction zones. These results reconcile observed lateral variability in 410 topography and sharpness across hotspots and subduction zones, establishing the 410 as a sensitive seismological probe of phase transition kinetics. More broadly, they demonstrate that incorporating kinetics into compressible geodynamic simulations is essential for connecting mineral-scale processes to mantle-scale dynamics, and for constraining phase transitions kinetics from seismic data.

\cleardoublepage

# Definition of Symbols {.unnumbered #sec:symbols}

|Parameter|Symbol|Unit|Equations|
|:--------------|:------|:----|:----|
|Activation enthalpy|$H^{\ast}$|J mol$^{-1}$|\ref{eq:growth-rate}|
|Activation volume|$V^{\ast}$|m$^3$ mol$^{-1}$|\ref{eq:growth-rate}|
|Compressibility (reference)|$\bar{\beta}$|Pa$^{-1}$|\ref{eq:density-ala}|
|Density|$\rho$|kg m$^{-3}$|\ref{eq:navier-stokes-no-inertia}–\ref{eq:continuity-expanded}, \ref{eq:density-ala-expansion}–\ref{eq:density-ala}|
|Density (reference)|$\bar{\rho}$|kg m$^{-3}$|\ref{eq:adiabatic-pressure}–\ref{eq:density-ala}|
|Density (dynamic)|$\hat{\rho}$|kg m$^{-3}$|-|
|Deviatoric stess tensor|$\sigma^{\prime}$|Pa|\ref{eq:navier-stokes-no-inertia}, \ref{eq:energy}|
|Deviatoric strain rate tensor|$\dot{\epsilon}^{\prime}$|s$^{-1}$|\ref{eq:energy}|
|Dislocation density|$D$|m$^{-2}$|\ref{eq:growth-rate}|
|Gas constant|$R$|J K$^{-1}$ mol$^{-1}$|\ref{eq:growth-rate}|
|Grain size|$d$|1 m|\ref{eq:growth-rate}|
|Gravitational acceleration|$g$|m s$^{-2}$|\ref{eq:navier-stokes-no-inertia}, \ref{eq:adiabatic-temperature}–\ref{eq:adiabatic-pressure}|
|Growth rate|$\dot{x}$|m s$^{-1}$|\ref{eq:growth-rate}|
|Kinetic prefactor|$A$|m s$^{-1}$ K$^{-1}$ ppm$_{OH}^{-n}$|\ref{eq:growth-rate}|
|Latent heat|$Q_L$|J kg$^{-1}$|\ref{eq:energy}|
|Molar entropy|$\bar{S}$|J mol$^{-1}$ K$^{-1}$|\ref{eq:excess-gibbs}|
|Molar Gibbs free-energy|$\bar{G}$|J mol$^{-1}$|\ref{eq:excess-gibbs}|
|Molar volume|$\bar{V}$|m$^{3}$ mol$^{-1}$|\ref{eq:excess-gibbs}|
|Nucleation site factor|$S$|m$^{-1}$|\ref{eq:growth-rate}|
|Phase transition rate|$\dot{X}$|s$^{-1}$|\ref{eq:phase-transition-rate}|
|Pressure|$P$|Pa|\ref{eq:navier-stokes-no-inertia}, \ref{eq:energy}, \ref{eq:growth-rate}|
|Pressure (reference)|$\bar{P}$|K|\ref{eq:adiabatic-pressure}|
|Pressure (dynamic)|$\hat{P}$|Pa|\ref{eq:density-ala}, \ref{eq:excess-gibbs}|
|Specific heat capacity (reference)|$\bar{C}_p$|J kg$^{-1}$ K$^{-1}$|\ref{eq:energy}, \ref{eq:adiabatic-temperature}|
|Temperature|$T$|K|\ref{eq:energy}, \ref{eq:growth-rate}|
|Temperature (reference)|$\bar{T}$|K|\ref{eq:adiabatic-temperature}, \ref{eq:rheological-model}|
|Temperature (dynamic)|$\hat{T}$|K|\ref{eq:density-ala}, \ref{eq:excess-gibbs}, \ref{eq:rheological-model}|
|Thermal conductivity (reference)|$\bar{k}$|W m$^{-1}$ K$^{-1}$|\ref{eq:energy}|
|Thermal expansivity (reference)|$\bar{\alpha}$|Pa$^{-1}$|\ref{eq:energy}, \ref{eq:adiabatic-temperature}, \ref{eq:density-ala}|
|Time|$t$|s|\ref{eq:continuity-compressible}–\ref{eq:continuity-expanded}, \ref{eq:volume-fraction}, \ref{eq:phase-transition-rate}|
|Velocity|$\vec{u}$|m s$^{-1}$|\ref{eq:continuity-compressible}–\ref{eq:continuity-expanded}, \ref{eq:composition}|
|Viscosity|$\eta$|Pa s|\ref{eq:rheological-model}|
|Viscosity (reference)|$\bar{\eta}$|Pa s|\ref{eq:rheological-model}|
|Viscosity exponent|$B$|-|\ref{eq:rheological-model}|
|Volume fraction|$X$|-|\ref{eq:volume-fraction},  \ref{eq:phase-transition-rate}–\ref{eq:composition}|
|Water content|$C_{OH}$|ppm|\ref{eq:growth-rate}|
|Water content exponent|$n$|-|\ref{eq:growth-rate}|

\cleardoublepage

# Introduction {#sec:introduction}

The mantle transition zone (MTZ) hosts two major seismic discontinuities near 410 and 660 km depth, linked to polymorphic phase transitions of olivine and garnet [@ringwood1975; @katsura2004]. Global seismic studies reveal that the depth, sharpness, and amplitude of the 410 vary widely: some regions display sharp, high-amplitude discontinuities, whereas others exhibit broad, displaced, or weakened signals [@deuss2009; @lawrence2008; @schmerr2007]. Such variability cannot be explained by equilibrium thermodynamics alone, which predicts discontinuity topography primarily from temperature and Clapeyron slopes. Additional factors---notably phase transition kinetics, water content, and compressibility---likely play central roles [@rubie1994; @perrillat2016; @fukao2013].

Mineral-physics experiments demonstrate that the olivine $\Leftrightarrow$ wadsleyite transition is strongly rate-limited, with growth and nucleation controlled by temperature, pressure, water, and deformation state [@kubo2004; @ohuchi2022; @ledoux2023]. Under cold slab conditions, sluggish kinetics can permit metastable persistence of olivine well below 410 km, influencing slab stagnation and deep seismicity [@fukao2009]. In hot plume environments, slow kinetics can broaden and uplift the discontinuity, potentially explaining reduced amplitudes beneath hotspots [@chambers2005; @jenkins2016]. Yet kinetic parameters remain poorly constrained, spanning several orders of magnitude, and feedbacks introduced by coupling kinetic rate laws to compressible mantle flow simulations are incompletely understood.

Here we address this gap by embedding a growth-controlled kinetic model within compressible geodynamic simulations of plumes and slabs using the numerical geodynamic software ASPECT. We ask: 1) how do feedbacks among thermodynamics, kinetics, and compressibility shape the 410? 2) What are the characteristic timescales and lengthscales of these feedbacks? 3) What are the implications for mantle flow and seismic expression? And 4) can seismic observations be used to constrain effective kinetic parameters?

We explore these questions through a systematic suite of numerical simulations that vary kinetic parameters across the experimentally constrained ranges of uncertainties. We evaluate the resulting phase transition zone (PTZ) displacements, widths, and seismic velocity contrasts. Our results show that plumes display predictable power-law scalings between kinetics and PTZ structure, while slabs exhibit abrupt thresholds separating metastable ponding from penetrative behaviour, and broad PTZs from sharp seismic discontinuities. Together, these findings demonstrate that phase transition kinetics provide a first-order control on MTZ structure, help reconcile observed 410 variability, and offer a pathway to constrain phase transition rates and kinetic parameters from seismic data.

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

where $\sigma^{\prime}$ is the deviatoric stress tensor, $\rho$ is density, $g$ is gravitational acceleration, $t$ is time, $\bar{C}_p$, $\bar{k}$, $\bar{\alpha}$ are the reference specific heat capacity, thermal conductivity, and thermal expansivity, respectively (see Section \ref{sec:adiabatic-reference-conditions}), and $Q_L$ is the latent heat released or absorbed during phase transitions. Equations \ref{eq:navier-stokes-no-inertia} and \ref{eq:continuity-compressible} together describe the buoyancy-driven flow of a nearly-incompressible isotropic fluid with negligible inertia and Equation \ref{eq:energy} describes the conduction, advection, and production (or consumption) of thermal energy [@schubert2001]. Note that the pressure $P$ in this context is equal to the mean normal stress and is positive under compression: $P = - \frac{\sigma_{xx} + \sigma_{yy}}{2}$ (see Appendix \ref{sec:momentum-derivation}).

The fully compressible form of the continuity equation above (Equation \ref{eq:continuity-compressible}) can result in a dynamic feedback between density changes and pressures that cause numerical oscillations if the the advection timestep is less than the viscous relaxation timescale [@curbelo2019]. Therefore the continuity equation is reformulated using the *projected density approximation* [PDA, @gassmoller2020] by applying the product rule to $\nabla \cdot (\rho\, \vec{u})$ and multiplying both sides of Equation \ref{eq:continuity-compressible} by $\frac{1}{\rho}$ to get to the following expression:

\begin{equation}
  \frac{1}{\rho} \frac{\partial \rho}{\partial t} + \nabla \cdot \vec{u} + \left(\frac{1}{\rho} \nabla \rho \right) \cdot \vec{u} = 0
  \label{eq:continuity-expanded}
\end{equation}

Adopting the PDA allows us to retain the coupling between density, pressure, temperature, and phase transitions, while mitigating the numerical instabilities associated with the fully compressible form of the continuity equation [@gassmoller2020]. The PDA formulation is therefore particularly appropriate for our models, which include dynamic thermal effects as well as an olivine $\Leftrightarrow$ wadsleyite phase transition that would otherwise be neglected by incompressible forms of the continuity equation (e.g., Boussinesq). This formulation unique to ASPECT provides a stable and internally consistent framework for solving the governing equations---ensuring that compressibility-related feedbacks were accurately represented without compromising numerical robustness.

## Numerical Setup {#sec:numerical-setup}

### Adiabatic Reference Conditions {#sec:adiabatic-reference-conditions}

In order to effectively converge on a solution, we needed to initialize our ASPECT simulations with reasonable guesses for the pressure-temperature (PT) fields and material properties in Earth's upper mantle [@aspectmanual; @kronbichler2012; @heister2017]. For this purpose, we began by evaluating entropy changes over a PT range of 1573–1973 K and 0.001–25 GPa (Figure \ref{fig:isentrope}) using the Gibbs free-energy minimization software Perple_X [v.7.0.9, @connolly2009]. We assumed a dry pyrolitic bulk composition after @green1979 and phase equilibria were evaluated in the Na$_2$O‐CaO‐FeO‐MgO‐Al$_2$O$_3$‐SiO$_2$ (NCFMAS) chemical system with thermodynamic data and solution models of @stixrude2022. Equations of state were included for solid solution phases: olivine, plagioclase, spinel, clinopyroxene, wadsleyite, ringwoodite, perovskite, ferropericlase, high‐pressure C2/c pyroxene, orthopyroxene, akimotoite, post‐perovskite, Ca‐ferrite, garnet, and Na‐Al phase.

![Entropy (a) and density (b) changes in Earth's upper mantle under thermodynamic equilibrium. Material properties were computed with Perple_X using the equations of state and thermodynamic data of @stixrude2022. The black box indicates the approximate PT range of our ASPECT simulations, while the white line indicates the isentropic adiabat used to calculate reference material properties.](../figs/PYR-material-table.png){#fig:isentrope width=75%}

We then determined the mantle adiabat by applying the Newton–Raphson algorithm to find corresponding temperatures for each pressure such that entropy remains constant (white line in Figure \ref{fig:isentrope}). Material properties were evaluated at each PT point along the isentrope to construct the adiabatic reference conditions shown in Figure \ref{fig:material-property-profiles}. These reference conditions serve three main purposes: 1) initializing the PT fields and material properties in our ASPECT simulations (see Section \ref{sec:initialization-and-boundary-conditions}), 2) updating the material model during the simulations (see Section \ref{sec:material-model}), and 3) serving as a basis for computing "dynamic" quantities, such as the dynamic temperature $\hat{T} = T - \bar{T}$, dynamic pressure $\hat{P} = P - \bar{P}$, and dynamic density $\hat{\rho} = \rho - \bar{\rho}$, that quantify how much the approximate numerical solution deviates from the reference adiabatic conditions (i.e., a non-convecting ambient mantle).

![Reference material properties used in our ASPECT simulations. Profiles were computed using the BurnMan software [@cottaar2014; @myhill2023] and were based on the equations of state and thermodynamic data of @stixrude2022 for pure Mg$_{0.9}$Fe$_{0.01}$ olivine (ol) and wadsleyite (wad).](../figs/material-property-profiles.png){#fig:material-property-profiles}

### Initialization and Boundary Conditions {#sec:initialization-and-boundary-conditions}

Our ASPECT simulations were initialized with pure Mg$_{0.9}$Fe$_{0.1}$ olivine and wadsleyite within a 396 $\times$ 264 km rectangular model domain (Figure \ref{fig:initial-setup}). "Surface" PT conditions of 10 GPa and 1706 K were applied at the top boundary such that the olivine $\Leftrightarrow$ wadsleyite transition occurs approximately half-way down the model. The initial PT fields were then computed by numerically integrating the following equations:

\begin{equation}
  \frac{d\bar{T}}{dy} = \frac{\bar{\alpha}\, \bar{T}\, g}{\bar{C}_p}
  \label{eq:adiabatic-temperature}
\end{equation}

\begin{equation}
  \frac{d\bar{P}}{dy} = \bar{\rho}\, g
  \label{eq:adiabatic-pressure}
\end{equation}

where $d\bar{P}/dy$ and $d\bar{T}/dy$ are adiabatic PT profiles applied uniformly across the model domain, $g$ is gravitational acceleration, and the density $\bar{\rho}$, thermal expansivity $\bar{\alpha}$, and specific heat capacity $\bar{C}_p$ were determined from the adiabatic reference conditions shown in Figure \ref{fig:material-property-profiles}. Normal Guassian-shaped thermal anomalies of $\leq$ 500 K were then applied at the top and bottom boundaries for slab and plume simulations, respectively.

![Initial setup and boundary conditions for slab (top) and plume (bottom) simulations. The thermal anomalies and prescribed inflows for slab and plume models were essentially mirrored, except for a horizontal boundary velocity of $\vec{u}_x$ = 0.5 cm/yr applied to slab models. Traction boundary conditions ($\sigma_{xy}$ and $\sigma_{yy}$) ensure that outflows must be driven by dynamic pressures. The top boundary has a constant "surface" PT of 10 GPa and 1706 K.](../figs/initial-setup.png){#fig:initial-setup width=70%}

Velocity and traction boundary conditions were set to ensures that outflow of material from the model must be driven by dynamic pressures generated by convection and/or volume changes due to the olivine $\Leftrightarrow$ wadsleyite phase transition. For slab models, a prescribed inflow velocity of $\vec{u}_x$ = 0.5 and $\vec{u}_y$ = -1.5 cm/yr was applied at the top boundary. All other boundaries (right, bottom, left) were prescribed a constant horizontal velocity of $\vec{u}_x$ = 0.5 cm/yr and constant vertical stress equal to the initial hydrostatic pressure $\bar{P}$ determined by numerical integration of Equation \ref{eq:adiabatic-pressure}. Velocity and traction boundary conditions for plume models essentially mirror those of slab models with one key change: plume models have zero horizontal velocity at all boundaries.

### Material Model {#sec:material-model}

#### Material Properties {#sec:material-properties}

Material properties were updated during our ASPECT simulations by referencing the adiabatic reference conditions shown in Figure \ref{fig:material-property-profiles}. Besides density, no PT corrections were applied to material properties, effectively assuming that deviations in material properties from a non-convecting ambient mantle were negligible. For density, however, we applied a dynamic PT correction through a first-order Taylor expansion [@jarvis1980; @gassmoller2020]:

\begin{equation}
  \rho \approx \bar{\rho} + \left(\frac{\partial \bar{\rho}}{\partial P} \right)_T \Delta P + \left(\frac{\partial \bar{\rho}}{\partial T} \right)_P \Delta T
  \label{eq:density-ala-expansion}
\end{equation}

Equation \ref{eq:density-ala-expansion} is rewritten using standard thermodynamic relations $\beta = \frac{1}{\rho} \left(\frac{\partial \rho}{\partial P}\right)_T$ and $\alpha = -\frac{1}{\rho} \left(\frac{\partial \rho}{\partial T}\right)_P$ to get the expression:

\begin{equation}
  \rho = \bar{\rho} \left(1 + \bar{\beta}\, \hat{P} - \bar{\alpha}\, \hat{T} \right)
  \label{eq:density-ala}
\end{equation}

where $\bar{\rho}$, $\, \bar{\beta}$, $\, \bar{\alpha}$, are the adiabatic reference density, compressibility, and thermal expansivity, respectively, and $\Delta P = \hat{P} = P - \bar{P}$ and $\Delta T = \hat{T} = T - \bar{T}$ are the dynamic PT. Note that the reference thermal conductivity $\bar{k}$ = 4.0 Wm$^{-1}$K$^{-1}$ is constant in all our numerical experiments.

#### Phase Transition Kinetics {#sec:phase-transition-kinetics}

The kinetics of the olivine $\Leftrightarrow$ wadsleyite phase transition were governed entirely by interface growth, as nucleation was assumed to saturate rapidly and did not limit the phase transition [@cahn1956; @hosoya2005]. Following @faccenda2017, the transformed volume fraction is given by:

\begin{equation}
  X = 1 - \exp\left(-S\, \dot{x}\, t \right)
  \label{eq:volume-fraction}
\end{equation}

where $X$ is the volume fraction of the product phase (olivine or wadsleyite), $S$ is a geometric factor that accounts for nucleation sites, $\dot{x}$ is the interface growth rate, and $t$ is the elapsed time after site saturation. For inter-crystalline grain-boundary controlled growth, $S = 6.67/d$, where $d$ is grain size. For intra-crystalline dislocation-controlled growth, $S = 2\sqrt{D}$, where $D$ is dislocation density [after @mohiuddin2018].

Since we assumed growth-controlled kinetics, the expression for the interface growth rate $\dot{x}$ determined the overall phase transition rate [@turnbull1956]:

\begin{equation}
  \dot{x} = A\, T\, C_{OH}^n \exp\left(-\frac{H^{\ast} + P V^{\ast}}{R\, T}\right) \left(1 - \exp\left[-\frac{\Delta G}{R\, T}\right] \right)
  \label{eq:growth-rate}
\end{equation}

where $A$ is a kinetic prefactor, $C_{OH}$ is the concentration of water in the reactant phase, $n$ is the water content exponent, $H^{\ast}$ is activation enthalphy, $V^{\ast}$ is activation volume, $P$ is pressure, $T$ is temperature, $R$ is the gas constant, and $\Delta G$ is the excess Gibbs free-energy difference between olivine and wadsleyite. The excess Gibbs free-energy $\Delta G$ is approximated by:

\begin{equation}
  \Delta G \approx \Delta \bar{G} + \hat{P}\, \Delta \bar{V} - \hat{T}\, \Delta \bar{S}
  \label{eq:excess-gibbs}
\end{equation}

where $\Delta \bar{G}$, $\Delta \bar{V}$, and $\Delta \bar{S}$ are the Gibbs free-energy, volume, and entropy differences between olivine and wadsleyite along the adiabatic reference profile (Figure \ref{fig:phase-transition-kinetics-profile}), and $\hat{P}$ and $\hat{T}$ are the dynamic PT.

In this formulation, the time evolution of the olivine $\Leftrightarrow$ wadsleyite phase transition is fully described by the interplay of pressure, temperature, and growth parameters (Table \ref{tbl:growth-parameters}), without explicit consideration of nucleation kinetics [@cahn1956; @rubie1994; @hosoya2005; @faccenda2017]. The macro-scale olivine $\Leftrightarrow$ wadsleyite phase transition rate was therefore computed by taking the time derivative of Equation \ref{eq:volume-fraction}:

\begin{equation}
  \frac{dX}{dt} = \dot{X} = S\, \dot{x}\, \left(1 - X \right)
  \label{eq:phase-transition-rate}
\end{equation}

#### Operator Splitting {#sec:operator-splitting}

Since the phase transition rate $\dot{X}$ is slower than the advection timescale in our ASPECT simulations, we employ a first-order operator splitting scheme to decouple advection from phase-change kinetics. In this approach, the phase fraction is updated in two sequential steps within each overall time step $\Delta t$:

\begin{equation}
  \frac{\partial X}{\partial t} + \vec{u} \cdot \nabla X = 0
  \label{eq:composition}
\end{equation}

  1. **Advection step:** Solve the transport of material $\left(\vec{u} \cdot \nabla X \right)$ without phase changes over the time interval $\Delta t$ to yield an intermediate composition $X^\ast$
  2. **Reaction step:** starting from $X^\ast$, integrate Equation \ref{eq:phase-transition-rate} over the same time interval using a smaller sub-step $\delta t \le \Delta t$ to obtain the updated composition $X^{n+1}$

This operator splitting scheme ensures numerical stability while accurately capturing slow kinetics without restricting the convective timestep [@aspectmanual].

![Reference thermodynamic properties used in our ASPECT simulations. Profiles were computed using the BurnMan software [@cottaar2014; @myhill2023] and are based on the equations of state and thermodynamic data of @stixrude2022 for pure Mg$_{0.9}$Fe$_{0.01}$ olivine (ol) and wadsleyite (wad).](../figs/phase-transition-kinetics-profile.png){#fig:phase-transition-kinetics-profile width=80%}

### Rheological Model {#sec:rheological-model}

Our ASPECT simulations use a simple rheological model where mantle viscosity is modified by a depth-dependent piecewise function:

\begin{equation}
  \bar{\eta} =
  \begin{cases}
    1\, \eta_0, & y \leq 132 \text{km} \\
    1\, \eta_0, & y > 132 \text{km}
  \end{cases}
  \label{eq:piecewise-viscosity-function}
\end{equation}

where $\eta_0$ is the nominal background viscosity of the upper mantle [@ranalli1995; @karato2008]. A thermal dependency is then implemented through an exponential term:

\begin{equation}
  \eta = \bar{\eta} \exp \left(-B\, \frac{\hat{T}}{\bar{T}} \right)
  \label{eq:rheological-model}
\end{equation}

where $B$ is the thermal viscosity exponent factor, and $\hat{T}$ is the dynamic temperature.

| $A$       | $H^{\ast}$ | $V^{\ast}$ | $d$   | $C_{OH}$  | $n$ | $\eta_0$ | $B$ |
|:----------|:-----------|:-----------|:------|:----------|:----|:---------|:----|
| $e^{-18}$ | 274        | 3.3e-6     | 0.005 | 50        | 2.6 | 1e21     | 1   |
| $e^{-18}$ | 274        | 3.3e-6     | 0.001 | 50        | 3.2 | 1e21     | 1   |
| $e^{-18}$ | 274        | 3.3e-6     | 0.005 | 50        | 3.2 | 1e21     | 1   |
| $e^{-18}$ | 274        | 3.3e-6     | 0.01  | 50        | 3.2 | 1e21     | 1   |
| $e^{-18}$ | 274        | 3.3e-6     | 0.005 | 50        | 3.8 | 1e21     | 1   |
| $e^{-18}$ | 274        | 3.3e-6     | 0.005 | 500       | 2.6 | 1e21     | 1   |
| $e^{-18}$ | 274        | 3.3e-6     | 0.001 | 500       | 3.2 | 1e21     | 1   |
| $e^{-18}$ | 274        | 3.3e-6     | 0.005 | 500       | 3.2 | 1e21     | 1   |
| $e^{-18}$ | 274        | 3.3e-6     | 0.01  | 500       | 3.2 | 1e21     | 1   |
| $e^{-18}$ | 274        | 3.3e-6     | 0.005 | 500       | 3.8 | 1e21     | 1   |
| $e^{-18}$ | 274        | 3.3e-6     | 0.005 | 5000      | 2.6 | 1e21     | 1   |
| $e^{-18}$ | 274        | 3.3e-6     | 0.001 | 5000      | 3.2 | 1e21     | 1   |
| $e^{-18}$ | 274        | 3.3e-6     | 0.005 | 5000      | 3.2 | 1e21     | 1   |
| $e^{-18}$ | 274        | 3.3e-6     | 0.01  | 5000      | 3.2 | 1e21     | 1   |
| $e^{-18}$ | 274        | 3.3e-6     | 0.005 | 5000      | 3.8 | 1e21     | 1   |

Table: Full list of kinetic and rheological parameters used in plume and slab simulations. Kinetic parameter values are consistent with the range of experimental data from @hosoya2005. Units are $A$: m s$^{-1}$ K$^{-1}$ ppm$_{OH}^{-n}$, $H^{\ast}$: J mol$^{-1}$, $V^{\ast}$: m$^3$ mol$^{-1}$, $d$: m, $D$: m$^{-2}$, $C_{OH}$: ppm, $n$: –, $\eta_0$: Pa s, $B$: –. {#tbl:growth-parameters}

\cleardoublepage

# Results {#sec:results}

## Simulation Snapshots: Slabs and Plumes {#sec:simulation-snapshots}

Figures \ref{fig:slab-comp-set2} and \ref{fig:plume-comp-set2} illustrate how thermodynamics and phase transition kinetics collude to complicate seismic structures in a convecting mantle. These “snapshots”, taken after 100 Myr of evolution, provide context for the quantitative analysis presented below.

In slab models (Figure \ref{fig:slab-comp-set2}), sluggish kinetics allow extensive metastable olivine to persist, producing displaced PTZs with faded $\hat{T}$, $\hat{\rho}$, and $V_p$ gradients. As kinetics accelerate, the olivine $\Leftrightarrow$ wadsleyite transition progresses more rapidly within the slab interior, sharpening the PTZ and generating narrower, higher-contrast $V_p$ gradients closer to the nominal (equilibrium) transition depth. Faster kinetics also promotes continuous slab descent through the transition zone, while sluggish kinetics favor stagnation and ponding.

Plume models (Figure \ref{fig:plume-comp-set2}) exhibit complementary behavior. Sluggish kinetics delay the wadsleyite $\Leftrightarrow$ olivine reaction within hot upwellings, yielding broad PTZs displaced upward by tens of kilometers. Moderately fast kinetics concentrate the phase transition into a narrower vertical interval, strengthening $\hat{\rho}$ anomalies and $V_p$ gradients. Under fast kinetic regimes, plumes display sharp reaction fronts confined to thin PTZ interfaces, with the largest density anomalies and most strongly focused seismic contrasts. Rapid conversion also enhances the vertical velocity of plume heads and penetration through the PTZ.

These visualizations demonstrate how the kinetics of the olivine $\Leftrightarrow$ wadsleyite transition govern whether slabs and plumes generate diffuse, low-amplitude anomalies or sharp, high-contrast seismic signatures. They also demonstrate that phase transition kinetics---in addition to standard thermodynamic considerations (i.e., Clapeyron slopes)---play a crucial role in developing PTZ topography. Larger sets of visualized simulation outputs for the models in Figures \ref{fig:slab-comp-set2}–\ref{fig:plume-comp-set2} are shown in Appendix \ref{sec:simulation-snapshots-continued}.

![Slab model with sluggish (a–c: slab-lnk18-Ha300-d5mm-060ppm), moderately-sluggish (d–f: slab-lnk18-Ha274-d5mm-060ppm), and fast (g–i: slab-lnk18-Ha274-D1e12-060ppm) kinetics after 100 Ma evolution. Panels show dynamic temperature $\hat{T}$ (left column), dynamic density $\hat{\rho}$ (middle column), and pressure-wave velocity $V_p$ (right column).](../figs/simulation/compositions/slab-slow-default-fast-set2-composition-0100.png){#fig:slab-comp-set2}

![Plume model with sluggish (a–c: plume-lnk18-Ha300-d5mm-060ppm), moderately-sluggish (d–f: plume-lnk18-Ha274-d5mm-060ppm), and fast (g–i: plume-lnk18-Ha274-D1e12-060ppm) kinetics after 100 Ma evolution. Panels show dynamic temperature $\hat{T}$ (left column), dynamic density $\hat{\rho}$ (middle column), and pressure-wave velocity $V_p$ (right column).](../figs/simulation/compositions/plume-slow-default-fast-set2-composition-0100.png){#fig:plume-comp-set2}

\cleardoublepage

## Phase Transition Zone: Displacement and Width {#sec:ptz-displacement-width}

Figures \ref{fig:ptz-plumes}–\ref{fig:ptz-slabs} summarize the quantitative relationships among PTZ displacement, effective PTZ width, and maximum phase transition rate $\dot{X}_{\mathrm{max}}$ in plume and slab models after 100 Ma of evolution (Table \ref{tbl:centerline-profile-results}). PTZ width and displacement were evaluated from phase fraction field $X$ along a vertical profile through the center of the model domain, where width was defined as the difference between the depths at $X$ = 0.9 and $X$ = 0.1, and displacement was defined as the offset between the nominal equilibrium reaction depth and the depth at $X$ = 0.9. The maximum phase transition rate, $\dot{X}_{\mathrm{max}}$, was evaluated from the phase transition rate field $\dot{X}$ along the same vertical profile. The results highlight contrasting responses of the PTZ in plume and slab models to variations in phase transition kinetics.

Plume models (Figure \ref{fig:ptz-plumes}) show well-defined monotonic trends. PTZ width increases linearly with upward displacement, while both width and displacement decrease smoothly (as power-laws) with increasing $\dot{X}_{\mathrm{max}}$. These consistent relationships demonstrate simple scaling between PTZ structure and phase transition kinetics for plumes.

![Quantitative relationships between PTZ structure and phase transition rate in plume models after 100 Ma. Plume PTZs widen linearly with upward displacement (a), and both width (b) and displacement (c) decrease monotonically with increasing phase transition rate.](../figs/ptz-plumes.png){#fig:ptz-plumes}

In contrast to plume models, slab models (Figure \ref{fig:ptz-slabs}) show nonlinear and non-monotonic behavior. PTZ width first broadens with downward displacement but narrows again once displacement exceeds ~-50 km, forming a quadratic trend. The quadratic trend in Figure \ref{fig:ptz-slabs}a points to a threshold effect, where strong thermodynamic driving forces accelerate reactions and sharpen the PTZ at displacements $\leq$ -50 km. This threshold effect is also evident in Figures \ref{fig:ptz-slabs}b–\ref{fig:ptz-slabs}c, where PTZ width first increases and then decreases with increasing $\dot{X}_{\mathrm{max}}$ as a complex power law (quadratic in log-log space), and displacement versus $\dot{X}_{\mathrm{max}}$ shows a weaker power-law correlation than plume models.

![Quantitative relationships between PTZ structure and phase transition rate in slab models after 100 Ma. Slab PTZs broaden with moderate downward displacement but narrow again at displacements $\leq$ -50 km (a). Width increases, then decreases with phase transition rate, but not as a simple power law (b), while displacement versus phase transition rate shows a weaker correlation than plume models (c). Outliers (in red) deviate from the quadratic trend in (a). Solid lines exclude outliers, while dashed lines were fit to the entire dataset.](../figs/ptz-slabs.png){#fig:ptz-slabs}

In summary, plume PTZs follow linear and power-law scalings that are smooth and predictable, whereas slab PTZs exhibit quadratic and curved relationships shaped by threshold behavior. These contrasting trends underscore the different ways phase transition kinetics govern PTZ structure in upwellings versus subducting slabs.

\cleardoublepage

# Discussion {#sec:discussion}

Our numerical simulations show that growth-controlled phase-transition kinetics, when coupled to realistic thermodynamics and a compressible treatment of mantle flow, exert a first-order control on the geometry and seismic expression of the 410 km phase transition. Our results demonstrate that plume and slab flow regimes respond in systematically different ways: plumes produce monotonic, smoothly varying PTZ scalings with kinetics, whereas slabs show thresholded, non-monotonic behavior (Figures \ref{fig:ptz-plumes}–\ref{fig:ptz-slabs}; Table \ref{tbl:centerline-profile-results}). The implications of such contrasting relationships in slabs versus plumes are discussed below.

## Uncertainties and Model Limitations {#sec:uncertainties-and-model-limitations}

The kinetic parameters employed here span several orders of magnitude (Table \ref{tbl:growth-parameters}), reflecting the substantial uncertainties in experimental constraints on olivine $\Leftrightarrow$ wadsleyite transition rates. Laboratory studies yield activation enthalpies $H^{\ast}$ ranging from ~200–500 kJ/mol and kinetic prefactors $A$ varying by 6–8 orders of magnitude depending on water content, grain size, and deformation mechanism [@rubie1994; @kubo2004; @hosoya2005]. Moreover, our choice of grain-boundary versus dislocation-controlled growth mechanisms significantly impacts the nucleation site factor $S$ and thus the overall transition rate [@mohiuddin2018]. These parameter uncertainties are the principal source of quantitative uncertainty in PTZ widths, displacements and associated seismic contrasts in our numerical simulations.

We assumed rapid nucleation saturation and growth-limited kinetics (Equations \ref{eq:volume-fraction}–\ref{eq:phase-transition-rate}), which is appropriate for many experimental conditions [@cahn1956; @hosoya2005; @faccenda2017]. However, experiments and in-situ studies show that nucleation can be important---especially at low temperatures and in coarse-grained or dry lithologies---and can slow the net phase transition relative to a pure growth-limited formulation [@rubie1994; @kubo2004; @perrillat2016]. Recent in-situ X-ray and acoustic experiments and microstructure studies [@ohuchi2022; @ledoux2023] further document complex nucleation/growth microstructures (including nanocrystalline, incoherent products) that can limit effective phase transition rates under some PT–deformation paths. These results indicate that our assumption of saturated nucleation is an approximation that will tend to decrease metastability and sharpen PTZs in cold slabs.

Assuming pure Mg-rich end-members neglects Fe-partitioning and minor-element effects that can shift equilibrium depths by ~10–20 km and alter kinetics [@katsura2004; @perrillat2016]. Likewise, our assumption of a relatively dry mantle neglects any potential sharpening and/or shifting of the 410 km seismic discontinuity due to variable water contents in olivine and wadsleyite [@chen2002]. Our simple temperature-dependent viscosity also omits grain-size evolution, stress-dependent rheologies and anisotropic deformation, all of which can modify strain localization and hence the thermodynamic/kinetic feedbacks controlling PTZ structure [@karato2001]. These omissions should be addressed in future higher-fidelity models---for now they imply that the quantitative thresholds reported above are model-dependent and should be interpreted as first-order effects.

## Implications for Subduction Dynamics {#sec:implications-for-subduction-dynamics}

The quadratic relationship between PTZ width and displacement in slabs (Figure \ref{fig:ptz-slabs}a) suggests a critical threshold near -50 km displacement where thermodynamic forces overcome kinetic barriers. Above this threshold---where kinetics are sufficiently sluggish---metastable olivine persists extensively, consistent with deep earthquake observations attributed to transformational faulting in metastable wedges [@kirby1996; @green1995; @ishii2021; @sindhusuta2025]. Our results indicate that $\dot{X}_{\mathrm{max}}$ < 0.2 Ma$^{-1}$ produce extensive downward PTZ displacement that promote slab stagnation at the 410. In contrast, $\dot{X}_{\mathrm{max}}$ > 1 Ma$^{-1}$ show small PTZ displacements that enable continuous penetration into the lower mantle. Thus, in our models the transition between a metastable, ponding regime and effective penetration occurs over roughly one to two orders of magnitude in $\dot{X}_{\mathrm{max}}$ (Table \ref{tbl:centerline-profile-results}).

These mechanistic thresholds have observational relevance. Regions where seismic tomography and receiver functions show slab flattening or stagnation above the transition zone [e.g., @fukao2013] are compatible with low effective phase transition rates and extensive metastability in our models. At the same time, transformational faulting and deep earthquakes in metastable wedges remain consistent with substantial metastable volumes and sharp local stress concentrations produced during delayed phase transitions [@ohuchi2022]. Furthermore, the non-monotonic width-rate relationship in slabs (Figure \ref{fig:ptz-slabs}b) implies that intermediate kinetics (0.05–0.15  Ma$^{-1}$) generate the broadest PTZs, potentially explaining the variable thickness of the 410 discontinuity beneath different subduction zones [@lee2014; @vanstiphout2019; @jiang2015; @shen2020; @han2021]. Our work therefore demonstrates that kinetics provides an important control on the expression of 410 seismic discontinuity, particularly where cold slabs create conditions far from equilibrium.

## Seismic Structure and Discontinuity Sharpness {#sec:seismic-structure-and-discontinuity-sharpness}

The correlation between $\dot{X}_{\mathrm{max}}$ and PTZ structure has direct implications for interpreting seismic discontinuities. Sharp interfaces are most readily detected with SS precursors and receiver-function stacks. Global SS/PP precursor studies emphasize preferential sensitivity to very thin PTZs [a few km, @shearer2000; @chambers2005; @deuss2009], but actual resolution depends on frequency content, signal-to-noise, and processing strategy. High-frequency and regional receiver-function approaches demonstrate that structures as thin as ~5 km—and occasionally thinner—can be resolved under favorable conditions [@helffrich1996; @wei2017; @dokht2016; @frazer2023].

In our models, thin and easily detectable discontinuities correspond to $\dot{X}_{\mathrm{max}}$ > 4 Ma$^{-1}$ in plumes and > 1 Ma$^{-1}$ in slabs (Figures \ref{fig:ptz-plumes}–\ref{fig:ptz-slabs}; Table \ref{tbl:centerline-profile-results}). This implies that regions with sharp, high-amplitude 410 signals may reflect either fast kinetics that maintain reactions close to equilibrium or large effective displacements across the PTZ. Observationally, sharp 410s are reported beneath several plume regions, but sharpness is not uniquely diagnostic of kinetics: composition, anisotropy, and imaging method also play important roles [@lawrence2008; @deuss2009].

Discontinuities observed near 500–600 km depth provide an additional cautionary example. These features have been attributed to akimotoite formation or other mid-transition-zone reactions [@cottaar2016], and alternative explanations include $\beta \Leftrightarrow \gamma$ olivine transitions or compositional heterogeneity [@saikia2008; @deuss2001; @tauzin2017]. Because the olivine $\alpha \Leftrightarrow \beta$ transition occurs much shallower (~410 km), we do not interpret 500–600 km signals as delayed wadsleyite formation. Conversely, broad or absent discontinuities within the upper transition zone may indicate sluggish kinetics, in addition to thermal or compositional influences.

The systematic upward displacement and broadening of the 410 produced by sluggish kinetics in hot upwellings (Figure \ref{fig:ptz-plumes}) offers a kinetic complement to purely thermal explanations for reduced-amplitude or displaced discontinuities beneath hotspots. Global SS and receiver-function compilations show reduced amplitude and complex expressions of the 410 beneath many hotspots [@deuss2009], and high-resolution receiver-function and precursor studies have reported topography and thickness anomalies on scales of tens of kilometres [@chambers2005; @agius2017; @jenkins2016; @glasgow2024]. Accounting for kinetically controlled broadening and displacement (10–30 km in many of our plume cases) reduces the remaining discrepancy between thermal predictions and the observed total uplift or weakening.

Finally, observed lateral variations in apparent 410 thickness---ranging from ~5–30 km in some Pacific regions [@alex2004; g@schmerr2007]---can be produced either by thermal/compositional heterogeneity or by spatial changes in effective kinetics ($\dot{X}_{\mathrm{max}}$) or both. Because our model suite spans ~2–3 orders of magnitude in $\dot{X}_{\mathrm{max}}$ (Table \ref{tbl:centerline-profile-results}), the observed variability is broadly consistent with realistic variations in both kinetics and background mantle properties.

## Constraining Kinetic Parameters from Seismic Observations {#sec:constraining-kinetic-parameters-from-seismic-observations}

The distinct and relatively simple scaling laws we find for plumes (linear width–displacement and power-law width–$\dot{X}_{\mathrm{max}}$ scalings; Figure \ref{fig:ptz-plumes}) provide a potential pathway to infer effective phase transition rates from seismic observables, provided independent constraints on thermal structure are available (e.g., from tomography or surface heat-flow). In slabs the relation is more complex and thresholded (Figure \ref{fig:ptz-slabs}), but the threshold behaviour itself---the depth and abruptness of the catastrophic conversion---is a diagnostic that can distinguish slow vs. fast kinetic regimes. Operationally, two complementary approaches are promising: thickness inversion coupled to independent thermometry, and targeted regional tests.

High-resolution receiver-function or SS precursor maps [@lawrence2008; @schmerr2007; @chambers2005; @houser2010] can produce spatial maps of apparent PTZ thickness and displacement. Combining these with tomographic temperature estimates and forward mapping from our model scalings (Figures \ref{fig:ptz-plumes}–\ref{fig:ptz-slabs}) yields order-of-magnitude bounds on $\dot{X}_{\mathrm{max}}$ and discriminates growth mechanisms [grain-boundary vs. dislocation control, @hosoya2005]. Mineral-physics experiments that better quantify nucleation vs. growth rates, the role of water and Fe content, and microstructural inheritance [@perrillat2013; @perrillat2016; @ohuchi2022; @ledoux2023] are essential to shrink model uncertainty and make seismic inversions for kinetics quantitative. Combining such experiments with the forward model scalings we present here is a roadmap to constrain effective mantle kinetics from seismic data.

Moreover, regions where thermal structure is relatively constrained but seismic expression varies laterally (for example, across slabs with different ages or hydration states, or across hotspot swells vs. ambient mantle) are high-value tests. Systematic surveys that compare PTZ thickness/topography in tectonically similar thermal settings but different deformation histories (or water contents) can isolate kinetic effects from purely thermal/compositional signals [e.g., @agius2017; @vanstiphout2019; @schmandt2012; @perrillat2022].

\cleardoublepage

# Conclusions {#sec:conclusions}

Seismic structure in Earth’s upper mantle reflects a balance between equilibrium thermodynamics and phase transition kinetics. This study set out to quantify how these factors contribute to the expression of the 410 km seismic discontinuity (olivine $\Leftrightarrow$ wadsleyite transition) in plumes and slabs.

By coupling equilibrium thermodynamics to a growth-controlled kinetic formulation within compressible mantle flow models, we systematically explored the sensitivity of PTZ structure to kinetic rate parameters across a broad experimental uncertainty range. Our plume simulations reveal smooth, power-law scaling: sluggish kinetics broaden and uplift the 410, while faster kinetics sharpen it to seismically detectable thicknesses. Slab simulations, in contrast, display thresholded dynamics: moderately slow kinetics yield metastable wedges, broad PTZs, and slab stagnation at the 410 km PTZ, whereas both ultra-slow ($\dot{X}_{\mathrm{max}}$ < 0.02 Ma$^{-1}$) and ultra-fast ($\dot{X}_{\mathrm{max}}$ > 1 Ma$^{-1}$) kinetics permit sharp seismic contrasts.

These results show that phase transition kinetics in compressibile mantle flow can reproduce the observed diversity of 410 km discontinuity sharpness and topography across the globe. They further indicate that kinetic thresholds of $\dot{X}_{\mathrm{max}}$ < 1 Ma$^{-1}$ separate stagnant, metastable regimes from penetrating, equilibrated ones. Such thresholds link mantle dynamics to effective phase transition rates, suggesting that regional seismic observations can be inverted---with supporting thermal constraints---to constrain micro-scale kinetic parameters.

In summmary, our results demonstrate the utility of the 410 km discontinuity as a potential seismological probe of kinetics in Earth's upper mantle. Integrating seismic imaging, mineral physics, and forward geodynamic models offers a path toward quantifying the role of kinetics in mantle dynamics and in shaping the seismic expression of the MTZ.

\cleardoublepage

# Acknowledgements {.unnumbered}

This work was funded by the UKRI NERC Large Grant no. NE/V018477/1 awarded to John Wheeler at the University of Liverpool. All computations were undertaken on Barkla2, part of the High Performance Computing facilities at the University of Liverpool, who graciously provided expert support. We thank the Computational Infrastructure for Geodynamics ([geodynamics.org](geodynamics.org)) which is funded by the National Science Foundation under award EAR-0949446 and EAR-1550901 for supporting the development of ASPECT.

# Data Availability {.unnumbered}

All data, code, and relevant information for reproducing this work can be found at [https://github.com/buchanankerswell/kerswell_et_al_dynp](https://github.com/buchanankerswell/kerswell_et_al_dynp), and at ..., the official Open Science Framework data repository. All code is MIT Licensed and free for use and distribution (see license details). ASPECT version 3.0.0, [@aspect-doi-v3.0.0; @aspectmanual; @heister2017; @kronbichler2012; @gassmoller2018; @clevenger2021; @fraters2019; @fraters2020] used in these computations is freely available under the GPL v2.0 or later license through its software landing page [https://geodynamics.org/resources/aspect](https://geodynamics.org/resources/aspect) or [https://aspect.geodynamics.org](https://aspect.geodynamics.org) and is being actively developed on GitHub and can be accessed via [https://github.com/geodynamics/aspect](https://github.com/geodynamics/aspect).

\cleardoublepage

# References {.unnumbered}

::: {#refs}
:::

\cleardoublepage

# Appendix {.unnumbered}

## Stress, Pressure, and The Momentum Equation {.unnumbered #sec:momentum-derivation}

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

where $\dot{\epsilon}^{\prime} = \frac{1}{2} \left(\nabla \vec{u} + \left(\nabla \vec{u} \right)^\intercal \right) - \frac{1}{3} \left(\nabla \cdot \vec{u} \right) I$ is the deviatoric strain rate tensor. In full component form the deviatoric stress tensor is:

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

\cleardoublepage

## Simulation Snapshots: Slabs and Plumes Continued {.unnumbered #sec:simulation-snapshots-continued}

![Slab model with sluggish kinetics (slab-lnk18-Ha300-d5mm-060ppm) after 100 Ma evolution. Panels show dynamic temperature $\hat{T}$ (a), dynamic pressure $\hat{P}$ (b), dynamic density $\hat{\rho}$ (c), thermodynamic term $\left(1 - \left[\Delta G/R\,T\right]\right)$ (d), growth rate $\dot{x}$ (e), phase transition rate $\dot{X}$ (f), volume fraction of wadsleyite $X$ (g), pressure-wave velocity $V_p$ (h), and shear-wave velocity $V_s$ (i).](../figs/simulation/compositions/slab-lnk18-Ha300-d5mm-060ppm-full-set-composition-0100.png){#fig:slab-comp-slow}

![Slab model with moderately-sluggish kinetics (slab-lnk18-Ha274-d5mm-060ppm) after 100 Ma evolution. Panels show dynamic temperature $\hat{T}$ (a), dynamic pressure $\hat{P}$ (b), dynamic density $\hat{\rho}$ (c), thermodynamic term $\left(1 - \left[\Delta G/R\,T\right]\right)$ (d), growth rate $\dot{x}$ (e), phase transition rate $\dot{X}$ (f), volume fraction of wadsleyite $X$ (g), pressure-wave velocity $V_p$ (h), and shear-wave velocity $V_s$ (i).](../figs/simulation/compositions/slab-lnk18-Ha274-d5mm-060ppm-full-set-composition-0100.png){#fig:slab-comp-moderately-slow}

![Slab model with fast kinetics (slab-lnk18-Ha274-D1e12-060ppm) after 100 Ma evolution. Panels show dynamic temperature $\hat{T}$ (a), dynamic pressure $\hat{P}$ (b), dynamic density $\hat{\rho}$ (c), thermodynamic term $\left(1 - \left[\Delta G/R\,T\right]\right)$ (d), growth rate $\dot{x}$ (e), phase transition rate $\dot{X}$ (f), volume fraction of wadsleyite $X$ (g), pressure-wave velocity $V_p$ (h), and shear-wave velocity $V_s$ (i).](../figs/simulation/compositions/slab-lnk18-Ha274-D1e12-060ppm-full-set-composition-0100.png){#fig:slab-comp-fast}

![Plume model with sluggish kinetics (plume-lnk18-Ha300-d5mm-060ppm) after 100 Ma evolution. Panels show dynamic temperature $\hat{T}$ (a), dynamic pressure $\hat{P}$ (b), dynamic density $\hat{\rho}$ (c), thermodynamic term $\left(1 - \left[\Delta G/R\,T\right]\right)$ (d), growth rate $\dot{x}$ (e), phase transition rate $\dot{X}$ (f), volume fraction of olivine $X$ (g), pressure-wave velocity $V_p$ (h), and shear-wave velocity $V_s$ (i).](../figs/simulation/compositions/plume-lnk18-Ha300-d5mm-060ppm-full-set-composition-0100.png){#fig:plume-comp-slow}

![Plume model with moderately-sluggish kinetics (plume-lnk18-Ha274-d5mm-060ppm) after 100 Ma evolution. Panels show dynamic temperature $\hat{T}$ (a), dynamic pressure $\hat{P}$ (b), dynamic density $\hat{\rho}$ (c), thermodynamic term $\left(1 - \left[\Delta G/R\,T\right]\right)$ (d), growth rate $\dot{x}$ (e), phase transition rate $\dot{X}$ (f), volume fraction of olivine $X$ (g), pressure-wave velocity $V_p$ (h), and shear-wave velocity $V_s$ (i).](../figs/simulation/compositions/plume-lnk18-Ha274-d5mm-060ppm-full-set-composition-0100.png){#fig:plume-comp-moderately-slow}

![Plume model with fast kinetics (plume-lnk18-Ha274-D1e12-060ppm) after 100 Ma evolution. Panels show dynamic temperature $\hat{T}$ (a), dynamic pressure $\hat{P}$ (b), dynamic density $\hat{\rho}$ (c), thermodynamic term $\left(1 - \left[\Delta G/R\,T\right]\right)$ (d), growth rate $\dot{x}$ (e), phase transition rate $\dot{X}$ (f), volume fraction of olivine $X$ (g), pressure-wave velocity $V_p$ (h), and shear-wave velocity $V_s$ (i).](../figs/simulation/compositions/plume-lnk18-Ha274-D1e12-060ppm-full-set-composition-0100.png){#fig:plume-comp-fast}

\cleardoublepage

## Phase Transition Zone: Displacement and Width Continued {.unnumbered #sec:ptz-displacement-width-continued}

| Model ID | PTZ Displacement | PTZ Width | $\dot{X}_{\mathrm{max}}$ |
|:--|--:|--:|--:|
| plume-lnk21-Ha274-d5mm-060ppm | 75.86 | -73.69 | 0.231 |
| plume-lnk21-Ha274-d5mm-060ppm | 75.86 | -73.69 | 0.231 |
| plume-lnk20-Ha274-d5mm-060ppm | 41.25 | -46.17 | 0.424 |
| plume-lnk20-Ha274-d5mm-060ppm | 41.25 | -46.17 | 0.424 |
| plume-lnk18-Ha300-d5mm-060ppm | 30.94 | -38.01 | 0.538 |
| plume-lnk18-Ha300-d5mm-060ppm | 30.94 | -38.01 | 0.538 |
| plume-lnk19-Ha274-d5mm-060ppm | 18.07 | -27.67 | 0.751 |
| plume-lnk19-Ha274-d5mm-060ppm | 18.07 | -27.67 | 0.751 |
| plume-lnk18-Ha274-d1cm-060ppm | 12.95 | -23.73 | 0.889 |
| plume-lnk18-Ha274-d1cm-060ppm | 12.95 | -23.73 | 0.889 |
| plume-lnk18-Ha274-d5mm-060ppm | 4.12 | -16.50 | 1.292 |
| plume-lnk18-Ha274-d5mm-060ppm | 4.12 | -16.50 | 1.292 |
| plume-lnk18-Ha274-D1e06-060ppm | 0.65 | -14.05 | 1.604 |
| plume-lnk18-Ha274-D1e06-060ppm | 0.65 | -14.05 | 1.604 |
| plume-lnk18-Ha274-d5mm-080ppm | -3.09 | -11.16 | 2.095 |
| plume-lnk18-Ha274-d5mm-080ppm | -3.09 | -11.16 | 2.095 |
| plume-lnk17-Ha274-d5mm-060ppm | -4.12 | -10.25 | 2.188 |
| plume-lnk17-Ha274-d5mm-060ppm | -4.12 | -10.25 | 2.188 |
| plume-lnk18-Ha274-d1mm-060ppm | -7.22 | -7.83 | 2.973 |
| plume-lnk18-Ha274-d1mm-060ppm | -7.22 | -7.83 | 2.973 |
| plume-lnk18-Ha274-d5mm-100ppm | -7.22 | -7.84 | 3.017 |
| plume-lnk18-Ha274-d5mm-100ppm | -7.22 | -7.84 | 3.017 |
| plume-lnk16-Ha274-d5mm-060ppm | -8.46 | -6.83 | 3.623 |
| plume-lnk16-Ha274-d5mm-060ppm | -8.46 | -6.83 | 3.623 |
| plume-lnk18-Ha274-d5mm-120ppm | -9.39 | -6.08 | 4.006 |
| plume-lnk18-Ha274-d5mm-120ppm | -9.39 | -6.08 | 4.006 |
| plume-lnk18-Ha274-D1e08-060ppm | -11.17 | -4.52 | 5.146 |
| plume-lnk18-Ha274-D1e08-060ppm | -11.17 | -4.52 | 5.146 |
| plume-lnk18-Ha274-d5mm-140ppm | -11.19 | -4.51 | 5.176 |
| plume-lnk18-Ha274-d5mm-140ppm | -11.19 | -4.51 | 5.176 |
| plume-lnk15-Ha274-d5mm-060ppm | -11.34 | -4.12 | 5.859 |
| plume-lnk15-Ha274-d5mm-060ppm | -11.34 | -4.12 | 5.859 |
| plume-lnk18-Ha274-D1e10-060ppm | -14.44 | -2.06 | 14.335 |
| plume-lnk18-Ha274-D1e10-060ppm | -14.44 | -2.06 | 14.335 |
| plume-lnk18-Ha274-D1e12-060ppm | -15.69 | -2.87 | 31.381 |
| plume-lnk18-Ha274-D1e12-060ppm | -15.69 | -2.87 | 31.381 |
| slab-lnk21-Ha274-d5mm-060ppm | -109.82 | 37.55 | 0.012 |
| slab-lnk21-Ha274-d5mm-060ppm | -109.82 | 37.55 | 0.012 |
| slab-lnk20-Ha274-d5mm-060ppm | -86.41 | 32.46 | 0.019 |
| slab-lnk20-Ha274-d5mm-060ppm | -86.41 | 32.46 | 0.019 |
| slab-lnk18-Ha300-d5mm-060ppm | -81.30 | 29.59 | 0.020 |
| slab-lnk18-Ha300-d5mm-060ppm | -81.30 | 29.59 | 0.020 |
| slab-lnk19-Ha274-d5mm-060ppm | -66.26 | 35.32 | 0.029 |
| slab-lnk19-Ha274-d5mm-060ppm | -66.26 | 35.32 | 0.029 |
| slab-lnk18-Ha274-d1cm-060ppm | -61.03 | 37.31 | 0.033 |
| slab-lnk18-Ha274-d1cm-060ppm | -61.03 | 37.31 | 0.033 |
| slab-lnk18-Ha274-d5mm-060ppm | -51.10 | 42.85 | 0.045 |
| slab-lnk18-Ha274-d5mm-060ppm | -51.10 | 42.85 | 0.045 |
| slab-lnk18-Ha274-D1e06-060ppm | -46.44 | 44.38 | 0.055 |
| slab-lnk18-Ha274-D1e06-060ppm | -46.44 | 44.38 | 0.055 |
| slab-lnk18-Ha274-d5mm-080ppm | -41.25 | 44.34 | 0.072 |
| slab-lnk18-Ha274-d5mm-080ppm | -41.25 | 44.34 | 0.072 |
| slab-lnk17-Ha274-d5mm-060ppm | -40.57 | 44.25 | 0.075 |
| slab-lnk17-Ha274-d5mm-060ppm | -40.57 | 44.25 | 0.075 |
| slab-lnk18-Ha274-d1mm-060ppm | -34.69 | 42.11 | 0.108 |
| slab-lnk18-Ha274-d1mm-060ppm | -34.69 | 42.11 | 0.108 |
| slab-lnk18-Ha274-d5mm-100ppm | -34.45 | 41.99 | 0.109 |
| slab-lnk18-Ha274-d5mm-100ppm | -34.45 | 41.99 | 0.109 |
| slab-lnk16-Ha274-d5mm-060ppm | -30.94 | 40.43 | 0.137 |
| slab-lnk16-Ha274-d5mm-060ppm | -30.94 | 40.43 | 0.137 |
| slab-lnk18-Ha274-d5mm-120ppm | -28.98 | 39.41 | 0.159 |
| slab-lnk18-Ha274-d5mm-120ppm | -28.98 | 39.41 | 0.159 |
| slab-lnk18-Ha274-D1e08-060ppm | -23.80 | 36.18 | 0.226 |
| slab-lnk18-Ha274-D1e08-060ppm | -23.80 | 36.18 | 0.226 |
| slab-lnk18-Ha274-d5mm-140ppm | -23.73 | 36.19 | 0.227 |
| slab-lnk18-Ha274-d5mm-140ppm | -23.73 | 36.19 | 0.227 |
| slab-lnk15-Ha274-d5mm-060ppm | -20.62 | 34.03 | 0.281 |
| slab-lnk15-Ha274-d5mm-060ppm | -20.62 | 34.03 | 0.281 |
| slab-lnk18-Ha274-D1e10-060ppm | 3.09 | 14.34 | 1.144 |
| slab-lnk18-Ha274-D1e10-060ppm | 3.09 | 14.34 | 1.144 |
| slab-lnk18-Ha274-D1e12-060ppm | 14.44 | 4.12 | 4.388 |
| slab-lnk18-Ha274-D1e12-060ppm | 14.44 | 4.12 | 4.388 |

Table: PTZ displacement, width, and maximum phase transition rate $\dot{X}_{\mathrm{max}}$ evaluated in plume and slab simulations after 100 Ma of evolution. Units are PTZ displacement: km, PTZ width: km, $\dot{X}_{\mathrm{max}}$: Ma$^{-1}$. {#tbl:centerline-profile-results}
