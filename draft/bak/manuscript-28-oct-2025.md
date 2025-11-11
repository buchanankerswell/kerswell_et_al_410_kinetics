---
title: "Displaced and Faded 410s"
subtitle: "How micro-scale kinetics complicate mantle-scale seismic structures and flow dynamics"
author: ["Kerswell B.", "Wheeler J.", "Gassmöller R."]
date: "28 October 2025"
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
header-includes:
  - \makeatletter
  - \newcounter{none}
  - \renewcommand{\thenone}{}
  - \makeatother
---

# Abstract {.unnumbered #sec:abstract}

The seismic expression of Earth's 410 km discontinuity (the "410") varies substantially across different tectonic settings, from sharp, high-amplitude interfaces to broad, diffuse transitions---patterns that cannot be explained by equilibrium thermodynamics alone. Laboratory experiments demonstrate that the olivine $\Leftrightarrow$ wadsleyite phase transition responsible for the 410 is strongly rate-limited, yet quantitative links between micro-scale kinetics and mantle-scale seismic structures and flow dynamics remain poorly understood. Here we systematically investigate these relationships by coupling an interface-controlled growth model to compressible simulations of mantle plumes and subducting slabs using ASPECT. We vary kinetic parameters across seven orders of magnitude and quantify the resulting 410 displacements and widths. Our results reveal a fundamental asymmetry between hot and cold mantle environments. In plumes, high temperatures suppress olivine metastability, producing consistently sharp 410s (2–5 km wide) that remain insensitive to kinetic variations across the explored parameter space. In slabs, by contrast, kinetics exert first-order control on 410 structure and flow dynamics. We identify three distinct kinetic regimes in these cold environments: (1) quasi-equilibrium behavior at high phase transition rates ($\dot{X}$ > 10$^1$ Ma$^{-1}$) enabling continuous slab penetration with narrow, positively displaced 410s; (2) intermediate phase transition rates (10$^{-1.5}$ < $\dot{X}$ < 10$^1$ Ma$^{-1}$) generating progressively broader, deeper 410s and metastable olivine wedges that resist slab descent without preventing it; and (3) ultra-sluggish phase transition rates ($\dot{X}$ < 10$^{-1.5}$ Ma$^{-1}$) causing complete slab stagnation with re-sharpened but deeply displaced 410s (> 100 km). These findings demonstrate that phase transition rates strongly influence 410 structure in subduction zones and establish the 410 as a potential seismological constraint on kinetic processes operating in Earth's upper mantle, particularly in cold environments where disequilibrium effects are amplified.

\cleardoublepage

# Definition of Symbols {.unnumbered #sec:symbols}

|Parameter|Symbol|Unit|Equations|
|:--------------|:------|:----|:----|
|Activation enthalpy|$H^{\ast}$|J mol$^{-1}$|\ref{eq:growth-rate}, \ref{eq:phase-transition-rate}|
|Activation volume|$V^{\ast}$|m$^3$ mol$^{-1}$|\ref{eq:growth-rate}, \ref{eq:phase-transition-rate}|
|Compressibility (reference)|$\bar{\beta}$|Pa$^{-1}$|\ref{eq:density-ala}|
|Density|$\rho$|kg m$^{-3}$|\ref{eq:navier-stokes-no-inertia}–\ref{eq:continuity-expanded}, \ref{eq:density-ala-expansion}–\ref{eq:density-ala}|
|Density (reference)|$\bar{\rho}$|kg m$^{-3}$|\ref{eq:adiabatic-pressure}–\ref{eq:density-ala}|
|Density (dynamic)|$\hat{\rho}$|kg m$^{-3}$|-|
|Deviatoric stress tensor|$\sigma^{\prime}$|Pa|\ref{eq:navier-stokes-no-inertia}, \ref{eq:energy}|
|Deviatoric strain rate tensor|$\dot{\epsilon}^{\prime}$|s$^{-1}$|\ref{eq:energy}|
|Gas constant|$R$|J K$^{-1}$ mol$^{-1}$|\ref{eq:growth-rate}, \ref{eq:phase-transition-rate}|
|Grain size|$d$|1 m|\ref{eq:growth-rate}, \ref{eq:phase-transition-rate}|
|Gravitational acceleration|$g$|m s$^{-2}$|\ref{eq:navier-stokes-no-inertia}, \ref{eq:adiabatic-temperature}–\ref{eq:adiabatic-pressure}|
|Growth rate|$\dot{x}$|m s$^{-1}$|\ref{eq:volume-fraction}, \ref{eq:growth-rate}|
|Kinetic prefactor|$A$|m s$^{-1}$ K$^{-1}$ ppm$_\mathrm{OH}^{-n}$|\ref{eq:growth-rate}|
|Latent heat|$Q_L$|J kg$^{-1}$|\ref{eq:energy}|
|Molar entropy|$\bar{S}$|J mol$^{-1}$ K$^{-1}$|\ref{eq:gibbs}|
|Molar Gibbs free energy|$\bar{G}$|J mol$^{-1}$|\ref{eq:gibbs}|
|Molar volume|$\bar{V}$|m$^{3}$ mol$^{-1}$|\ref{eq:gibbs}|
|Nucleation site factor|$N$|m$^{-1}$|\ref{eq:volume-fraction}, \ref{eq:phase-transition-rate-short}|
|Phase transition rate|$\frac{dX}{dt}$, $\frac{\partial X}{\partial t}$, $\dot{X}$|s$^{-1}$|\ref{eq:phase-transition-rate-short}–\ref{eq:composition}|
|Pressure|$P$|Pa|\ref{eq:navier-stokes-no-inertia}, \ref{eq:energy}, \ref{eq:growth-rate}, \ref{eq:phase-transition-rate}|
|Pressure (reference)|$\bar{P}$|K|\ref{eq:adiabatic-pressure}|
|Pressure (dynamic)|$\hat{P}$|Pa|\ref{eq:density-ala}, \ref{eq:gibbs}|
|Specific heat capacity (reference)|$\bar{C}_p$|J kg$^{-1}$ K$^{-1}$|\ref{eq:energy}, \ref{eq:adiabatic-temperature}|
|Temperature|$T$|K|\ref{eq:energy}, \ref{eq:growth-rate}, \ref{eq:phase-transition-rate}|
|Temperature (reference)|$\bar{T}$|K|\ref{eq:adiabatic-temperature}, \ref{eq:rheological-model}|
|Temperature (dynamic)|$\hat{T}$|K|\ref{eq:density-ala}, \ref{eq:gibbs}, \ref{eq:rheological-model}|
|Thermal conductivity (reference)|$\bar{k}$|W m$^{-1}$ K$^{-1}$|\ref{eq:energy}|
|Thermal expansivity (reference)|$\bar{\alpha}$|Pa$^{-1}$|\ref{eq:energy}, \ref{eq:adiabatic-temperature}, \ref{eq:density-ala}|
|Thermal viscosity exponent factor|$B$|-|\ref{eq:rheological-model}|
|Time|$t$|s|\ref{eq:continuity-compressible}–\ref{eq:continuity-expanded}, \ref{eq:volume-fraction}, \ref{eq:phase-transition-rate-short}–\ref{eq:composition}|
|Velocity|$\vec{u}$|m s$^{-1}$|\ref{eq:continuity-compressible}–\ref{eq:continuity-expanded}, \ref{eq:composition}|
|Viscosity|$\eta$|Pa s|\ref{eq:rheological-model}|
|Viscosity (reference)|$\bar{\eta}$|Pa s|\ref{eq:rheological-model}|
|Volume fraction|$X$|-|\ref{eq:volume-fraction},  \ref{eq:phase-transition-rate-short}–\ref{eq:composition}|
|Water content|$C_\mathrm{OH}$|ppm|\ref{eq:growth-rate}|
|Water content exponent|$n$|-|\ref{eq:growth-rate}|

\cleardoublepage

# Introduction {#sec:introduction}

Earth's mantle transition zone hosts two prominent seismic discontinuities near 410 and 660 km depth, attributed to polymorphic phase transitions of olivine [@ringwood1975; @katsura2004]. While these discontinuities are observed globally, their detailed seismic characteristics—depth, sharpness, amplitude, and lateral continuity—vary substantially between tectonic settings [@deuss2009; @lawrence2008; @schmerr2007; @fukao2013]. Some regions display sharp, high-amplitude reflectors consistent with abrupt mineralogical boundaries, while others exhibit broad, weakened, or laterally variable signals. Such heterogeneity cannot be explained by equilibrium thermodynamics alone, which relates discontinuity topography mainly to temperature-dependent phase boundaries defined by Clapeyron slopes. Additional physical processes---including phase transition kinetics and dynamic pressure effects---likely contribute to the observed variability [@rubie1994; @faccenda2017].

Laboratory studies provide crucial insights into the mechanisms controlling the olivine $\Leftrightarrow$ wadsleyite phase transition at 410 km depth (referred to hereafter as the nominal "410"). Mineral physics experiments consistently demonstrate that this transition is strongly rate-limited, with kinetics governed by temperature, pressure, water content, bulk chemical composition, grain size, and microstructural evolution [@rubie1994; @kubo2004; @hosoya2005; @perrillat2013; @ledoux2023]. In cold subducting slabs, sluggish reaction rates can allow metastable olivine to persist tens of kilometers below its thermodynamic stability limit, promoting slab stagnation and triggering deep earthquakes via transformational faulting [@rubie1994; @green1995; @kirby1996; @ishii2021; @ohuchi2022; @sindhusuta2025]. In hot upwellings beneath mantle plumes, slow kinetics may broaden and uplift the discontinuity, possibly explaining reduced seismic amplitudes observed beneath some hotspots [@chambers2005]. However, because published kinetic models remain poorly constrained, with parameters spanning several orders of magnitude [e.g., @hosoya2005], the effects of micro-scale kinetic processes on mantle-scale dynamics and seismic observables are ambiguous.

Bridging the gap between laboratory-derived kinetic rate laws and mantle-scale seismic observations requires numerical models that couple phase transition kinetics to realistic treatments of mantle convection. Previous modeling efforts have demonstrated qualitatively that kinetics can strongly influence mantle flow [@agrusta2017; @faccenda2017], but systematic investigations quantifying the sensitivity of 410 structure to kinetic parameters remain limited. Moreover, most prior studies employ simplified treatments of mantle compressibility that may inadequately capture feedbacks among density changes, pressure perturbations, and phase transition rates.

This study seeks to clarify these unresolved issues by implementing an interface-controlled growth model [after @hosoya2005] within the geodynamic software ASPECT. We systematically explore how thermodynamics and kinetics interact to control the structure of the 410 within a compressible treatment of mantle flow. Specifically, we address the following questions:

1. How do coupled feedbacks among thermodynamic driving forces, phase transition rates, and compressible flow shape 410 structure?
2. What are the characteristic timescales and length scales over which these feedbacks operate?
3. How do kinetic effects differ between plume and slab environments?
4. Can seismic observations of 410 structure constrain effective kinetic parameters in Earth's mantle?

To investigate these questions, we analyze a suite of numerical experiments that vary kinetic parameters across seven orders of magnitude, encompassing a wide range of water contents and grain sizes. For each experiment, we quantify 410 displacement and width, enabling direct comparisons with seismological observations. Our simulations reveal that plumes and slabs respond differently to kinetic variations and establish quantitative relationships between phase transition rates, flow dynamics, and 410 structure. More broadly, our results demonstrate that realistic treatment of phase transition kinetics is essential for accurately modeling subduction dynamics and interpreting seismic structures.

\cleardoublepage

# Methods {#sec:methods}

## Governing Equations for Compressible Mantle Flow {#sec:governing-equations}

Mantle flow is simulated using the finite-element geodynamic code ASPECT [v3.0.0, @aspect-doi-v3.0.0; @aspectmanual; @heister2017; @kronbichler2012; @gassmoller2018; @clevenger2021; @fraters2019; @fraters2020] to find the velocity $\vec{u}$, pressure $P$, and temperature $T$ fields that satisfy the following equations:

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

where $\sigma^{\prime}$ is the deviatoric stress tensor, $\rho$ is density, $g$ is gravitational acceleration, $t$ is time, $\bar{C}_p$, $\bar{k}$, $\bar{\alpha}$ are the reference specific heat capacity, thermal conductivity, and thermal expansivity, respectively (see Section \ref{sec:adiabatic-reference-conditions}), and $Q_L$ is the latent heat released or absorbed during phase transitions. Equations \ref{eq:navier-stokes-no-inertia} and \ref{eq:continuity-compressible} together describe the buoyancy-driven flow of an isotropic fluid with negligible inertia and Equation \ref{eq:energy} describes the conduction, advection, and production (or consumption) of thermal energy [@schubert2001]. Note that the pressure $P$ in this context is equal to the mean normal stress and is positive under compression: $P = - \frac{\sigma_{xx} + \sigma_{yy}}{2}$ (see Appendix \ref{sec:momentum-derivation}).

The compressible form of the continuity equation (Equation \ref{eq:continuity-compressible}) can result in dynamic feedbacks between density changes and pressures that cause numerical oscillations when the advection timestep is less than the viscous relaxation timescale [@curbelo2019]. Therefore the continuity equation is reformulated using the *projected density approximation* [PDA, @gassmoller2020] by applying the product rule to $\nabla \cdot (\rho\, \vec{u})$ and multiplying both sides of Equation \ref{eq:continuity-compressible} by $\frac{1}{\rho}$ to obtain the following expression:

\begin{equation}
  \frac{1}{\rho} \frac{\partial \rho}{\partial t} + \nabla \cdot \vec{u} + \left(\frac{1}{\rho} \nabla \rho \right) \cdot \vec{u} = 0
  \label{eq:continuity-expanded}
\end{equation}

This PDA formulation maintains the inherent coupling among density, pressure, temperature, and phase transitions while mitigating the numerical instabilities associated with compressible flow [@gassmoller2020]. These coupled feedbacks would otherwise be neglected by incompressible formulations of the continuity equation (e.g., Boussinesq). This makes the PDA approach particularly suitable for our numerical experiments, which incorporate density changes due to dynamic PT effects and the olivine $\Leftrightarrow$ wadsleyite phase transition.

## Numerical Setup {#sec:numerical-setup}

### Adiabatic Reference Conditions {#sec:adiabatic-reference-conditions}

To ensure numerical convergence, we initialized our ASPECT simulations with reasonable initial estimates of the pressure-temperature (PT) fields and material properties in Earth's upper mantle [@aspectmanual; @kronbichler2012; @heister2017]. We began by evaluating entropy changes over a PT range of 1573–1973 K and 0.001–25 GPa (Figure \ref{fig:isentrope}) using the Gibbs free energy minimization software Perple_X [v.7.0.9, @connolly2009]. We assumed a dry pyrolitic bulk composition after @green1979 and phase equilibria were evaluated in the Na$_2$O‐CaO‐FeO‐MgO‐Al$_2$O$_3$‐SiO$_2$ (NCFMAS) chemical system with thermodynamic data and solution models of @stixrude2022. Equations of state were included for solid solution phases: olivine, plagioclase, spinel, clinopyroxene, wadsleyite, ringwoodite, perovskite, ferropericlase, high‐pressure C2/c pyroxene, orthopyroxene, akimotoite, post‐perovskite, Ca‐ferrite, garnet, and Na‐Al phase.

![Entropy (a) and density (b) changes in Earth's upper mantle under thermodynamic equilibrium and hydrostatic stress conditions. Material properties were computed with Perple_X using the equations of state and thermodynamic data of @stixrude2022. The black box indicates the approximate PT range of our ASPECT simulations, while the white line indicates the isentropic adiabat used to calculate reference material properties.](../figs/PYR-material-table.png){#fig:isentrope width=70%}

We then determined the mantle adiabat by applying the Newton–Raphson algorithm to find temperatures corresponding to each pressure that maintain constant entropy (white line in Figure \ref{fig:isentrope}). Material properties were evaluated at each PT point along the isentrope to construct the adiabatic reference conditions shown in Figure \ref{fig:material-property-profile}. These reference conditions serve three main purposes: 1) initializing the PT fields and material properties in our ASPECT simulations (see Section \ref{sec:initialization-and-boundary-conditions}), 2) updating the material model during the simulations (see Section \ref{sec:material-model}), and 3) serving as a basis for computing "dynamic" quantities, such as the dynamic temperature $\hat{T} = T - \bar{T}$, dynamic pressure $\hat{P} = P - \bar{P}$, and dynamic density $\hat{\rho} = \rho - \bar{\rho}$, that quantify how much the approximate numerical solution deviates from the reference adiabatic conditions (i.e., a non-convecting ambient mantle).

![Reference material properties used in our ASPECT simulations. Profiles were computed using the BurnMan software [@cottaar2014; @myhill2023] and were based on the equations of state and thermodynamic data of @stixrude2022 for pure Mg olivine (ol) and wadsleyite (wd).](../figs/material-property-profile.png){#fig:material-property-profile}

### Initialization and Boundary Conditions {#sec:initialization-and-boundary-conditions}

Our ASPECT simulations used a 900 $\times$ 600 km rectangular model domain initialized with pure Mg olivine and wadsleyite (Figure \ref{fig:initial-setup}). "Surface" PT conditions of 10 GPa and 1706 K were applied at the top boundary such that the olivine $\Leftrightarrow$ wadsleyite transition occurs at approximately 130 km from the top boundary. The initial PT fields were then computed by numerically integrating the following equations:

\begin{equation}
  \frac{d\bar{T}}{dy} = \frac{\bar{\alpha}\, \bar{T}\, g}{\bar{C}_p}
  \label{eq:adiabatic-temperature}
\end{equation}

\begin{equation}
  \frac{d\bar{P}}{dy} = \bar{\rho}\, g
  \label{eq:adiabatic-pressure}
\end{equation}

where $d\bar{P}/dy$ and $d\bar{T}/dy$ are adiabatic PT profiles applied uniformly across the model domain, $g$ is gravitational acceleration, and the density $\bar{\rho}$, thermal expansivity $\bar{\alpha}$, and specific heat capacity $\bar{C}_p$ were determined from the adiabatic reference conditions shown in Figure \ref{fig:material-property-profile}. Gaussian-shaped thermal anomalies of $\pm$ 500 K were then applied along lines of length $L$ extending from the top and bottom boundaries for slab and plume simulations, respectively.

![Initial setup for slab (top) and plume (bottom) simulations. The free-slip (left and right), open (top or bottom), and prescribed inflow (top or bottom) boundary conditions for slab and plume simulations were essentially mirrored. Prescribed inflow velocities ($\vec{v}$ = 5 cm/yr) act parallel to the thermal anomalies, which are held at a fixed temperature at the boundary. The normal hydrostatic stress component ($\sigma_{yy}$) at the open boundary ensures that outflows are driven by dynamic pressures. The top boundary has a constant "surface" PT of 10 GPa and 1706 K such that the olivine $\Leftrightarrow$ wadsleyite phase transition (dashed line) occurs at 130 km from the top boundary.](../figs/initial-setup.png){#fig:initial-setup width=50%}

Velocity and stress boundary conditions were set to ensure that constant prescribed inflows of slab or plume material were balanced by free outflows at the opposite open boundaries. Prescribed inflow velocities of $\vec{u}$ = 5 cm/yr were applied along the top (slabs) or bottom (plumes) boundaries parallel to the thermal anomalies, decaying smoothly to $\vec{u}$ = 0 cm/yr at $\pm$ 15 km from the thermal anomaly centers. The left and right boundaries are free-slip ($\sigma_{xy}$ = $\bar{P}$, $\vec{u}_x$ = 0), and the top (slab) or bottom (plume) boundaries have constant normal stress component $\sigma_{yy}$ equal to the initial hydrostatic pressure $\bar{P}$ as determined by numerical integration of Equation \ref{eq:adiabatic-pressure}. The hydrostatic stress conditions at the open boundaries (bottom for slabs; top for plumes) ensure that outflows are driven by dynamic pressures due ot convection and/or volume changes during the olivine $\Leftrightarrow$ wadsleyite phase transition.

### Material Model {#sec:material-model}

#### Material Properties {#sec:material-properties}

Material properties were updated during our ASPECT simulations by referencing the adiabatic reference conditions shown in Figure \ref{fig:material-property-profile}. Except for density, material properties received no PT corrections, effectively assuming that deviations from a non-convecting ambient mantle were negligible. For density, however, we applied a dynamic PT correction through a first-order Taylor expansion [@jarvis1980; @gassmoller2020]:

\begin{equation}
  \rho \approx \bar{\rho} + \left(\frac{\partial \bar{\rho}}{\partial P} \right)_T \Delta P + \left(\frac{\partial \bar{\rho}}{\partial T} \right)_P \Delta T
  \label{eq:density-ala-expansion}
\end{equation}

Equation \ref{eq:density-ala-expansion} is rewritten using standard thermodynamic relations $\beta = \frac{1}{\rho} \left(\frac{\partial \rho}{\partial P}\right)_T$ and $\alpha = -\frac{1}{\rho} \left(\frac{\partial \rho}{\partial T}\right)_P$ to obtain the expression:

\begin{equation}
  \rho = \bar{\rho} \left(1 + \bar{\beta}\, \hat{P} - \bar{\alpha}\, \hat{T} \right)
  \label{eq:density-ala}
\end{equation}

where $\bar{\rho}$, $\, \bar{\beta}$, $\, \bar{\alpha}$, are the adiabatic reference density, compressibility, and thermal expansivity, respectively, and $\Delta P = \hat{P} = P - \bar{P}$ and $\Delta T = \hat{T} = T - \bar{T}$ are the dynamic PT. Note that the reference thermal conductivity $\bar{k}$ = 4.0 Wm$^{-1}$K$^{-1}$ is constant in all our numerical experiments.

#### Phase Transition Kinetics {#sec:phase-transition-kinetics}

The kinetics of the olivine $\Leftrightarrow$ wadsleyite phase transition were governed entirely by interface-controlled growth, as nucleation was assumed to saturate rapidly and did not limit the phase transition [@cahn1956]. Following @faccenda2017, the transformed volume fraction is given by:

\begin{equation}
  X = 1 - \exp\left(-N\, \dot{x}\, t \right)
  \label{eq:volume-fraction}
\end{equation}

where $X$ is the volume fraction of the product phase (olivine or wadsleyite), $N$ is a geometric factor that accounts for nucleation sites, $\dot{x}$ is the growth rate, and $t$ is the elapsed time after site saturation. For inter-crystalline grain-boundary controlled growth, $N = 6.67/d$, where $d$ is grain size.

Since we assumed interface-controlled growth kinetics, the following expression determined the overall phase transition rate [@hosoya2005]:

\begin{equation}
  \dot{x} = A\, T\, C_\mathrm{OH}^n\, \exp\left(-\frac{H^{\ast} + P V^{\ast}}{R\, T}\right) \left(1 - \exp\left[-\frac{\Delta G}{R\, T}\right] \right)
  \label{eq:growth-rate}
\end{equation}

where $A$ is a kinetic prefactor, $C_\mathrm{OH}$ is the concentration of water in the reactant phase, $n$ is the water content exponent, $H^{\ast}$ is activation enthalpy, $V^{\ast}$ is activation volume, $P$ is pressure, $T$ is temperature, $R$ is the gas constant, and $\Delta G$ is the Gibbs free energy difference between olivine and wadsleyite, which is approximated by:

\begin{equation}
  \Delta G \approx \Delta \bar{G} + \hat{P}\, \Delta \bar{V} - \hat{T}\, \Delta \bar{S}
  \label{eq:gibbs}
\end{equation}

where $\Delta \bar{G}$, $\Delta \bar{V}$, and $\Delta \bar{S}$ are the molar Gibbs free energy, volume, and entropy differences between olivine and wadsleyite along the adiabatic reference profile (Figure \ref{fig:thermodynamic-property-profile}), respectively, and $\hat{P}$ and $\hat{T}$ are the dynamic PT.

![Reference thermodynamic properties used in our ASPECT simulations. Profiles were computed using the same methods as described in Figure \ref{fig:material-property-profile} (see Section \ref{sec:adiabatic-reference-conditions}).](../figs/thermodynamic-property-profile.png){#fig:thermodynamic-property-profile width=80%}

In this formulation, the time evolution of the olivine $\Leftrightarrow$ wadsleyite phase transition is fully described by the interplay of pressure, temperature, and kinetic parameters applied to the interface-controlled growth model, without explicit consideration of nucleation kinetics [@hosoya2005; @faccenda2017]. The macro-scale olivine $\Leftrightarrow$ wadsleyite phase transition rate was therefore computed by taking the time derivative of Equation \ref{eq:volume-fraction}:

\begin{equation}
  \frac{dX}{dt} = \dot{X} = N\, \dot{x}\, \left(1 - X \right)
  \label{eq:phase-transition-rate-short}
\end{equation}

To simplify our numerical implementation of Equation \ref{eq:phase-transition-rate-short}, we combined the parameters $N$, $A$, and $C_\mathrm{OH}^n$ into a single kinetic factor $Z = \frac{6.67}{d}\, A\, C_\mathrm{OH}^n$. Thus, the full expression for the phase transition rate is:

\begin{equation}
  \dot{X} = Z\, T\, \exp\left(-\frac{H^{\ast} + P V^{\ast}}{R\, T}\right) \left(1 - \exp\left[-\frac{\Delta G}{R\, T}\right] \right)\, \left(1 - X \right)
  \label{eq:phase-transition-rate}
\end{equation}

The range of kinetic factors $Z$ used in our numerical experiments (3.0e0–7.0e7 K s$^{-1}$) is consistent with the kinetic growth rate parameters reported in @hosoya2005. We determined the range of $Z$ by holding $A$ = $e^{-18}$ m s$^{-1}$ K$^{-1}$ ppm$_\mathrm{OH}^{-n}$, $H^{\ast}$ = 274 kJ mol$^{-1}$, $V^{\ast}$ = 3.0e-6 m$^3$ mol$^{-1}$, and $n$ = 3.2 constant, while varying water content $C_\mathrm{OH}$ from 50–5000ppm and grain size $d$ from 1–10 mm.

#### Operator Splitting {#sec:operator-splitting}

Since the phase transition rate $\dot{X}$ is faster than the advection timescale in our ASPECT simulations, we employ a first-order operator splitting scheme to decouple advection from interface-controlled growth kinetics. In this approach, the phase fraction $X$ is updated in two sequential steps within each overall time step $\Delta t$:

\begin{equation}
  \frac{\partial X}{\partial t} + \vec{u} \cdot \nabla X = 0
  \label{eq:composition}
\end{equation}

  1. **Advection step:** Solve the transport of material $\left(\vec{u} \cdot \nabla X \right)$ without phase changes over the time interval $\Delta t$ to yield an intermediate composition $X^\ast$
  2. **Reaction step:** Starting from $X^\ast$, integrate Equation \ref{eq:phase-transition-rate} over the same time interval using a smaller sub-step $\delta t \le \Delta t$ to obtain the updated composition $X^{n+1}$

This operator splitting scheme ensures numerical stability while accurately capturing fast kinetics without restricting the convective timestep [@aspectmanual].

### Rheological Model {#sec:rheological-model}

Our ASPECT simulations use a simple rheological model where mantle viscosity is modified by a depth-dependent piecewise function:

\begin{equation}
  \bar{\eta} =
  \begin{cases}
    1\, \eta_0, & y \leq 130\, \text{km} \\
    1\, \eta_0, & y > 130\, \text{km}
  \end{cases}
  \label{eq:piecewise-viscosity-function}
\end{equation}

where $\eta_0$ = 1e21 is the nominal background viscosity of the upper mantle [@ranalli1995; @karato2008], which we assume is identical for pure olivine and wadsleyite. A thermal dependency is implemented through an exponential term:

\begin{equation}
  \eta = \bar{\eta} \exp \left(-B\, \frac{\hat{T}}{\bar{T}} \right)
  \label{eq:rheological-model}
\end{equation}

where $B$ = 1 is the thermal viscosity exponent factor, and $\hat{T}$ is the dynamic temperature.

\cleardoublepage

# Results {#sec:results}

## Simulation Snapshots: Slabs and Plumes {#sec:simulation-snapshots}

Figures \ref{fig:slab-composition-set2} and \ref{fig:plume-composition-set2} illustrate how thermodynamic driving forces, phase transition rates, and compressible flow interact to shape the 410 discontinuity. These snapshots, taken after 100 Ma of evolution, provide visual context for the quantitative analysis in Section \ref{sec:410-displacement-width}.

In slab simulations, ultra-sluggish kinetics (Figure \ref{fig:slab-composition-set2}a–c) allow metastable olivine to persist deep into the transition zone. This inhibition causes the slab to stagnate and pond, depressing the 410. Within the cold, metastable olivine region, Gibbs free energy accumulates and wadsleyite saturation remains low until the thermodynamic driving force overcomes kinetic barriers. Once this threshold is reached, the olivine $\Leftrightarrow$ wadsleyite reaction rapidly completes, producing a sharp 410 that is displaced downwards by tens of kilometers.

At intermediate phase transition rates (Figure \ref{fig:slab-composition-set2}d–f), the olivine $\Leftrightarrow$ wadsleyite reaction still lags but is fast enough to limit widespread olivine metastability and avoid total slab stagnation. The resulting 410 is broad and diffuse, as density and seismic velocity contrasts gradually fade with depth. This moderately-sluggish kinetic regime produces complex 410 structures through intermediate phase transition rates, incomplete slab stagnation, and deflected flow patterns. However, when phase transition rates are sufficiently fast to maintain quasi-equilibrium conditions (Figure \ref{fig:slab-composition-set2}g–i), the 410 sharpens and rapid wadsleyite growth within the slab allows continuous slab descent through the 410 without hesitation.

In plume simulations, thermodynamics dominate mantle flow dynamics and 410 structure. Even under ultra-sluggish kinetics (Figure \ref{fig:plume-composition-set2}a–c), the high temperatures of upwellings prevent significant olivine metastablility. The olivine $\Leftrightarrow$ wadsleyite transition proceeds rapidly, maintaining thin, sharp 410 interfaces and strong density and seismic contrasts. Although ultra-sluggish kinetics slightly broaden and uplift the 410, reducing buoyancy contrasts, plume structures remain vertically coherent across the full range of tested kinetic factors (Figure \ref{fig:plume-composition-set2}).

Altogether, these simulations demonstrate that in cold environments, kinetics strongly influence slab dynamics and control whether the 410 appears as a diffuse, low-amplitude feature or as a sharp, high-contrast seismic boundary. In contrast, thermodynamics dominate in hot plume environments, producing stable, sharply defined 410s that are largely independent of kinetic factors.

![Slab simulation within ultra-sluggish (a–c: $Z$ = 3.0e0 K s$^{-1}$), intermediate (d–f: $Z$ = 4.7e2 K s$^{-1}$), and quasi-equilibrium (g–i: $Z$ = 7.0e7 K s$^{-1}$) kinetic regimes after 100 Ma evolution. Panels show dynamic temperature $\hat{T}$ (left column), dynamic density $\hat{\rho}$ (middle column), and pressure-wave velocity $V_p$ (right column). Additional sets of visualized outputs are shown in Appendix \ref{sec:simulation-snapshots-continued}.](../figs/simulation/compositions/slab-3.0e0-4.7e2-7.0e7-set2-composition-0100.png){#fig:slab-composition-set2}

![Plume simulation within ultra-sluggish (a–c: $Z$ = 3.0e0 K s$^{-1}$), intermediate-sluggish (d–f: $Z$ = 4.7e2 K s$^{-1}$), and quasi-equilibrium (g–i: $Z$ = 7.0e7 K s$^{-1}$) kinetic regimes after 100 Ma evolution. Panels show dynamic temperature $\hat{T}$ (left column), dynamic density $\hat{\rho}$ (middle column), and pressure-wave velocity $V_p$ (right column). Additional sets of visualized outputs are shown in Appendix \ref{sec:simulation-snapshots-continued}.](../figs/simulation/compositions/plume-3.0e0-4.7e2-7.0e7-set2-composition-0100.png){#fig:plume-composition-set2}

\cleardoublepage

## Structure of the 410: Displacement and Width {#sec:410-displacement-width}

Figure \ref{fig:410-structure} summarizes the quantitative relationships between 410 structure and the maximum phase transition rate $\dot{X}_{\mathrm{max}}$ evaluated in slab and plume simulations after 100 Ma of evolution. The results reveal fundamentally different responses of plumes and slabs to thermodynamic driving forces, phase transition rates, and compressible flow.

In plume simulations, the 410 shows little systematic variation with $\dot{X}_{\mathrm{max}}$. Its structure remains nearly constant across seven orders of magnitude variation in $\dot{X}_{\mathrm{max}}$, with consistent displacements of -26 km and widths between 2–5 km. The only exception is a few ultra-sluggish kinetic models where both displacement and width increase slightly to -16 and 19 km, respectively (Table \ref{tbl:depth-profile-summary}). The weak dependence of 410 structure on $\dot{X}_{\mathrm{max}}$ reflects the strong thermal control of the reaction front in upwellings, where high temperatures promote rapid wadsleyite $\Leftrightarrow$ olivine transition---maintaining a sharp discontinuity regardless of the kinetic factor $Z$ applied to the interface-controlled growth model (Equation \ref{eq:phase-transition-rate}).

In slab simulations, the 410 exhibits distinct structural changes across three kinetic regimes. At high phase transition rates ($Z$ = 2.4e6 K s$^{-1}$; $\dot{X}_{\mathrm{max}}$ > 10$^1$ Ma$^{-1}$), the olivine $\Leftrightarrow$ wadsleyite transition remains near thermodynamic equilibrium, producing a narrow 410 (< 5 km), displaced 36–39 km upwards within the slab's inner core. As $\dot{X}_{\mathrm{max}}$ decreases (4.7e2 < $Z$ < 2.4e6 K s$^{-1}$; 10$^{-1.5}$ < $\dot{X}_{\mathrm{max}}$ < 10$^1$ Ma$^{-1}$), the 410 deepens and widens, forming a log-linear relationship where reductions in $\dot{X}_{\mathrm{max}}$ progressively broaden the reaction front (Table \ref{tbl:depth-profile-summary}). This intermediate kinetic regime corresponds to a partially inhibited olivine $\Leftrightarrow$ phase transition that proceeds slowly, hindering downward flow without complete slab stagnation.

At the lowest phase transition rates ($Z$ < 4.7e2 K s$^{-1}$; $\dot{X}_{\mathrm{max}}$ < 10$^{-1.5}$ Ma$^{-1}$), a third kinetic regime emerges. While the 410 is displaced downwards, its width narrows with further reductions in $\dot{X}_{\mathrm{max}}$. This ultra-sluggish kinetic regime reflects a transition to strong disequilibrium conditions and complete slab stagnation, where reaction progress becomes localized within the cold, high-pressure regions of the ponding slab. The apparent sharpening of the 410 under these conditions is due to a narrow portion of the phase transition zone remaining reactive, while surrounding regions are kinetically frozen.

In summary, 410 structure near plumes is regulated by thermal effects near thermodynamic equilibrium, whereas 410 structure near slabs exhibits distinct kinetic thresholds and non-linear scaling between its width, displacement, and the phase transition rate $\dot{X}$. These contrasting behaviors underscore the differing roles of kinetics in hot versus cold mantle environments and imply that the 410 beneath slabs can transition abruptly between thermodynamically- and kinetically controlled regimes as phase transition rates decrease.

![Quantitative relationships between 410 structure and maximum phase transition rates evaluated in plume and slab simulations after 100 Ma. Structure of the 410 near plumes (left column) shows minimal dependence on $\dot{X}_{\mathrm{max}}$, with both displacement and width remaining nearly constant across seven orders of magnitude variation in $\dot{X}_{\mathrm{max}}$. In contrast, 410 structure near slabs (right column) changes distinctly across three kinetic regimes: (1) quasi-equilibrium at high $\dot{X}_{\mathrm{max}}$, where 410 widths are narrow and displacements positive; (2) an intermediate regime where decreasing phase transition rates $\dot{X}_{\mathrm{max}}$ progressively widen and deepen the 410; and (3) an ultra-sluggish regime at low $\dot{X}_{\mathrm{max}}$, where the 410 narrows while deepening, and slabs completely stall and pond above -100 km displacement.](../figs/410-structure.png){#fig:410-structure width=70%}

\cleardoublepage

# Discussion {#sec:discussion}

Our numerical simulations show that interface-controlled growth kinetics, when coupled to a compressible treatment of mantle flow, exert a first-order control on the geometry and seismic expression of the 410 discontinuity. The results presented in Section \ref{sec:results} demonstrate that plume and slab flow dynamics respond in systematically different ways: plumes are insensitive to kinetics due to high temperatures, whereas slabs show three distinct kinetic regimes with thresholded behavior (Figure \ref{fig:410-structure} and Table \ref{tbl:depth-profile-summary}). The implications of such contrasting relationships in slabs versus plumes are discussed below.

## Uncertainties and Model Limitations {#sec:uncertainties-and-model-limitations}

The primary quantitative uncertainty in our analysis arises from the kinetic factor $Z$ in the interface-controlled growth model governing the olivine $\Leftrightarrow$ wadsleyite phase transition (Equation \ref{eq:phase-transition-rate}). This factor spans several orders of magnitude, reflecting a large range of water contents (50–5000 ppm) and grain sizes (1–10 mm) that impact phase transition rates. These ranges reflect experimental conditions from @hosoya2005, previous numerical studies of metastable olivine wedges [@rubie1994], and typical grain sizes of upper mantle xenoliths [~3–10 mm, @karato1984; @karato2008]. Although we explicitly hold the other kinetic parameters in Equation \ref{eq:growth-rate} constant, laboratory studies show large uncertainties for $n$, $A$, $H^{\ast}$, and $V^{\ast}$ that depend strongly on water content, grain size, Mg-Fe composition, and microstructural evolution [@rubie1994; @kubo2004; @hosoya2005; @perrillat2013; @ledoux2023]. Our simulations therefore explore only a limited subset of the potential phase transition rates in Earth's upper mantle.

Our kinetic model also relies on a key simplification: we assume instantaneous nucleation site saturation followed by interface-controlled growth (Equations \ref{eq:volume-fraction}–\ref{eq:phase-transition-rate}), thereby neglecting nucleation kinetics. This approach is often justified because nucleation rates typically occur too rapidly to be measured reliably [@hosoya2005; @kubo2004; @faccenda2017; @perrillat2016]. However, recent in-situ X-ray and acoustic studies [@ohuchi2022; @ledoux2023] document complex nucleation-growth microstructures that can limit net phase transition rates under some PT conditions. Therefore, our saturated nucleation assumption generally overestimates phase transition rates and consequently underestimates olivine metastability and its effects on flow dynamics and 410 structure.

Finally, compositional and rheological simplifications make our simulations less representative of natural behavior. Assuming pure Mg-rich end-members neglects Fe-partitioning and minor-element effects that can shift equilibrium depths by ~10–20 km and alter kinetics [@katsura2004; @perrillat2013; @perrillat2016]. We also neglect the direct effects of non-hydrostatic stress on microstructures and solid-state reactions [@wheeler2014; @wheeler2018; @wheeler2020]. Moreover, our use of a simple temperature-dependent viscosity neglects crucial rheological complexities such as grain-size evolution, plastic deformation, and stress-dependent rheologies---factors known to modify strain localization [@karato2001] and thus impact dynamic feedbacks controlling 410 structure. These omissions merit future investigation but suggest that the quantitative thresholds reported here capture the primary, first-order effects.

## Implications for Subduction Dynamics {#sec:implications-for-subduction-dynamics}

The three kinetic regimes identified in our slab simulations---quasi-equilibrium, intermediate, and ultra-sluggish---provide a framework for understanding how phase transition kinetics control slab penetration and deep earthquake activity in the mantle transition zone.

The absence of widespread slab ponding at the 410 in seismic tomography [@fukao2013] constrains the permissible kinetic conditions in Earth's mantle. The ultra-sluggish kinetic regime ($\dot{X}$ < 10$^{-1.5}$ Ma$^{-1}$), which produces complete stagnation in our simulations, appears inconsistent with tomographic images that show continuous slab descent through the 410. Most subduction zones must therefore experience sufficiently rapid olivine $\Leftrightarrow$ wadsleyite phase transition rates to avoid complete stagnation, placing an upper bound on the degree of metastability that can develop during typical subduction.

Yet this constraint alone cannot fully explain subduction zone behavior. Deep earthquakes attributed to transformational faulting [@green1995; @kirby1996; @ishii2021; @ohuchi2022; @sindhusuta2025] require olivine persistence well into the wadsleyite stability field---evidence that metastability does occur to a significant degree. Our simulations suggest this behavior is consistent with the intermediate kinetic regime (10$^{-1.5}$ < $\dot{X}$ < 10$^1$ Ma$^{-1}$), which generates localized buoyant regions of metastable olivine that resist, but do not prevent, downward flow. The coexistence of deep seismicity with continued slab penetration therefore requires a delicate balance: phase transition rates must be slow enough to sustain metastable volumes sufficient for transformational faulting, yet fast enough to permit overall slab descent through the 410.

The dynamic threshold near $\dot{X}$ ~ 10$^{-1.5}$ Ma$^{-1}$ therefore represents a critical boundary for mantle convection. Below this threshold, buoyancy forces from incomplete phase transition overwhelm slab pull, arresting descent. Above it, phase transitions proceed rapidly enough to permit 410 penetration despite temporary kinetic resistance. The sensitivity of this threshold---occurring over less than one order of magnitude in $\dot{X}$---suggests that individual subduction zones could oscillate between penetration and temporary stagnation as thermal or kinematic conditions evolve, potentially explaining temporarily stalled slabs that subsequently resume descent [@agrusta2017].

Observed variability in slab behavior beneath different subduction zones, ranging from rapid penetration to temporary stalling [@fukao2013], likely reflects regional differences influencing effective phase transition rates. The nonlinear scaling between 410 structure and $\dot{X}$ in our simulations implies that relatively modest variations in slab thermal structure, descent velocity, water content, or grain size could position different subduction zones in different kinetic regimes. Young, hot slabs descending slowly may maintain quasi-equilibrium conditions, while old, cold slabs descending rapidly may approach the intermediate regime where metastability becomes significant. This framework suggests that subduction zone diversity arises not only from differences in plate age and convergence rate but also from kinetically controlled feedbacks between phase transitions and flow dynamics.

The fundamental asymmetry between plume and slab responses to kinetics has important implications for numerical modeling. Thermal dominance in hot upwellings renders plume dynamics insensitive to kinetics, whereas the low temperatures in slabs amplify kinetic effects to the point of controlling flow regime. This asymmetry implies that numerical geodynamic models assuming thermodynamic equilibrium systematically underestimate the role of metastability in subduction. Realistic simulations of whole-mantle convection and long-term plate tectonic evolution must therefore account for kinetically controlled phase transitions, particularly in cold downwelling environments where disequilibrium effects are most pronounced.

## Implications for 410 Detectability {#sec:410-detectability}

The relationships between phase transition kinetics and 410 structure established in our simulations provide testable predictions for seismic detectability. Sharp interfaces with widths of a few kilometers are readily detected with SS precursors and receiver function stacks [@shearer2000; @chambers2005; @deuss2009], though regional high-frequency receiver function approaches can resolve structures as thin as ~5 km under favorable conditions [@helffrich1996; @wei2017; @dokht2016; @frazer2023].

Our plume simulations predict consistently thin, easily detectable discontinuities in hot upwelling environments. Across seven orders of magnitude variation in $\dot{X}$, plumes maintain 410 widths of 2–5 km with displacements of approximately -26 km, except under ultra-sluggish kinetics where widths broaden to ~19 km (Table \ref{tbl:depth-profile-summary}). These sharp discontinuities beneath plume regions should be readily detectable with standard seismic methods, consistent with observations of well-defined 410s beneath hotspots [@deuss2009; @lawrence2008]. However, composition, anisotropy, and imaging methods also influence observed sharpness, complicating unique interpretations.

Slab simulations predict more variable seismic detectability across the three kinetic regimes. In the quasi-equilibrium regime ($\dot{X}$ > 10$^1$ Ma$^{-1}$), narrow 410s with positive displacements produce sharp, high-amplitude, readily detectible seismic signals. In the intermediate regime (10$^{-1.5}$ < $\dot{X}$ < 10$^1$ Ma$^{-1}$), progressive broadening systematically reduces detectability as the phase transition front becomes diffuse, potentially rendering the 410 invisible to lower-frequency seismic methods. In the ultra-sluggish regime ($\dot{X}$ < 10$^{-1.5}$ Ma$^{-1}$), re-sharpening produces detectable 410s despite substantial deepening (> 100 km displacement), reflecting localized phase transitioning in ponded slabs.

These predictions offer a quantitative framework for interpreting observed seismic heterogeneity. Reported 410 thickness variations ranging from ~5–30 km in Pacific regions [@alex2004; @schmerr2007] can arise from spatial changes in effective phase transition rates superimposed on thermal and compositional heterogeneity. Broad, weakened 410 signals beneath some subduction zones [@lee2014; @vanstiphout2019; @jiang2015; @shen2020; @han2021] are consistent with intermediate kinetic conditions where partial metastability produces diffuse reaction fronts. Conversely, sharp 410s in cold slabs suggest either quasi-equilibrium maintained by rapid kinetics or localized phase transitioning in ultra-sluggish conditions---though the latter should be accompanied by slab stagnation, which is rarely observed at the 410 [@fukao2013].

The intermediate kinetic regime presents a particular challenge for seismic detection. Where 410 widths exceed ~10–20 km, the gradual density and velocity gradients may produce weak or absent reflections in SS precursor and receiver function studies, even during substantial phase transitioning. Such "invisible" 410s could be misinterpreted as evidence for compositional anomalies or unusual thermal structures when they actually reflect kinetically controlled broadening. High-resolution tomographic studies that image continuous velocity gradients rather than discrete discontinuities may better detect these diffuse transition zones.

Within this kinetic framework, discontinuities observed near 500–600 km depth [@cottaar2016; @saikia2008; @deuss2001; @tauzin2017] are not interpreted as metastable olivine persisting past the 410. These features likely reflect akimotoite formation, the wadsleyite $\Leftrightarrow$ ringwoodite transition, or more complex compositional layering. Broad or absent discontinuities within the mantle transition zone may indicate sluggish kinetics in addition to thermal or compositional influences, but distinguishing these effects requires independent temperature and composition constraints from complementary geophysical observations.

## Implications for Constraining Kinetic Parameters from Seismic Observations {#sec:constraining-kinetic-parameters-from-seismic-observations}

The contrasting sensitivities of plume and slab 410 structures to phase transition kinetics (Figure \ref{fig:410-structure}) suggest different strategies for extracting kinetic constraints from seismic observations.

For plumes, the near-independence of 410 structure from $\dot{X}$ limits the utility of seismic observations for constraining kinetics. The consistent 410 widths of 2–5 km and displacements of approximately -26 km observed in our simulations are primarily controlled by thermal structure rather than kinetic rates. Only in the ultra-sluggish regime, where 410 widths broaden to ~19 km and displacements decrease to -16 km, do kinetic effects become seismically distinguishable. However, this regime represents extreme kinetic inhibition unlikely to be widespread in hot upwelling environments.

For slabs, the three distinct kinetic regimes (Figure \ref{fig:410-structure}) provide exploitable diagnostic signatures. The threshold behavior near $\dot{X}$ ~ 10$^{-1.5}$ Ma$^{-1}$, which separates slab penetration from ponding, offers a critical constraint: widespread observation of slabs penetrating the 410 in global tomography [@fukao2013] require effective phase transition rates exceeding this threshold in most subduction zones. Regional variations in 410 topography and sharpness can then be interpreted within the intermediate kinetic regime (10$^{-1.5}$ < $\dot{X}$ < 10$^1$ Ma$^{-1}$), where the log-linear relationship between 410 width and $\dot{X}$ provides a potential inversion target.

Operationally, constraining kinetic parameters requires combining high-resolution seismic imaging with independent thermal constraints. Regional receiver function or SS precursor studies mapping spatial variations in 410 thickness and displacement [@chambers2005; @deuss2009; @lawrence2008; @schmerr2007; @houser2010] can be compared against tomographic temperature estimates and our simulation results to infer order-of-magnitude bounds on $\dot{X}$. Such inversions can potentially discriminate between grain sizes or dry versus wet environments following the formulation of @hosoya2005, though quantitative constraints require reducing uncertainties in kinetic parameters through targeted mineral physics experiments [e.g., @hosoya2005; @kubo2004; @perrillat2013; @perrillat2016; @ledoux2023].

Regions where thermal structure is relatively well-constrained but seismic expression varies systematically are particularly valuable---for example, across slabs with different ages, descent velocities, or hydration states [@agius2017; @vanstiphout2019; @schmandt2012]. In such settings, thermal and compositional effects can be approximately controlled, allowing lateral variations in 410 structure to isolate kinetic influences. Similarly, comparing 410 structure in hotspots versus ambient mantle, where temperature contrasts are independently estimated, can test predictions of kinetically controlled broadening and displacement in plumes [@chambers2005; @jenkins2016; @glasgow2024].

The primary limitation of seismic inversions for kinetics remains the uncertainties in the kinetic parameters themselves (Equation \ref{eq:growth-rate}). Our simulations span a broad range of $Z$ values reflecting a range of water contents and grain sizes, but the other critical kinetic parameters $n$, $A$, $H^{\ast}$, and $V^{\ast}$ remain largely unconstrained and were not tested here. Until mineral physics experiments reduce these uncertainties, seismic observations can provide order-of-magnitude estimates of effective $\dot{X}$ rather than precise determinations of micro-scale kinetic parameters. Nevertheless, the threshold behavior and scaling relationships emerging from slab simulations clarify how phase transition kinetics can account for seismic heterogeneity.

\cleardoublepage

# Conclusions {#sec:conclusions}

The olivine $\Leftrightarrow$ wadsleyite phase transition in Earth's upper mantle is strongly influenced by coupled feedbacks among thermodynamic driving forces, phase transition rates, and compressible flow. We quantified how these factors contribute to 410 expression by coupling an interface-controlled growth model to compressible simulations of mantle plumes and slabs, systematically exploring 410 structure sensitivity to kinetic factors spanning seven orders of magnitude.

Our results reveal fundamentally different responses in hot versus cold mantle environments. Plume simulations produce consistently sharp discontinuities, implying that seismic observations of 410s beneath hotspots primarily constrain thermal structure near thermodynamic equilibrium rather than phase transition rates. Slab simulations, in contrast, exhibit distinct threshold behavior across three kinetic regimes: quasi-equilibrium, intermediate, and ultra-sluggish. Widespread observations of slabs penetrating the 410 in seismic tomography suggest effective phase transition rates exceed the ultra-sluggish threshold ($\dot{X}$ > 10$^{-1.5}$ Ma$^{-1}$), while regional variations in 410 topography and sharpness likely reflect intermediate kinetic conditions (10$^{-1.5}$ Ma$^{-1}$ < $\dot{X}$ < 10$^1$ Ma$^{-1}$) modulated by local thermal structure, water content, and grain size.

These findings demonstrate that phase transition kinetics exert first-order control on slab dynamics and 410 structure beneath subduction zones. The threshold separating penetrating from ponded slabs near $\dot{X}$ ~10$^{-1.5}$ Ma$^{-1}$ directly links mantle flow to micro-scale kinetic parameters and suggests that modest variations in effective growth rates can produce substantial diversity in observed 410 topography. The 410 can therefore serve as a seismological probe of upper mantle kinetics, particularly in cold subduction environments where kinetic effects are amplified. Realizing this potential requires reducing uncertainties in kinetic parameters through targeted mineral physics experiments that better quantify nucleation versus growth mechanisms, water and compositional effects, and microstructural evolution during deformation. Integrating such experimental constraints with high-resolution seismic imaging and the forward modeling framework presented here offers a pathway toward understanding how micro-scale kinetic processes shape Earth's large-scale interior structure.

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

where $\sigma_{ij}$ is total stress, $\sigma^{\prime}_{ij}$ is the deviatoric (non-hydrostatic) component of stress, $P = - \frac{\sigma_{xx} + \sigma_{yy} + \sigma_{zz}}{3}$ is the hydrostatic component of stress, and $\delta_{ij}$ is the Kronecker delta:

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

![Slab simulation with ultra-sluggish kinetics ($Z$ = 3.0e0 K s$^{-1}$) after 100 Ma evolution. Panels show dynamic temperature $\hat{T}$ (a), dynamic pressure $\hat{P}$ (b), dynamic density $\hat{\rho}$ (c), Arrhenius term $\exp\left(-H^{\ast} + PV^{\ast}/RT\right)$ (d), thermodynamic term $\left(1 - \left[\Delta G/R\,T\right]\right)$ (e), phase transition rate $\dot{X}$ (f), volume fraction of wadsleyite $X$ (g), pressure-wave velocity $V_p$ (h), and shear-wave velocity $V_s$ (i).](../figs/simulation/compositions/slab-3.0e0-full-set-composition-0100.png){#fig:slab-composition-slow}

![Slab simulation with moderately-sluggish kinetics ($Z$ = 4.7e2 K s$^{-1}$) after 100 Ma evolution. Panels show dynamic temperature $\hat{T}$ (a), dynamic pressure $\hat{P}$ (b), dynamic density $\hat{\rho}$ (c), Arrhenius term $\exp\left(-H^{\ast} + PV^{\ast}/RT\right)$ (d), thermodynamic term $\left(1 - \left[\Delta G/R\,T\right]\right)$ (e), phase transition rate $\dot{X}$ (f), volume fraction of wadsleyite $X$ (g), pressure-wave velocity $V_p$ (h), and shear-wave velocity $V_s$ (i).](../figs/simulation/compositions/slab-4.7e2-full-set-composition-0100.png){#fig:slab-composition-moderately-slow}

![Slab simulation with quasi-equilibrium kinetics ($Z$ = 7.0e7 K s$^{-1}$) after 100 Ma evolution. Panels show dynamic temperature $\hat{T}$ (a), dynamic pressure $\hat{P}$ (b), dynamic density $\hat{\rho}$ (c), Arrhenius term $\exp\left(-H^{\ast} + PV^{\ast}/RT\right)$ (d) thermodynamic term $\left(1 - \left[\Delta G/R\,T\right]\right)$ (e), phase transition rate $\dot{X}$ (f), volume fraction of wadsleyite $X$ (g), pressure-wave velocity $V_p$ (h), and shear-wave velocity $V_s$ (i).](../figs/simulation/compositions/slab-7.0e7-full-set-composition-0100.png){#fig:slab-composition-fast}

![Plume simulation with ultra-sluggish kinetics ($Z$ = 3.0e0 K s$^{-1}$) after 100 Ma evolution. Panels show dynamic temperature $\hat{T}$ (a), dynamic pressure $\hat{P}$ (b), dynamic density $\hat{\rho}$ (c), Arrhenius term $\exp\left(-H^{\ast} + PV^{\ast}/RT\right)$ (d) thermodynamic term $\left(1 - \left[\Delta G/R\,T\right]\right)$ (e), phase transition rate $\dot{X}$ (f), volume fraction of olivine $X$ (g), pressure-wave velocity $V_p$ (h), and shear-wave velocity $V_s$ (i).](../figs/simulation/compositions/plume-3.0e0-full-set-composition-0100.png){#fig:plume-composition-slow}

![Plume simulation with moderately-sluggish kinetics ($Z$ = 4.7e2 K s$^{-1}$) after 100 Ma evolution. Panels show dynamic temperature $\hat{T}$ (a), dynamic pressure $\hat{P}$ (b), dynamic density $\hat{\rho}$ (c), Arrhenius term $\exp\left(-H^{\ast} + PV^{\ast}/RT\right)$ (d) thermodynamic term $\left(1 - \left[\Delta G/R\,T\right]\right)$ (e), phase transition rate $\dot{X}$ (f), volume fraction of olivine $X$ (g), pressure-wave velocity $V_p$ (h), and shear-wave velocity $V_s$ (i).](../figs/simulation/compositions/plume-4.7e2-full-set-composition-0100.png){#fig:plume-composition-moderately-slow}

![Plume simulation with quasi-equilibrium kinetics ($Z$ = 7.0e7 K s$^{-1}$) after 100 Ma evolution. Panels show dynamic temperature $\hat{T}$ (a), dynamic pressure $\hat{P}$ (b), dynamic density $\hat{\rho}$ (c), Arrhenius term $\exp\left(-H^{\ast} + PV^{\ast}/RT\right)$ (d) thermodynamic term $\left(1 - \left[\Delta G/R\,T\right]\right)$ (e), phase transition rate $\dot{X}$ (f), volume fraction of olivine $X$ (g), pressure-wave velocity $V_p$ (h), and shear-wave velocity $V_s$ (i).](../figs/simulation/compositions/plume-7.0e7-full-set-composition-0100.png){#fig:plume-composition-fast}

\cleardoublepage

## Structure of the 410: Displacement and Width Continued {.unnumbered #sec:410-displacement-width-continued}

The 410 width and displacement were evaluated from the phase fraction field $X$ along a vertical profile that intersected the widest 410 found within the central third of model domain ($x$ = 450 $\pm$ 140 km), where 410 width was defined as the difference between the depths at $X$ = 0.9 and $X$ = 0.1, and 410 displacement was defined as the offset between the nominal equilibrium reaction depth and the depth at $X$ = 0.9. The maximum phase transition rate, $\dot{X}_{\mathrm{max}}$, was evaluated from the phase transition rate field $\dot{X}$ along the same vertical profile.

| Simulation   | $Z$   | Displacement | Width | Max $\vec{u}_y$ | Log$_{10}$ $\dot{X}_{\mathrm{max}}$ |
|:-------------|------:|-------------:|------:|----------------:|------------------------------------:|
| plume        | 3.0e0 |          -16 |    19 |            1.92 |                               -2.87 |
| plume        | 7.0e0 |          -22 |    14 |            4.69 |                               -2.22 |
| plume        | 1.6e1 |          -23 |    10 |            4.69 |                               -1.65 |
| plume        | 3.7e1 |          -24 |     8 |            4.69 |                               -0.98 |
| plume        | 8.7e1 |          -25 |     5 |            4.69 |                               -0.31 |
| plume        | 2.0e2 |          -25 |     4 |            4.64 |                                0.26 |
| plume        | 4.7e2 |          -25 |     4 |            4.69 |                                0.81 |
| plume        | 1.1e3 |          -26 |     4 |            4.49 |                                1.39 |
| plume        | 2.6e3 |          -26 |     3 |            4.64 |                                2.49 |
| plume        | 6.0e3 |          -26 |     3 |            4.69 |                                3.37 |
| plume        | 1.4e4 |          -26 |     2 |            4.69 |                                3.84 |
| plume        | 3.3e4 |          -26 |     3 |            4.49 |                                3.46 |
| plume        | 7.8e4 |          -26 |     2 |            4.49 |                                4.59 |
| plume        | 1.8e5 |          -26 |     2 |            4.69 |                                4.95 |
| plume        | 4.3e5 |          -26 |     3 |            3.19 |                                4.84 |
| plume        | 1.0e6 |          -26 |     2 |            4.69 |                                5.60 |
| plume        | 2.4e6 |          -26 |     2 |            4.69 |                                5.97 |
| plume        | 5.6e6 |          -26 |     3 |            4.69 |                                5.48 |
| plume        | 1.3e7 |          -26 |     3 |            4.69 |                                5.73 |
| plume        | 3.0e7 |          -26 |     3 |            4.69 |                                6.27 |
| plume        | 7.0e7 |          -26 |     3 |            4.69 |                                6.93 |
| slab         | 3.0e0 |          -83 |     7 |            1.68 |                               -2.07 |
| slab         | 7.0e0 |          -66 |     8 |            1.66 |                               -1.74 |
| slab         | 1.6e1 |          -49 |    10 |            2.17 |                               -1.59 |
| slab         | 3.7e1 |          -32 |    16 |            0.61 |                               -1.46 |
| slab         | 8.7e1 |          -42 |    24 |            2.33 |                               -1.42 |
| slab         | 2.0e2 |          -52 |    63 |            0.61 |                               -1.55 |
| slab         | 4.7e2 |         -145 |   168 |            1.77 |                               -1.42 |
| slab         | 1.1e3 |         -120 |   142 |            1.91 |                               -0.86 |
| slab         | 2.6e3 |          -76 |   100 |            2.06 |                               -0.53 |
| slab         | 6.0e3 |          -42 |    68 |            2.21 |                               -0.39 |
| slab         | 1.4e4 |          -17 |    47 |            2.43 |                               -0.21 |
| slab         | 3.3e4 |            1 |    32 |            2.66 |                               -0.02 |
| slab         | 7.8e4 |           14 |    22 |            2.82 |                                0.14 |
| slab         | 1.8e5 |           23 |    15 |            2.95 |                                0.33 |
| slab         | 4.3e5 |           29 |     9 |            2.99 |                                0.64 |
| slab         | 1.0e6 |           33 |     8 |            2.98 |                                0.66 |
| slab         | 2.4e6 |           36 |     5 |            3.01 |                                0.89 |
| slab         | 5.6e6 |           37 |     4 |            3.10 |                                1.10 |
| slab         | 1.3e7 |           37 |     3 |            3.11 |                                1.91 |
| slab         | 3.0e7 |           39 |     2 |            3.23 |                                3.00 |
| slab         | 7.0e7 |           39 |     2 |            3.25 |                                3.37 |

Table: Summary of the kinetic factor $Z$, 410 displacement, 410 width, maximum vertical velocity $\vec{u}_y$, and maximum phase transition rate $\dot{X}$ evaluated in plume and slab simulations after 100 Ma of evolution. Units are $Z$: K s$^{-1}$, 410 displacement: km, 410 width: km, $\vec{u}_y$: cm/yr, Log$_{10}$ $\dot{X}_{\mathrm{max}}$: Log$_{10}$ Ma$^{-1}$. {#tbl:depth-profile-summary}
