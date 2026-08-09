<!--
# Abstract {.unnumbered #sec:abstract}

The seismic expression of Earth's 410 km discontinuity varies across tectonic settings, from sharp, high-amplitude interfaces to broad transitions---patterns that cannot be explained by equilibrium thermodynamics without invoking large-scale thermal or compositional heterogeneities. Laboratory experiments show the olivine $\Leftrightarrow$ wadsleyite transition responsible for the 410 is rate-limited, yet previous numerical studies have not directly evaluated the sensitivity of 410 structure to kinetic and rheological factors. Here we investigate these relationships by coupling a grain-scale, interface-controlled olivine $\Leftrightarrow$ wadsleyite growth model to compressible simulations of mantle plumes and subducting slabs. We vary kinetic parameters across seven orders of magnitude and quantify the resulting 410 displacements and widths. Our results reveal an asymmetry between hot and cold environments. In plumes, high temperatures produce sharp 410s (2--3 km wide) regardless of kinetics. In slabs, kinetics exert first-order control on 410 structure through three regimes: (1) quasi-equilibrium conditions producing narrow, uplifted 410s and continuous slab descent; (2) intermediate reaction rates generating broader, deeper 410s with metastable olivine wedges resisting downward slab motion; and (3) ultra-sluggish reaction rates causing slab stagnation with re-sharpened, deeply displaced 410s ($\lesssim$ 100 km). Rheological contrasts modulate these kinetic effects by controlling slab geometry and residence time in the phase transition zone. These findings demonstrate that reaction rates strongly influence 410 structure in subduction zones, establishing the 410 as a potential seismological constraint on upper mantle kinetic processes, particularly in cold environments where disequilibrium effects are amplified.

# Plain Language Summary {.unnumbered #sec:plain-language-summary}

Seismic imaging reveals that Earth's 410 km discontinuity---a boundary in Earth's mantle where olivine transforms into a denser form called wadsleyite---varies significantly, sometimes appearing sharp and other times as a broad zone. While often explained through temperature and compositional differences, laboratory experiments show that this mineral transformation takes time to complete, suggesting that reaction speed (kinetics) is important. We simulated rising mantle plumes and descending slabs to explore how kinetics and rock strength affect the 410 km discontinuity. Our models reveal striking differences. In hot plumes, high temperatures allow fast reactions, maintaining a consistently sharp discontinuity (2--3 km thick). In cold subducting slabs, however, kinetics are critical. Fast reactions produce a narrow, uplifted discontinuity as slabs sink smoothly. Moderate reaction rates create broader, deeper discontinuities with wedges of unreacted olivine that slow downward slab motion. Extremely slow reactions cause slabs to stagnate with sharp 410 km discontinuities that are shifted downwards by up to 100 km. These findings suggest that observations of the 410 km discontinuity, particularly in subduction zones, could provide valuable constraints on reaction rates deep within Earth's mantle.

# Keypoints {.unnumbered #sec:key-points}

- Plumes produce sharp 410s regardless of reaction rates; slabs show three distinct kinetic regimes controlling discontinuity structure
- Rheology modulates kinetic effects: strong slabs descend slowly allowing complete reaction; weak slabs descend fast amplifying metastability
- Seismic observations of 410 structure in subduction zones can constrain reaction rates but require independent rheological constraints
-->

\clearpage

# Definition of Symbols {.unnumbered #sec:symbols}

|Parameter|Symbol|Unit|Equations|
|:--------------|:------|:----|:----|
|Activation energy|$E^{\ast}$|J mol$^{-1}$|\ref{eq:arrhenius-viscosity}|
|Activation enthalpy|$H^{\ast}$|J mol$^{-1}$|\ref{eq:growth-rate}, \ref{eq:reaction-rate}|
|Activation factor (rheology)|$B$|-|\ref{eq:rheological-model}|
|Activation volume|$V^{\ast}$|m$^3$ mol$^{-1}$|\ref{eq:growth-rate}, \ref{eq:reaction-rate}|
|Compressibility (reference)|$\bar{\beta}$|Pa$^{-1}$|\ref{eq:density-ala}|
|Density|$\rho$|kg m$^{-3}$|\ref{eq:navier-stokes-no-inertia}--\ref{eq:continuity-expanded}, \ref{eq:density-ala-expansion}--\ref{eq:density-ala}|
|Density (reference)|$\bar{\rho}$|kg m$^{-3}$|\ref{eq:adiabatic-pressure}--\ref{eq:density-ala}|
|Density (dynamic)|$\hat{\rho}$|kg m$^{-3}$|-|
|Deviatoric stress tensor|$\sigma^{\prime}$|Pa|\ref{eq:navier-stokes-no-inertia}, \ref{eq:energy}|
|Deviatoric strain rate tensor|$\dot{\epsilon}^{\prime}$|s$^{-1}$|\ref{eq:energy}|
|Gas constant|$R$|J mol$^{-1}$ K$^{-1}$|\ref{eq:growth-rate}, \ref{eq:reaction-rate}, \ref{eq:arrhenius-viscosity}--\ref{eq:rheological-model}|
|Grain size|$d$|m|-|
|Gravitational acceleration|$g$|m s$^{-2}$|\ref{eq:navier-stokes-no-inertia}, \ref{eq:adiabatic-temperature}--\ref{eq:adiabatic-pressure}|
|Growth rate|$\dot{x}$|m s$^{-1}$|\ref{eq:volume-fraction}--\ref{eq:growth-rate}, \ref{eq:reaction-rate-short}|
|Latent heat|$Q_L$|J kg$^{-1}$|\ref{eq:energy}|
|Molar entropy|$\bar{S}$|J mol$^{-1}$ K$^{-1}$|\ref{eq:gibbs}|
|Molar Gibbs free energy|$\bar{G}$|J mol$^{-1}$|\ref{eq:gibbs}|
|Molar volume|$\bar{V}$|m$^{3}$ mol$^{-1}$|\ref{eq:gibbs}|
|Nucleation site factor|$N$|m$^{-1}$|\ref{eq:volume-fraction}, \ref{eq:reaction-rate-short}|
|Prefactor (growth rate)|$\Gamma$|m s$^{-1}$ K$^{-1}$ ppm$_\mathrm{OH}^{-n}$|\ref{eq:growth-rate}|
|Prefactor (kinetic)|$Z$|K$^{-1}$ s$^{-1}$|\ref{eq:reaction-rate}|
|Prefactor (viscosity)|$A$|Pa s|\ref{eq:arrhenius-viscosity}--\ref{eq:background-viscosity}|
|Pressure|$P$|Pa|\ref{eq:navier-stokes-no-inertia}, \ref{eq:energy}, \ref{eq:growth-rate}, \ref{eq:reaction-rate}|
|Pressure (reference)|$\bar{P}$|Pa|\ref{eq:adiabatic-pressure}|
|Pressure (dynamic)|$\hat{P}$|Pa|\ref{eq:density-ala}, \ref{eq:gibbs}|
|Reaction rate|$\frac{dX}{dt}$, $\frac{\partial X}{\partial t}$, $\dot{X}$|s$^{-1}$|\ref{eq:reaction-rate-short}--\ref{eq:composition}|
|Specific heat capacity (reference)|$\bar{C}_p$|J kg$^{-1}$ K$^{-1}$|\ref{eq:energy}, \ref{eq:adiabatic-temperature}|
|Temperature|$T$|K|\ref{eq:energy}, \ref{eq:growth-rate}, \ref{eq:reaction-rate}, \ref{eq:arrhenius-viscosity}, \ref{eq:rheological-model}|
|Temperature (reference)|$\bar{T}$|K|\ref{eq:adiabatic-temperature}, \ref{eq:arrhenius-viscosity-expanded}--\ref{eq:rheological-model}|
|Temperature (dynamic)|$\hat{T}$|K|\ref{eq:density-ala}, \ref{eq:gibbs}, \ref{eq:arrhenius-viscosity-expanded}, \ref{eq:rheological-model}|
|Thermal conductivity (reference)|$\bar{k}$|W m$^{-1}$ K$^{-1}$|\ref{eq:energy}|
|Thermal expansivity (reference)|$\bar{\alpha}$|K$^{-1}$|\ref{eq:energy}, \ref{eq:adiabatic-temperature}, \ref{eq:density-ala}|
|Time|$t$|s|\ref{eq:continuity-compressible}--\ref{eq:continuity-expanded}, \ref{eq:volume-fraction}, \ref{eq:reaction-rate-short}--\ref{eq:composition}|
|Velocity|$\vec{u}$|m s$^{-1}$|\ref{eq:continuity-compressible}--\ref{eq:continuity-expanded}, \ref{eq:composition}|
|Viscosity|$\eta$|Pa s|\ref{eq:arrhenius-viscosity}--\ref{eq:arrhenius-viscosity-expanded}, \ref{eq:rheological-model}|
|Viscosity (reference)|$\bar{\eta}$|Pa s|\ref{eq:background-viscosity}--\ref{eq:rheological-model}|
|Volume fraction|$X$|-|\ref{eq:volume-fraction},  \ref{eq:reaction-rate-short}--\ref{eq:composition}|
|Water content|$C_\mathrm{OH}$|ppm|\ref{eq:growth-rate}|
|Water content exponent|$n$|-|\ref{eq:growth-rate}|

\clearpage

# Introduction {#sec:introduction}

Earth's mantle transition zone hosts two prominent seismic discontinuities near 410 and 660 km depth, attributed to polymorphic phase transitions of olivine [@ringwood1975; @katsura2004]. While these discontinuities are observed globally, their detailed seismic characteristics---depth, sharpness, amplitude, and lateral continuity---vary substantially between tectonic settings [@deuss2009; @lawrence2008; @schmerr2007; @fukao2013]. Some regions display sharp, high-amplitude reflectors consistent with abrupt mineralogical boundaries, while others exhibit broad, weakened, or laterally variable signals. Such heterogeneity cannot be explained by equilibrium thermodynamics alone, which relates discontinuity topography mainly to temperature-dependent phase boundaries [e.g., @cottaar2016; @jenkins2016]. Additional physical processes, including reaction kinetics and dynamic pressure effects [@rubie1994; @faccenda2017], or compositional heterogeneities, including variations in water content [@karato2011; @smyth1987; @smyth2002] or bulk composition [@glasgow2024; @goes2022; @saikia2008; @schmerr2007; @tauzin2017], likely contribute to observed variability.

Laboratory studies demonstrate that the olivine $\Leftrightarrow$ wadsleyite transition at 410 km depth (the "410") occurs over finite time scales rather than instantaneously, with kinetics governed by temperature, pressure, water content, bulk chemical composition, grain size, and microstructural evolution [@rubie1994; @liu1998; @kubo2004; @hosoya2005; @perrillat2013; @ledoux2023; @mohiuddin2018]. In cold subducting slabs, sluggish reaction rates can allow metastable olivine to persist tens of kilometers below its thermodynamic stability limit, promoting slab stagnation and potentially triggering deep earthquakes via transformational faulting [@rubie1994; @green1995; @kirby1996; @ishii2021; @ohuchi2022; @sindhusuta2025]. In hot upwellings, slow kinetics may broaden and uplift the discontinuity, possibly explaining reduced seismic amplitudes beneath some hotspots [@chambers2005]. However, published kinetic models remain poorly constrained, with parameters spanning several orders of magnitude [e.g., @hosoya2005; @kerschhofer1998; @kubo2004; @mosenfelder2001], leaving the effects of micro-scale kinetic processes on flow dynamics and seismic observables ambiguous.

Bridging the gap between laboratory-derived kinetic rate laws and mantle-scale seismic observations requires numerical models that couple reaction kinetics to realistic treatments of mantle convection. Previous modeling efforts have demonstrated that kinetics can strongly influence mantle flow [@dassler1996a; @dassler1996b; @schmeling1999; @guest2004; @agrusta2017; @faccenda2017], but investigations quantifying the sensitivity of 410 structure to kinetic parameters remain limited. Moreover, most prior studies impose kinematic restrictions and/or employ simplified treatments of compressibility or kinetic rate laws that may inadequately capture feedbacks between kinetically controlled phase transitions and flow dynamics. Additionally, the role of rheology in modulating kinetic effects through its control on slab geometry and descent rate has not been thoroughly explored.

This study aims to clarify these issues by implementing a grain-scale, interface-controlled olivine $\Leftrightarrow$ wadsleyite growth model [after @hosoya2005] within compressible mantle flow simulations using the open-source geodynamic modeling software ASPECT. We systematically explore how reaction kinetics and rheological strength jointly influence 410 structure and address three specific questions:

1. How do kinetic factors impact flow dynamics and shape the 410 in hot versus cold environments?
2. How do viscosity contrasts modulate these kinetic effects?
3. Can seismic observations of 410 structure constrain effective kinetic and rheological parameters?

To investigate these questions, we analyze a suite of numerical experiments varying kinetic parameters across seven orders of magnitude. For each experiment, we test a range of viscosity contrasts and quantify 410 displacement and width, enabling direct comparisons with seismological observations. Our results establish quantitative relationships between reaction rates, rheological strength, flow dynamics, and 410 structure, demonstrating that realistic treatment of reaction kinetics is essential for accurately modeling subduction dynamics and interpreting seismic structures.

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

where $\sigma^{\prime}$ is the deviatoric stress tensor, $\rho$ is density, $g$ is gravitational acceleration, $t$ is time, $\bar{C}_p$, $\bar{k}$, $\bar{\alpha}$ are the reference specific heat capacity, thermal conductivity, and thermal expansivity, respectively (see Section \ref{sec:adiabatic-reference-conditions}), and $Q_L$ is the latent heat released or absorbed during phase transitions. Equations \ref{eq:navier-stokes-no-inertia} and \ref{eq:continuity-compressible} together describe the buoyancy-driven flow of an isotropic fluid with negligible inertia, and Equation \ref{eq:energy} describes the conduction, advection, and production (or consumption) of thermal energy [@schubert2001]. Note that the pressure $P$ in this context is equal to the mean normal stress and is positive under compression: $P = - \frac{\sigma_{xx} + \sigma_{yy} + \sigma_{zz}}{3}$.

The compressible form of the continuity equation (Equation \ref{eq:continuity-compressible}) is essential for capturing the full coupling between density changes, pressure and temperature (PT) variations, and kinetically controlled phase transitions. This is in contrast to simplified formulations such as the Boussinesq approximation or anelastic liquid approximation, which either neglect the time derivative of density ($\frac{\partial \rho}{\partial t}$) entirely or impose restrictions on where density variations can appear in the governing equations [@gassmoller2020]. By retaining $\frac{\partial \rho}{\partial t}$, the compressible continuity equation enables density changes from kinetically controlled phase transitions to influence the flow field anywhere within the model domain---not just through equilibrium thermodynamics, but through time-dependent reaction progress. This bidirectional coupling between flow dynamics and reaction kinetics is critical for accurately simulating systems where phase transitions occur over finite timescales comparable to advective timescales.

To solve Equation \ref{eq:continuity-compressible} numerically, we adopt the *projected density approximation* [PDA, @gassmoller2020], which reformulates the continuity equation by applying the product rule to $\nabla \cdot (\rho\, \vec{u})$ and multiplying both sides by $\frac{1}{\rho}$:

\begin{equation}
  \frac{1}{\rho} \frac{\partial \rho}{\partial t} + \nabla \cdot \vec{u} + \left(\frac{1}{\rho} \nabla \rho \right) \cdot \vec{u} = 0
  \label{eq:continuity-expanded}
\end{equation}

The projected density field $\rho(T, P, X)$ varies with temperature, pressure, and reaction progress, ensuring that local density changes arising from both PT variations and phase transitions influence the flow through buoyancy as well as volumetric expansion and compression. The phase transitions themselves are modeled using a separate kinetic rate law (described in Section \ref{sec:reaction-kinetics}), which determines the reaction progress $X$ based on local thermodynamic and kinetic conditions. This makes the PDA ideally suited for our numerical experiments, which incorporate density changes due to the olivine $\Leftrightarrow$ wadsleyite transition.

## Numerical Setup {#sec:numerical-setup}

### Adiabatic Reference Conditions {#sec:adiabatic-reference-conditions}

To ensure numerical convergence, we initialized our ASPECT simulations with reasonable estimates of the pressure-temperature (PT) fields and material properties in Earth's upper mantle. We began by evaluating entropy changes over a PT range of 1573--1973 K and 0.001--25 GPa (Figure \ref{fig:isentrope}) using the Gibbs free energy minimization software Perple_X [v.7.0.9, @connolly2009]. We assumed a dry pyrolitic bulk composition after @green1979, and phase equilibria were evaluated in the Na$_2$O-CaO-FeO-MgO-Al$_2$O$_3$-SiO$_2$ (NCFMAS) chemical system with thermodynamic data and solution models of @stixrude2022. Equations of state were included for solid solution phases: olivine, plagioclase, spinel, clinopyroxene, wadsleyite, ringwoodite, perovskite, ferropericlase, high-pressure C2/c pyroxene, orthopyroxene, akimotoite, post-perovskite, Ca-ferrite, garnet, and Na-Al phase.

![Entropy (left) and density (right) changes in Earth's upper mantle under thermodynamic equilibrium and hydrostatic stress conditions. Material properties were computed with Perple_X using the equations of state and thermodynamic data of @stixrude2022. The black box indicates the approximate PT range of our ASPECT simulations, while the white line indicates the isentropic adiabat used to calculate reference material properties.](../figs/PYR-material-table.png){#fig:isentrope width=65%}

We then determined the mantle adiabat by applying the Newton-Raphson algorithm to find temperatures corresponding to each pressure that maintain constant entropy (white line in Figure \ref{fig:isentrope}). Material properties were evaluated at each PT point along the isentrope to construct the adiabatic reference conditions shown in Figure \ref{fig:material-property-profile}. These reference conditions serve three main purposes: 1) initializing the PT fields and material properties in our ASPECT simulations (see Section \ref{sec:initialization-and-boundary-conditions}), 2) updating the material model during the simulations (see Section \ref{sec:material-model}), and 3) serving as a basis for computing "dynamic" quantities, such as the dynamic temperature $\hat{T} = T - \bar{T}$, dynamic pressure $\hat{P} = P - \bar{P}$, and dynamic density $\hat{\rho} = \rho - \bar{\rho}$, that quantify how much the approximate numerical solution deviates from the adiabatic reference conditions.

![Reference material properties used in our ASPECT simulations. Profiles were computed using the BurnMan software [@cottaar2014; @myhill2023] and were based on the equations of state and thermodynamic data of @stixrude2022 for pure Mg olivine (ol) and wadsleyite (wd).](../figs/material-property-profile.png){#fig:material-property-profile width=100%}

### Initialization and Boundary Conditions {#sec:initialization-and-boundary-conditions}

We set up our ASPECT simulations within a 900 $\times$ 600 km rectangular model domain, initialized with pure Mg olivine and wadsleyite (Figure \ref{fig:initial-setup}). "Surface" PT conditions of 10 GPa and 1706 K were applied at the top boundary such that the olivine $\Leftrightarrow$ wadsleyite transition occurs at approximately 130--140 km from the top boundary. The initial adiabatic PT profiles were computed by numerically integrating the following equations:

\begin{equation}
  \frac{d\bar{T}}{dy} = \frac{\bar{\alpha}\, \bar{T}\, g}{\bar{C}_p}
  \label{eq:adiabatic-temperature}
\end{equation}

\begin{equation}
  \frac{d\bar{P}}{dy} = \bar{\rho}\, g
  \label{eq:adiabatic-pressure}
\end{equation}

where the material properties $\bar{\rho}$, $\bar{\alpha}$, and $\bar{C}_p$ were determined from the adiabatic reference conditions shown in Figure \ref{fig:material-property-profile}.

![Initial setup for slab (top) and plume (bottom) simulations. The top boundary has a constant "surface" PT of 10 GPa and 1706 K such that the olivine $\Leftrightarrow$ wadsleyite transition (dashed line) occurs at 130--140 km from the top boundary. Thermal anomalies with Gaussian profiles were superimposed on top of the initial adiabatic temperature profile and remained fixed at the inflow boundary. Constant inflow velocities of 5 cm/yr were prescribed parallel to the thermal anomalies. Normal tractions equal to the initial lithostatic pressure profile were enforced on the left, right, and open boundaries to ensure that outflows are driven only by dynamic pressures.](../figs/initial-setup.png){#fig:initial-setup width=45%}

For the initial temperature field, thermal anomalies with amplitudes of $\pm$ 500 K were superimposed on the adiabatic temperature profile. These anomalies were defined as smooth linear features with Gaussian cross-sections (15 km half-width) and tanh-tapered ends (5 km taper length) to avoid sharp discontinuities. Slabs extended 100 km horizontally and 100 km downward from the top boundary; plumes extended 450 km upward from the bottom boundary. Velocity boundary conditions prescribed constant inflow of 5 cm/yr parallel to the thermal anomalies, tapering smoothly to zero at the thermal anomaly edges with the same Gaussian profile. Zero horizontal velocities were imposed at the side boundaries ($\vec{u}_x$ = 0). The full functional form of the Gaussian-tanh thermal anomalies and velocity boundary conditions is available within the accompanying data repository (see Data Availability statement for details).

Stress boundary conditions on the left and right boundaries enforced a normal traction equal to the lithostatic pressure profile computed in Equation \ref{eq:adiabatic-pressure}. No tangential (shear) stresses were applied to the side boundaries, so they approximated impermeable, free-slip surfaces under hydrostatic confinement. Open boundaries (bottom for slabs, top for plumes) were assigned a constant normal traction equal to the initial lithostatic pressure at the respective boundary $\sigma_{yy}$ = $\bar{P}(y)$. These stress conditions allow outflow to occur freely at the top or bottom boundaries (for plumes versus slabs), driven only by dynamic pressure variations associated with convection and/or volume changes during the olivine $\Leftrightarrow$ wadsleyite transition.

### Material Model {#sec:material-model}

#### Material Properties {#sec:material-properties}

Material properties were updated during our ASPECT simulations by referencing the adiabatic reference conditions shown in Figure \ref{fig:material-property-profile}. Except for density, material properties received no PT corrections, effectively assuming that deviations from the adiabatic reference conditions were negligible. For density, however, we applied a dynamic PT correction through a first-order Taylor expansion [@jarvis1980; @gassmoller2020]:

\begin{equation}
  \rho \approx \bar{\rho} + \left(\frac{\partial \bar{\rho}}{\partial P} \right)_T \Delta P + \left(\frac{\partial \bar{\rho}}{\partial T} \right)_P \Delta T
  \label{eq:density-ala-expansion}
\end{equation}

Equation \ref{eq:density-ala-expansion} is rewritten using standard thermodynamic relations $\beta = \frac{1}{\rho} \left(\frac{\partial \rho}{\partial P}\right)_T$ and $\alpha = -\frac{1}{\rho} \left(\frac{\partial \rho}{\partial T}\right)_P$ to obtain the expression:

\begin{equation}
  \rho \approx \bar{\rho} \left(1 + \bar{\beta}\, \hat{P} - \bar{\alpha}\, \hat{T} \right)
  \label{eq:density-ala}
\end{equation}

where $\bar{\rho}$, $\, \bar{\beta}$, and $\bar{\alpha}$ are the adiabatic reference density, compressibility, and thermal expansivity, respectively, and $\Delta P = \hat{P} = P - \bar{P}$ and $\Delta T = \hat{T} = T - \bar{T}$ are the dynamic PT. Note that the reference thermal conductivity $\bar{k}$ = 4.0 W m$^{-1}$K$^{-1}$ is constant in all our numerical experiments.

#### Reaction Kinetics {#sec:reaction-kinetics}

The kinetics of the olivine $\Leftrightarrow$ wadsleyite transition were governed entirely by interface-controlled growth, as nucleation was assumed to saturate rapidly and did not limit the reaction [@cahn1956]. Following @cahn1956, the transformed volume fraction is given by:

\begin{equation}
  X = 1 - \exp\!\left(-N\, \dot{x}\, t \right)
  \label{eq:volume-fraction}
\end{equation}

where $X$ is the volume fraction of the product phase (olivine or wadsleyite), $N$ is a geometric factor that accounts for nucleation sites, $\dot{x}$ is the growth rate, and $t$ is the elapsed time after site saturation. For inter-crystalline grain-boundary controlled growth, $N = 6.67/d$, where $d$ is grain size.

Since we assumed interface-controlled growth kinetics, the following expression determined the overall reaction rate [@hosoya2005]:

\begin{equation}
  \dot{x} = \Gamma\, T\, C_\mathrm{OH}^n\, \exp\!\left(-\frac{H^{\ast} + P V^{\ast}}{R\, T}\right) \left(1 - \exp\!\left[-\frac{\Delta G}{R\, T}\right] \right)
  \label{eq:growth-rate}
\end{equation}

where $\Gamma$ is the growth rate prefactor, $C_\mathrm{OH}$ is the concentration of water in the reactant phase, $n$ is the water content exponent, $H^{\ast}$ is activation enthalpy, $V^{\ast}$ is activation volume, $P$ is pressure, $T$ is temperature, $R$ is the gas constant, and $\Delta G$ is the Gibbs free energy difference between olivine and wadsleyite, which is approximated by:

\begin{equation}
  \Delta G \approx \Delta \bar{G} + \hat{P}\, \Delta \bar{V} - \hat{T}\, \Delta \bar{S}
  \label{eq:gibbs}
\end{equation}

where $\Delta \bar{G}$, $\Delta \bar{V}$, and $\Delta \bar{S}$ are the molar Gibbs free energy, volume, and entropy differences between olivine and wadsleyite along the adiabatic reference profile (Figure \ref{fig:thermodynamic-property-profile}), respectively, and $\hat{P}$ and $\hat{T}$ are the dynamic PT.

This formulation represents a significant thermodynamic simplification. We define Gibbs free energy based on mean normal stress (pressure), even though our simulations involve non-hydrostatic stress conditions where deviatoric stresses and normal stresses at individual grain interfaces directly influence reaction pathways and chemical potentials [@wheeler2015; @wheeler2020]. Mean stress can provide a reasonable approximation under certain conditions [@wheeler2015b]. A rigorous treatment will require consideration of stress-dependent chemical potentials that account for the full stress tensor's influence on phase stability and reaction kinetics (work in progress).

Setting aside these important microstructural effects for now, the time evolution of the olivine $\Leftrightarrow$ wadsleyite transition is fully described by the interplay of pressure, temperature, and kinetic parameters applied to the interface-controlled growth model, without explicit consideration of nucleation kinetics [@hosoya2005; @faccenda2017]. The macro-scale olivine $\Leftrightarrow$ wadsleyite reaction rate was therefore computed by taking the time derivative of Equation \ref{eq:volume-fraction}:

\begin{equation}
  \frac{dX}{dt} = \dot{X} = N\, \dot{x}\, \left(1 - X \right)
  \label{eq:reaction-rate-short}
\end{equation}

Since all of the kinetic parameters $N$, $\Gamma$, and $C_\mathrm{OH}^n$ ultimately scale the reaction rate (Equations \ref{eq:growth-rate}--\ref{eq:reaction-rate-short}), varying them independently adds little scientific value. Instead, we simplified our numerical implementation of Equation \ref{eq:reaction-rate-short} by combining the parameters $N$, $\Gamma$, and $C_\mathrm{OH}^n$ into a single kinetic prefactor $Z = \frac{6.67}{d}\, \Gamma\, C_\mathrm{OH}^n$. Thus, the full expression for the reaction rate became:

\begin{equation}
  \dot{X} = Z\, T\, \exp\!\left(-\frac{H^{\ast} + P V^{\ast}}{R\, T}\right) \left(1 - \exp\!\left[-\frac{\Delta G}{R\, T}\right] \right)\, \left(1 - X \right)
  \label{eq:reaction-rate}
\end{equation}

The range of kinetic prefactors $Z$ used in our numerical experiments (3.0e0--7.0e7 K$^{-1}$ s$^{-1}$) was determined by holding $\Gamma$ = $\exp\!\left(-18\right)$ m s$^{-1}$ K$^{-1}$ ppm$_\mathrm{OH}^{-n}$, $H^{\ast}$ = 274 kJ mol$^{-1}$, $V^{\ast}$ = 3.0e-6 m$^3$ mol$^{-1}$, and $n$ = 3.2 constant, while varying water content $C_\mathrm{OH}$ from 50--5000 ppm and grain size $d$ from 1--10 mm. These water contents and grain sizes are consistent with the experimental conditions of @hosoya2005, previous numerical studies of metastable olivine wedges [@rubie1994], and typical grain sizes of upper mantle xenoliths [$\sim$ 3--10 mm, @karato1984; @karato2008]. Thus, our experiments approximate kinetic conditions ranging from slow kinetics in dry rocks with large grain sizes (50 ppm OH; 10 mm; $Z$ = 3.0e0 K$^{-1}$ s$^{-1}$) to fast kinetics in hydrated rocks with small grain sizes (5000 ppm OH; 1 mm; $Z$ = 7.0e7 K$^{-1}$ s$^{-1}$).

![Reference thermodynamic properties used in our ASPECT simulations. Profiles were computed using the same methods as described in Figure \ref{fig:material-property-profile} (see Section \ref{sec:adiabatic-reference-conditions}).](../figs/thermodynamic-property-profile.png){#fig:thermodynamic-property-profile width=80%}

#### Operator Splitting {#sec:operator-splitting}

Since the reaction rate $\dot{X}$ was faster than the advection timescale in our ASPECT simulations, we employed a first-order operator splitting scheme to decouple advection from interface-controlled olivine $\Leftrightarrow$ wadsleyite growth kinetics. In this approach, the transformed volume fraction $X$ was updated in two sequential steps within each overall time step $\Delta t$:

\begin{equation}
  \frac{\partial X}{\partial t} + \vec{u} \cdot \nabla X = 0
  \label{eq:composition}
\end{equation}

  1. **Reaction step:** Starting from $X^{t}$, the reaction rate equation (Equation \ref{eq:reaction-rate}) was integrated forward in time over the full advection timestep $\Delta t$, yielding an intermediate composition $X^{\ast}$. This integration was performed using a sequence of internally determined reaction substeps $\delta t \le \Delta t$, such that the reaction operator was applied repeatedly until the cumulative integration time reached $\Delta t$.
  2. **Advection step:** Starting from $X^{\ast}$, the material transport term $\left(\vec{u} \cdot \nabla X \right)$ in Equation \ref{eq:composition} was solved over the advection time interval $\Delta t$, neglecting phase changes, to obtain the updated composition $X^{t+1}$.

In our simulations, the reaction step was solved using the ARKODE package of the SUNDIALS library [@reynolds2023; @hindmarsh2005], which employs an adaptive-step additive Runge-Kutta method [@aspectmanual]. The reaction substep size $\delta t$ was determined dynamically by ARKODE to satisfy a specified relative tolerance during the reaction step (set to $10^{-6}$ in this study). As a result, the reaction operator was subcycled as needed within each global timestep $\Delta t$, ensuring that fast solid-state reaction kinetics were accurately captured without imposing overly restrictive global timesteps.

### Rheological Model {#sec:rheological-model}

We use a temperature-dependent viscosity following an Arrhenius law:

\begin{equation}
  \eta = A \exp\!\left(\frac{E^{\ast}}{R\,T}\right)
  \label{eq:arrhenius-viscosity}
\end{equation}

with viscosity prefactor $A$, activation energy $E^{\ast}$, gas constant $R$, and absolute temperature $T$, which is decomposed into an adiabatic reference temperature $\bar{T}$ and a perturbation $\hat{T}$, $T=\bar{T}+\hat{T}$. Linearizing the Arrhenius relation through a first-order Taylor expansion about $\bar{T}$ yields:

\begin{equation}
  \eta \approx A \exp\!\left(\frac{E^{\ast}}{R\,\bar{T}}\right) \exp\!\left(-\frac{E^{\ast}}{R\,\bar{T}}\frac{\hat{T}}{\bar{T}}\right)
  \label{eq:arrhenius-viscosity-expanded}
\end{equation}

By defining a reference background viscosity as:

\begin{equation}
  \bar{\eta} = A \exp\!\left(\frac{E^{\ast}}{R\,\bar{T}}\right)
  \label{eq:background-viscosity}
\end{equation}

we arrive at a rheological model where viscosity varies exponentially with thermal perturbations about an adiabatic reference profile:

\begin{equation}
  \eta \approx \bar{\eta}\,\exp\!\left(-B\,\frac{\hat{T}}{\bar{T}}\right) \qquad B=\frac{E^{\ast}}{R\,\bar{T}}
  \label{eq:rheological-model}
\end{equation}

In our simulations, we prescribed a uniform background viscosity $\bar{\eta}$ = 10$^{21}$ Pa s throughout the upper mantle [e.g., @karato2008; @ranalli1995], implicitly assuming that $E^{\ast}$ was temperature-dependent and varied slightly with depth (Equation \ref{eq:background-viscosity}). We varied the rheological activation factor $B$ between 2 (low thermal sensitivity) and 10 (high thermal sensitivity). For our prescribed thermal anomalies of $\hat{T}$ = $\pm$ 500 K, these $B$ values correspond to approximately 2--20$\times$ stronger slabs or 2--20$\times$ weaker plumes relative to the background mantle for $B$ = 2 and 10, respectively. Thus, our numerical experiments explore a range of viscosity contrasts between thermal anomalies and background adiabatic reference conditions, with the exponential rheological model creating symmetric viscosity variations for heating and cooling.

### Numerical Stabilization of Dynamic Pressure Oscillations {#sec:numerical-stability}

A significant numerical limitation emerges from the coupled feedback between our kinetic rate law (Equations \ref{eq:volume-fraction}--\ref{eq:reaction-rate}) and the buildup of dynamic pressure in the fully compressible continuity equation (Equation \ref{eq:continuity-compressible}). At sufficiently high rheological contrasts ($B$ $\gtrsim$ 4), cold slabs develop internal dynamic pressures that can exceed several hundred MPa. These dynamic pressure perturbations accelerate the forward reaction (through the Gibbs free energy term in Equation \ref{eq:reaction-rate}), causing rapid local density changes. These density changes then alter the pressure field through the coupled continuity and momentum equations (Equations \ref{eq:navier-stokes-no-inertia}--\ref{eq:continuity-compressible}), which in turn drives the reverse reaction. This positive feedback loop manifests as spurious pressure waves propagating through the slab interior.

To address this issue, we adopt an approach similar to @gassmoller2020 and exclude the dynamic pressure contribution from the Gibbs free energy calculation (Equation \ref{eq:gibbs}) while retaining it in the Arrhenius term in Equation \ref{eq:reaction-rate} and density formulation (Equation \ref{eq:density-ala}). In practice, this means replacing $\Delta G \approx \Delta \bar{G} + \hat{P}\, \Delta \bar{V} - \hat{T}\, \Delta \bar{S}$ with $\Delta G \approx \Delta \bar{G} - \hat{T}\, \Delta \bar{S}$ in our kinetic rate law. As discussed by @gassmoller2020, this approximation is justified because density variations from dynamic pressure are typically small compared to those from compositional heterogeneities and temperature variations in mantle convection. However, this treatment does introduce a limitation. In scenarios where dynamic pressure effects dominate over thermal and compositional contributions, such as in regions with extreme viscosity contrasts or near phase boundaries with large $\Delta \bar{V}$, our kinetic model may underestimate the thermodynamic driving force for phase transitions. Nevertheless, this compromise ensures numerical stability while preserving the essential physics of kinetically controlled phase transitions coupled to compressible flow.

# Results {#sec:results}

## Simulation Snapshots: Slabs and Plumes {#sec:simulation-snapshots}

Figures \ref{fig:slab-composition-set2} and \ref{fig:plume-composition-set2} illustrate how reaction rates impact dynamic flow patterns and shape the 410 discontinuity. These snapshots, taken after 100 Ma of evolution, provide visual context for the quantitative analysis in Section \ref{sec:410-displacement-width}.

In slab simulations, ultra-sluggish kinetics (Figure \ref{fig:slab-composition-set2}: top row) allow metastable olivine to persist deep into the transition zone. This inhibition causes the slab to stagnate and pond, depressing the 410. Within the cold, metastable olivine region, Gibbs free energy accumulates and wadsleyite saturation remains low until the thermodynamic driving force overcomes kinetic barriers. Once this threshold is reached, the olivine $\Leftrightarrow$ wadsleyite reaction rapidly completes, producing a sharp 410 that is displaced downwards by tens of kilometers.

At intermediate reaction rates (Figure \ref{fig:slab-composition-set2}: middle row), the olivine $\Leftrightarrow$ wadsleyite reaction still lags but is fast enough to limit widespread olivine metastability and avoid total slab stagnation. The resulting 410 is broad and diffuse, as density and seismic velocity contrasts gradually fade with depth. This moderately sluggish kinetic regime produces complex 410 structures through intermediate reaction rates, incomplete slab stagnation, and deflected flow patterns. However, when reaction rates are sufficiently fast to maintain quasi-equilibrium conditions (Figure \ref{fig:slab-composition-set2}: bottom row), the 410 sharpens and shifts upwards as expected from equilibrium thermodynamics. Under this fast kinetic regime, rapid wadsleyite growth within the slab allows continuous slab descent through the 410 without hesitation.

In plume simulations, thermal effects dominate mantle flow dynamics and 410 structure. Even under ultra-sluggish kinetics (Figure \ref{fig:plume-composition-set2}: top row), the high temperatures of upwellings prevent significant wadsleyite metastability. The olivine $\Leftrightarrow$ wadsleyite transition proceeds rapidly, maintaining thin, sharp 410 interfaces and strong density and seismic contrasts. Although ultra-sluggish kinetics slightly broaden and uplift the 410---reducing buoyancy contrasts---plume structures remain vertically coherent across the full range of tested kinetic prefactors $Z$ (Figure \ref{fig:plume-composition-set2}).

Altogether, these simulations demonstrate that in cold environments, kinetics strongly influence slab dynamics and control whether the 410 appears as a diffuse, low-amplitude feature or as a sharp, high-contrast seismic boundary. In contrast, thermal effects dominate in hot plume environments, producing stable, sharply defined 410s that are largely independent of tested kinetic prefactors.

![Slab simulations with intermediate strength contrasts ($B$ = 4, viscosity contrast $\sim$ 3$\times$) showing ultra-sluggish (top row: $Z$ = 3.0e0 K$^{-1}$ s$^{-1}$), intermediate (middle row: $Z$ = 4.7e2 K$^{-1}$ s$^{-1}$), and quasi-equilibrium (bottom row: $Z$ = 7.0e7 K$^{-1}$ s$^{-1}$) kinetic regimes after 100 Ma evolution. Panels show dynamic temperature $\hat{T}$ (left column), dynamic density $\hat{\rho}$ (middle column), and pressure-wave velocity $V_p$ (right column). Thin lines highlight the 10% and 90% wadsleyite volume fraction contours ($X$ = 0.1 and 0.9). The 410 displacement is defined as the difference between the depth at X = 0.9 and the nominal equilibrium olivine $\Leftrightarrow$ wadsleyite transition depth, while the 410 width is defined as the difference between depths at X = 0.9 and X = 0.1 (see Supplementary Information for details). The white arrows (right column) indicate where the 410 structure was measured.](../figs/simulation/compositions/slab-Z3.0e0-B4-Z4.7e2-B4-Z7.0e7-B4-set2-composition-0010.png){#fig:slab-composition-set2 width=100%}

![Plume simulations with intermediate strength contrasts ($B$ = 4, viscosity contrast $\sim$ 1/3$\times$) showing ultra-sluggish (top row: $Z$ = 3.0e0 K$^{-1}$ s$^{-1}$), intermediate (middle row: $Z$ = 4.7e2 K$^{-1}$ s$^{-1}$), and quasi-equilibrium (bottom row: $Z$ = 7.0e7 K$^{-1}$ s$^{-1}$) kinetic regimes after 100 Ma evolution. Panels show dynamic temperature $\hat{T}$ (left column), dynamic density $\hat{\rho}$ (middle column), and pressure-wave velocity $V_p$ (right column). Thin lines highlight the 10% and 90% wadsleyite volume fraction contours ($X$ = 0.1 and 0.9). The 410 displacement is defined as the difference between the depth at X = 0.9 and the nominal equilibrium olivine $\Leftrightarrow$ wadsleyite transition depth, while the 410 width is defined as the difference between depths at X = 0.9 and X = 0.1 (see Supplementary Information for details). The white arrows (right column) indicate where the 410 structure was measured.](../figs/simulation/compositions/plume-Z3.0e0-B4-Z4.7e2-B4-Z7.0e7-B4-set2-composition-0010.png){#fig:plume-composition-set2 width=100%}

## Structure of the 410: Displacement and Width {#sec:410-displacement-width}

Figure \ref{fig:410-structure} summarizes the quantitative relationships between 410 structure and the maximum reaction rate $\dot{X}_{\mathrm{max}}$ evaluated in slab and plume simulations after 100 Ma of evolution. The results reveal fundamentally different responses of plumes and slabs to reaction kinetics. See the Supplementary Information for the full set of experimental results and technical details describing the 410 structural measurements.

![Measured 410 displacement and width versus maximum reaction rates $\dot{X}_{\mathrm{max}}$ in plume and slab simulations with intermediate strength contrasts ($B$ = 4) after 100 Ma. Structure of the 410 near plumes (left column) shows minimal dependence on $\dot{X}_{\mathrm{max}}$, with both displacement and width remaining nearly constant across seven orders of magnitude variation in $\dot{X}_{\mathrm{max}}$. In contrast, 410 structure near slabs (right column) changes distinctly across three kinetic regimes: (1) quasi-equilibrium at high $\dot{X}_{\mathrm{max}}$, where 410 widths are narrow and elevated; (2) an intermediate regime where decreasing reaction rates $\dot{X}_{\mathrm{max}}$ progressively widen and deepen the 410; and (3) an ultra-sluggish regime at low $\dot{X}_{\mathrm{max}}$, where the 410 narrows while deepening, and slabs completely stall and pond.](../figs/410-structure-B4.png){#fig:410-structure width=70%}

In plume simulations, the 410 shows little dependence on $\dot{X}_{\mathrm{max}}$. Its structure remains nearly constant across seven orders of magnitude variation in $\dot{X}_{\mathrm{max}}$, with consistent displacements of 24 km and widths between 2--9 km (Figure \ref{fig:410-structure}). The only exception is a few ultra-sluggish kinetic models where displacements decrease to 18--21 km and widths increase to 12--21 km (e.g., Figure \ref{fig:plume-composition-set2}: top row). The weak dependence of 410 structure on $\dot{X}_{\mathrm{max}}$ reflects the strong thermal control of the reaction front in upwellings, where high temperatures promote rapid wadsleyite $\Leftrightarrow$ olivine transition---maintaining a sharp discontinuity regardless of the kinetic prefactor $Z$ applied to the interface-controlled olivine $\Leftrightarrow$ wadsleyite growth model (Equation \ref{eq:reaction-rate}). Similarly, variations in rheological strength contrast (controlled by the $B$ parameter in Equation \ref{eq:rheological-model}) produce negligible changes to plume morphology or 410 structure, as the thermally dominated behavior is largely insensitive to viscosity variations (see Supplementary Information Figures S6--S8 for examples).

In slab simulations, the 410 exhibits distinct structural changes across three kinetic regimes (cf. regions (1)--(3) in Figure \ref{fig:410-structure}). At high reaction rates ($Z$ $\gtrsim$ 1.8e5 K$^{-1}$ s$^{-1}$; $\dot{X}_{\mathrm{max}}$ $\gtrsim$ 2.2 Ma$^{-1}$), the olivine $\Leftrightarrow$ wadsleyite transition remains near thermodynamic equilibrium, producing a narrow 410 ($\lesssim$ 10 km), displaced 27--39 km upwards within the slab's inner core (region (1) in Figure \ref{fig:410-structure}). As $\dot{X}_{\mathrm{max}}$ decreases (2.0e2 $\lesssim$ $Z$ $\lesssim$ 1.8e5 K$^{-1}$ s$^{-1}$; 0.07 $\lesssim$ $\dot{X}_{\mathrm{max}}$ $\lesssim$ 2.2 Ma$^{-1}$), the 410 deepens and widens, approximating a log-linear relationship where reductions in $\dot{X}_{\mathrm{max}}$ progressively broaden the reaction front (region (2) in Figure \ref{fig:410-structure}). This intermediate kinetic regime corresponds to a partially inhibited olivine $\Leftrightarrow$ wadsleyite transition that proceeds slowly, hindering downward flow without complete slab stagnation.

At the lowest reaction rates ($Z$ $\lesssim$ 2.0e2 K$^{-1}$ s$^{-1}$; $\dot{X}_{\mathrm{max}}$ $\lesssim$ 0.07 Ma$^{-1}$), a third kinetic regime emerges. While the 410 is displaced downwards, its width narrows with further reductions in $\dot{X}_{\mathrm{max}}$ (region (3) in Figure \ref{fig:410-structure}). This ultra-sluggish kinetic regime reflects a transition to strong disequilibrium conditions associated with complete slab stagnation and the persistence of metastable olivine to greater depths. As metastable olivine ponds and is pushed deeper, $\Delta G$ rapidly increases and the thermodynamic driving force term $\left(1 - \exp\left[- \Delta G / R\, T \right]\right)$ in Equation \ref{eq:reaction-rate} approaches unity, such that the olivine $\Leftrightarrow$ wadsleyite reaction rate becomes effectively independent of further increases in disequilibrium. In this limit, the kinetics are controlled primarily by the Arrhenius term $\exp\left(-(H^{\ast} + P\, V^{\ast}) / R\, T\right)$ together with the saturation term $\left(1 - X\right)$ in Equation \ref{eq:reaction-rate}. Consequently, the abrupt transformation of metastable olivine reflects the attainment of a critical temperature where the Arrhenius term increases sharply. Once this threshold is reached, reaction rates accelerate rapidly and the transformation proceeds over a narrow depth range while colder surrounding regions remain kinetically frozen (e.g., Figure \ref{fig:slab-composition-set2}: top row). The apparent 410 sharpening results from this Arrhenius-controlled reaction front, localized at the depth where temperature becomes sufficient to activate rapid kinetics under highly overstepped conditions.

Rheological strength contrasts substantially modulate slab dynamics and 410 structure through interactions with kinetic effects (Figure \ref{fig:410-structure-comp}). Higher $B$ values in Equation \ref{eq:rheological-model} increase the viscosity contrast between the cold slab and surrounding mantle, producing a more coherent slab that penetrates the 410 with less internal deformation. Stronger slabs sink more sub-horizontally through the 410, limiting their vertical descent rates. Therefore, as $B$ increases, the 410 sharpens because material traverses the phase transition zone more slowly, allowing greater time for reaction progress despite sluggish kinetics. In contrast, lower $B$ values permit greater internal deformation and faster vertical descent, which amplify kinetic effects and produce broader phase transition zones (see Supplementary Information Figures S3--S5 for examples).

![Variation in 410 structure and slab descent rate across kinetic and rheological regimes. Panels show maximum reaction rate (top left), maximum vertical velocity (top right), 410 width (bottom left), and 410 displacement (bottom right) as functions of the kinetic prefactor $Z$ (horizontal axis, log scale) and rheological activation factor $B$ (vertical axis). Each colored tile represents the measured value within a simulation after 100 Ma, and black/white lines delineate transitions between regime behaviors. The precise location of these transitions depends on $B$, where increasing rheological contrast progressively shifts the regime boundaries towards lower $Z$ (more sluggish kinetic conditions). Text labels represent the observed kinetic regimes: (1) quasi-equilibrium, (2) intermediate, and (3) ultra-sluggish. The complete dataset is given in the Supplementary Information.](../figs/410-structure-comp.png){#fig:410-structure-comp width=70%}

In summary, 410 structure near plumes is regulated by thermal effects near thermodynamic equilibrium, whereas 410 structure near slabs exhibits distinct kinetic thresholds and non-linear scaling between its width, displacement, and the reaction rate $\dot{X}$. These contrasting behaviors underscore the differing roles of kinetics in hot versus cold mantle environments and imply that the 410 beneath slabs can transition abruptly between thermodynamically and kinetically controlled regimes as reaction rates decrease.

# Discussion {#sec:discussion}

Our numerical simulations demonstrate that in cold subduction environments, the structure and seismic expression of the 410 discontinuity are strongly governed by the Arrhenius temperature-dependence of the olivine $\Leftrightarrow$ wadsleyite reaction rate. Moreover, the Arrhenius law in Equation \ref{eq:reaction-rate} explains the fundamental asymmetry between slabs and plumes observed in our results. Cold slab temperatures suppress reaction rates by many orders of magnitude, generating pronounced metastability, while hot plume temperatures drive rapid transformation toward quasi-equilibrium regardless of the kinetic prefactor $Z$. The results presented in Section \ref{sec:results} demonstrate these systematically different responses. Plumes maintain near-equilibrium conditions across large variations in $Z$, whereas slabs exhibit three distinct kinetic regimes with threshold behavior (Figure \ref{fig:410-structure}). Rheological strength contrasts, controlled by the activation factor $B$, further modulate these kinetic effects by altering slab geometry and descent rate (Figure \ref{fig:410-structure-comp}). The implications of these contrasting behaviors are discussed below.

Before a detailed discussion, however, an important conceptual distinction must be drawn between the "kinetic effects" described here and "thermal effects" as described in equilibrium-based models. Equilibrium-based interpretations of 410 structure relate discontinuity topography to temperature through the Clapeyron slope of the olivine $\Leftrightarrow$ wadsleyite phase boundary (+2.5--4 MPa/K), predicting approximately 3--9 km of 410 uplift per 100 K of cooling [@bina1994; @katsura2004; @flanagan1999]. Compositional effects further modify this interpretation. Fe-Mg partitioning can shift equilibrium depths by $\sim$ 10--20 km [@katsura2004; @helffrich1996; @perrillat2013], water content shifts the phase boundary to lower P and widens the transition [@smyth1987; @smyth2002; @vandermeijde2003], and the presence of harzburgitic and basaltic lithologies in subducting slabs introduces phase assemblages that depart substantially from the olivine system [@ringwood1991; @xu2008; @tauzin2022; @waszek2021; @glasgow2024]. These "thermal" and "compositional effects" on the 410 are well-established.

Our study addresses a complementary hypothesis: that olivine $\Leftrightarrow$ wadsleyite kinetics can independently account for observed diversity in 410 structure. In our kinetic framework, no Clapeyron slope is prescribed. Temperature enters Equation \ref{eq:reaction-rate} through a linear term $T$, the Arrhenius term $\exp\left(-(H^{\ast} + P\, V^{\ast}) / R\, T\right)$, and the thermodynamic driving force term $\left(1 - \exp\left[- \Delta G / R\, T \right]\right)$. As a result, "kinetic effects" and "thermal effects" cannot be meaningfully decoupled within our model. What appears as 410 topography is the direct manifestation of cold temperatures suppressing reaction rates, rather than shifting an equilibrium phase boundary. To isolate this kinetic signal without conflating it with the compositional, water, and grain-size effects explored elsewhere [@katsura2004; @hosoya2005; @glasgow2024], our simulations deliberately adopt a simplified single-component composition and temperature-dependent rheology. The limitations of applying these simplifications are discussed in the following section.

## Uncertainties and Model Limitations {#sec:uncertainties-and-model-limitations}

The primary quantitative uncertainty in our analysis stems from the kinetic prefactor $Z$ in the interface-controlled olivine $\Leftrightarrow$ wadsleyite growth model (Equation \ref{eq:reaction-rate}), which spans several orders of magnitude reflecting variable water contents (50--5000 ppm) and grain sizes (1--10 mm). Water content and grain size were intentionally absorbed into the single kinetic prefactor $Z$ in our formulation (Section \ref{sec:reaction-kinetics}), allowing characterization of first-order kinetic behavior in terms of "fast versus slow" kinetics without presupposing specific water or grain-size distributions. While the OH-dependence of the transformation rate [@hosoya2005] and the thermodynamic effect of water on the phase boundary [@smyth1987; @smyth2002] both motivate explicit treatment of water content, mantle water concentrations and their spatial distribution within subducted slabs remain uncertain [@houser2016; @karato2011; @ishii2021; @cerpa2022]. Laboratory studies also reveal large uncertainties in kinetic parameters $n$, $\Gamma$, $H^{\ast}$, and $V^{\ast}$ that depend strongly on water content, grain size, Mg-Fe composition, and microstructural evolution [@rubie1994; @liu1998; @kubo2004; @hosoya2005; @perrillat2013; @ledoux2023]. Our simulations therefore explore only a limited subset of potential reaction rates in Earth's upper mantle, and we interpret our results within the conceptual framework of "fast versus slow" kinetics rather than "wet versus dry" conditions and/or "fine versus coarse" grain sizes.

Our kinetic model also assumes instantaneous nucleation site saturation followed by interface-controlled growth (Equations \ref{eq:volume-fraction}--\ref{eq:reaction-rate}), thereby neglecting nucleation kinetics. While often justified because nucleation rates occur too rapidly for reliable measurement [@hosoya2005; @kubo2004; @faccenda2017; @perrillat2016], recent in-situ X-ray and acoustic studies [@ohuchi2022; @ledoux2023] document complex nucleation-growth microstructures that can limit net reaction rates under some PT conditions. By assuming saturated nucleation, we generally overestimate reaction rates and consequently underestimate the degree of olivine metastability and its influence on flow dynamics and 410 structure.

Beyond these kinetic uncertainties, our simulations incorporate several compositional simplifications. Assuming pure Mg endmembers neglects Fe-partitioning and minor element effects that can shift equilibrium depths by $\sim$ 10--20 km and alter kinetics [@katsura2004; @perrillat2013; @perrillat2016]. We also neglect non-hydrostatic stress effects on microstructures and solid-state reactions [@wheeler2014; @wheeler2018; @wheeler2020]. These simplifications are intentional. By adopting a single-component, monomineralic (Mg$_2$SiO$_4$) composition, we isolate the role of interface-controlled kinetics from the effects of Fe-Mg partitioning, water content, and lithological mixing that have been thoroughly characterized elsewhere [@helffrich1996; @katsura2004; @xu2008; @glasgow2024]. Our results represent a phenomenological assessment of kinetic contributions to the 410 structure, rather than direct quantitative predictions that can be applied to specific natural subduction zones.

Our simple temperature-dependent rheological model (Equations \ref{eq:arrhenius-viscosity}--\ref{eq:rheological-model}) omits several strain-dependent processes that could affect both flow dynamics and reaction kinetics. Notably, accumulated strain could influence 410 topography through multiple pathways. Grain size reduction would accelerate kinetics by increasing nucleation site density (affecting $N$ in Equation \ref{eq:reaction-rate-short}) while lowering local viscosity via diffusion creep [@karato1984; @karato1986; @hirth2003; @rosa2016; @mohiuddin2018; @mohiuddin2020]. Strain localization could similarly create zones of enhanced reaction rates, and strain-weakening could modify slab descent trajectory and residence time in the transition zone. While incorporating strain-weakening rheologies, grain-size evolution, and plastic deformation would provide a more complete treatment [@agrusta2017; @karato2001], such additions would substantially expand an already large parameter space, introduce additional uncertainties in constitutive relationships, and obscure the first-order kinetic relationships quantified here.

Our model also prescribes a uniform background viscosity across the upper-to-lower mantle boundary and omits the wadsleyite $\Leftrightarrow$ ringwoodite (520 km) and ringwoodite $\Leftrightarrow$ bridgmanite + periclase (660 km) reactions. While a viscosity jump of roughly 1--2 orders of magnitude is inferred across the 660 km discontinuity [@goes2017; @faccenda2017; @karato2008], no comparable jump is expected at the 410 [@ranalli1995]. Flow below the 410 in our simulations is governed by olivine-wadsleyite buoyancy contrasts alone, and the effects of the ringwoodite $\Leftrightarrow$ bridgmanite + periclase transition on slab dynamics are absent. This means we do not reproduce the slab stagnation and folding dynamics associated with the 660 [@schubert2001; @fukao2013; @garel2014; @agrusta2017], and slab descent rates reported here should be interpreted accordingly as upper bounds in the absence of deeper phase-transition hindrance.

Finally, our 2D simulations neglect 3D slab geometry and dynamics. While we expect the same kinetic regimes to occur in 3D, the transitions between regimes will likely shift as slabs behave differently in 3D versus 2D [@sime2024]. Despite these combined limitations in kinetic parameters, compositional complexity, rheological formulation, and dimensionality, our results capture first-order kinetic effects governing 410 structure (Figures \ref{fig:410-structure} and \ref{fig:410-structure-comp}).

## Implications for Subduction Dynamics {#sec:implications-for-subduction-dynamics}

The three kinetic regimes identified in our slab simulations---quasi-equilibrium, intermediate, and ultra-sluggish---emerge from feedbacks between kinetically controlled phase transitions and slab strength (Figure \ref{fig:410-structure-comp}). However, seismic tomography shows little to no evidence for widespread slab ponding at the 410 [@fukao2013], providing an observational constraint on plausible kinetic behavior in natural systems. In particular, the ultra-sluggish regime ($\dot{X}$ $\lesssim$ 0.07 Ma$^{-1}$), which leads to complete slab stagnation at the 410 in our simulations, is inconsistent with global seismic observations of continuous slab descent through the 410. This implies that most subduction zones experience sufficiently rapid olivine $\Leftrightarrow$ wadsleyite transformation to avoid long-lived stagnation and ponding at the 410.

These kinetic regimes complement, rather than replace, equilibrium-based interpretations of 410 topography. In the quasi-equilibrium kinetic regime ($\dot{X} \gtrsim$ 2.2 Ma$^{-1}$; region (1) in Figure \ref{fig:410-structure}), our simulations approach the thermodynamic predictions based on Clapeyron slopes [@katsura2004; @flanagan1999]. The intermediate and ultra-sluggish regimes ($\dot{X} \lesssim$ 2.2 Ma$^{-1}$; regions (2) and (3) in Figure \ref{fig:410-structure}), however, produce 410 structures that depart substantially from equilibrium predictions---wider, deeper, or paradoxically re-sharpened discontinuities that cannot be reproduced by temperature or water effects alone. A broad 410 is commonly attributed to elevated water content, which expands the two-phase coexistence interval and depresses the phase boundary [@smyth2002; @vandermeijde2003; @frost2007]. Our intermediate kinetic regime (region (2) in Figure \ref{fig:410-structure}) produces qualitatively similar broadening, but with deepening rather than shallowing of the 410, potentially offering a means to distinguish kinetic inhibition from hydration effects on the basis of observed topography. Distinguishing kinetic signatures from those arising from Fe-Mg partitioning, compositional heterogeneity [@katsura2004; @helffrich1996; @xu2008], or water content [@hosoya2005; @glasgow2024] in natural subduction zones requires independent constraints on slab thermal structure, hydration state, and mineralogy [@agius2017; @vanstiphout2019; @schmandt2012].

Meanwhile, the occurrence of deep earthquakes often attributed to transformational faulting could indicate metastable olivine persistence in various subduction zones [@green1995; @kirby1996; @ishii2021; @ohuchi2022; @sindhusuta2025] and thus be compatible with the intermediate kinetic regime in our simulations (0.07 $\lesssim$ $\dot{X}$ $\lesssim$ 2.2 Ma$^{-1}$; region (2) in Figure \ref{fig:410-structure}). This regime would generate localized buoyant regions that resist but do not prevent downward flow. Assuming that transformational faulting is primarily responsible for deep earthquakes, rather than alternative mechanisms [e.g., @zhan2020], our simulations suggest that coexistence of deep seismicity with continuous slab descent through the 410 requires a delicate balance: olivine $\Leftrightarrow$ wadsleyite reaction rates must be slow enough to sustain metastable volumes for transformational faulting, yet fast enough to permit 410 penetration.

Rheological contrasts further modulate these kinetic effects by controlling slab trajectory and descent rate. Strong slabs ($B$ $\gtrsim$ 6, viscosity contrast $\sim$ 6--20$\times$) maintain coherence and penetrate sub-horizontally, spending more time within the phase transition zone where olivine transforms more completely (Figure S4: middle and bottom rows). Weaker slabs ($B$ $\lesssim$ 6, viscosity contrast $\sim$ 2--6$\times$) deform internally and descend vertically, traversing the transition zone quickly with less time for sluggish reactions to complete (Figures \ref{fig:slab-composition-set2}: middle row; S4: top row). This morphological dependence on slab strength is consistent with previous numerical studies without reaction kinetics [@garel2014], indicating that rheological control of slab geometry operates independently of kinetic effects.

Considering both kinetic and rheological influences can produce counterintuitive subduction dynamics, however. Under equivalent kinetic conditions, our simulations show stronger slabs descending through the 410, while weaker slabs accumulate more metastable olivine, amplifying buoyancy forces that slow or arrest downward slab motion (Figures \ref{fig:410-structure-comp} and S4). This outcome appears contradictory to previous work linking deep earthquakes with colder, stronger slabs [@green1995; @houston2015; @zhan2020], suggesting that other factors---such as overriding plate forcing or subduction geometry [e.g., @agrusta2017]---may dominate over kinetics in natural systems.

Our slab simulations also reveal a systematic departure from observed descent rates across the current global subduction system. Despite the absence of any hindrance from the ringwoodite $\Leftrightarrow$ bridgmanite + periclase transition in our models, simulated slabs descend at 1.3--2.0 cm/yr (interquartile range excluding fully stagnated cases; see Supplementary Information Table S1), roughly half of the range compiled by @lallemand2005 for natural systems (interquartile range: 2.6--4.8 cm/yr). This discrepancy suggests that buoyancy forces from metastable olivine accumulation may be exaggerated in our simulations relative to Earth, although faster descent rates observed in natural subduction zones should amplify metastability if kinetic rates are sluggish.

Furthermore, our simulations point to a critical kinetic threshold near $\dot{X}$ $\sim$ 0.07 Ma$^{-1}$ (the boundary between regions (2) and (3) in Figure \ref{fig:410-structure}) that marks where buoyancy forces from incomplete olivine $\Leftrightarrow$ wadsleyite transformation either prevent or permit downward slab motion. This narrow threshold shifts with slab strength, where coherent slabs penetrate the 410 at reaction rates that would stall weaker slabs (Figure \ref{fig:410-structure-comp}). These dynamic sensitivities suggest individual subduction zones could oscillate between penetration and temporary stagnation at the 410 as kinematic and PT conditions evolve, with plate forcing and slab deformation behavior [e.g., @agrusta2017] potentially modulating this critical kinetic threshold.

Regional variability in slab behavior [@fukao2013] could therefore reflect diversity in both kinetics and rheology, in addition to slab age and convergence rate. In our simulations, slabs with lower viscosity contrasts descend steeply and accumulate moderate amounts of metastable olivine, while slabs with high viscosity contrasts flatten and sink slowly, allowing near-complete transformation under the same moderately sluggish kinetic conditions. These patterns distinguish slabs from plumes. In hot upwellings, elevated temperatures suppress sensitivity to both kinetics and rheology, while in cold slabs, kinetics govern reaction rates and rheology dictates reaction duration by regulating slab kinematics. Accurately modeling slab morphology, material exchange, and deep stress conditions in geodynamic simulations therefore requires incorporating both effects.

## Implications for 410 Detectability and Constraining Kinetics and Rheology from Seismic Observations {#sec:410-detectability}

The 410 discontinuity is readily detected using underside reflections (precursors to SS, PP, and ScS reverberated phases), as well as converted phases (receiver functions), and top-side triplicated phases [e.g., @flanagan1999; @shearer2000; @chambers2005; @schmerr2007; @deuss2009; @houser2010; @cottaar2016; @jenkins2016; @wang2017; @han2021; @glasgow2024; @miao2024; @parera2024; @waszek2024]. The majority of studies focus on mapping 410 topography and typically assume a sharp near-horizontal interface. Resolving the width of the phase transition zone requires more challenging analyses, including constraints on impedance contrast. The trade-off between phase transition width, impedance contrast, and small-scale topography depends strongly on the wavelengths of the seismic phases employed [@chambers2005].

Previous studies targeting the width of discontinuities have exploited the frequency dependence of observed amplitudes [e.g., @vandermeijde2003; @bonatto2020], the detectability of discontinuities at higher frequencies [e.g., @rost2002; @day2013], gradient-based inversions [@schmandt2012], and joint analyses of top- and bottom-side reflections [@vinnik2020]. In many of these approaches, the transition is parameterized as a simple linear gradient, whereas our models predict a more complex depth-dependent gradient, consistent with results from gradient-based inversions [@schmandt2012]. While many of these studies infer a sharp 410, locally broader phase transitions, up to 35 km, are observed and are commonly attributed to elevated water content in the mantle [e.g., @vandermeijde2003; @schmandt2012; @vinnik2020].

The presence of water is expected to significantly widen the olivine $\Leftrightarrow$ wadsleyite transition width, similar to effects observed in our simulations. However, a hydrous 410 transition would likely occur at lower pressures [@smyth2002], while our simulation results show a deepening of the discontinuity for the intermediate kinetic regime with an inhibited olivine $\Leftrightarrow$ wadsleyite transformation. Therefore, resolving topography could be a means to distinguish these two scenarios.

Our plume simulations predict consistently thin 410 discontinuities (2--3 km widths, $\sim$ 24 km displacements) across large rheological variations and seven orders of magnitude in $\dot{X}$, except under ultra-sluggish kinetics where widths broaden up to 21 km. These sharp discontinuities should be readily detectable, consistent with observations of well-defined 410s with strong deflections in various hotspot locations [@deuss2009; @lawrence2008].

Slab simulations, in contrast, imply that 410 detectability strongly depends on joint kinetic and rheological influences (Figure \ref{fig:410-structure-comp}). Strong slabs ($B$ $\gtrsim$ 6, viscosity contrast $\sim$ 6--20$\times$) descending slowly along sub-horizontal trajectories in our simulations maximize residence time in the phase transition zone, producing sharp, detectable 410 signals due to near-complete olivine $\Leftrightarrow$ wadsleyite transformation (Figure S4: middle and bottom rows). Weak slabs ($B$ $\lesssim$ 6, viscosity contrast $\sim$ 2--6$\times$) descending steeply with reduced residence times amplify kinetic inhibition, producing broader reaction fronts (Figures \ref{fig:slab-composition-set2}: middle row; S4: top row). Weakened 410 signals beneath some subduction zones [@vanstiphout2019], or observations of low-velocity regions possibly representing metastable olivine [@jiang2015; @shen2020; @han2021], align with intermediate kinetic conditions and relatively weak slabs, where steep, rapid descent amplifies metastability (Figures \ref{fig:slab-composition-set2}: middle row; S4: top row).

Understanding the joint influence of rheological and kinetic factors on the geometry and dynamics of slab descent and phase transition width will offer a quantitative framework for interpreting observed seismic structure. Kinetic factors could potentially lead to significantly extended phase transition zones, but these anomalies are predicted to occur over narrow regions where slabs interact with the 410 (Figure \ref{fig:slab-composition-set2}). Thus, finer-scale topographic mapping with higher frequency phases, such as receiver functions, is likely needed to localize sharp topography or non-detectability [@vanstiphout2019].

In principle, regional variations in both 410 topography and phase transition width could reflect variations in kinetics (through temperature, composition, water content, or grain size), variations in rheology (through additional effects like accumulated strain), or compensating variations in both. Without multiple independent constraints, seismic observations of 410 structure alone cannot uniquely separate these effects. Examples of independent constraints may include the presence of deep seismicity, as well as variations in slab thermal structure, age, and hydration across subduction zones [@agius2017; @vanstiphout2019; @schmandt2012].

Forward modeling of seismically observable 410 width and topography enables testing of whether specific parameter combinations match observed 410 structure from receiver functions or precursors. However, until mineral physics experiments [e.g., @hosoya2005; @kubo2004; @perrillat2013; @perrillat2016; @ledoux2023] reduce uncertainties and the issue of parameter degeneracy has been sufficiently mitigated through independent constraints, seismic observations provide order-of-magnitude estimates of effective $\dot{X}$ and relative rheological strength. Nevertheless, the threshold behavior and scaling relationships from our simulations demonstrate that comprehensive multi-observation approaches are essential for constraining micro-scale processes governing phase transitions and their geodynamic consequences.

# Conclusions {#sec:conclusions}

The olivine $\Leftrightarrow$ wadsleyite transition and resulting 410 structure are strongly influenced by the combined effects of reaction kinetics and rheological strength on flow dynamics. We quantified these effects by integrating an interface-controlled olivine $\Leftrightarrow$ wadsleyite growth model with compressible simulations of mantle plumes and slabs, systematically exploring kinetic factors spanning seven orders of magnitude. Each simulation was evaluated across a large range of viscosity contrasts, and 410 structure was determined after 100 Ma.

Our results reveal fundamentally different responses in hot versus cold environments. Plumes produce consistently sharp discontinuities (2--3 km wide) across the entire parameter space, implying that seismic observations beneath hotspots primarily constrain thermal structure near thermodynamic equilibrium rather than kinetic or rheological parameters. Slabs exhibit distinct threshold behavior across three kinetic regimes---quasi-equilibrium, intermediate, and ultra-sluggish---that are further modulated by viscosity contrasts controlling slab geometry and transit time through the phase transition zone.

Widespread slab penetration of the 410 in seismic tomography [@fukao2013] requires effective reaction rates exceeding the ultra-sluggish kinetic regime ($\dot{X}$ $\gtrsim$ 0.07 Ma$^{-1}$). This critical threshold shifts systematically with rheological strength, revealing that modest variations in either reaction rates or viscosity contrasts can produce substantial diversity in observed 410 topography. Uniquely constraining kinetic versus rheological contributions requires combining 410 structural observations with independent constraints from high-resolution tomography and kinematic data (providing descent geometry) and deep seismicity patterns (indicating metastable olivine volumes).

The 410 can therefore serve as a seismological probe of kinetic conditions in cold subduction environments where disequilibrium effects are amplified. Realizing this potential requires determining discontinuity width through targeted seismic studies as well as reducing uncertainties through targeted mineral physics experiments that better quantify nucleation versus growth mechanisms, water and compositional effects on reaction rates, and microstructural evolution during deformation. Integrating such constraints with high-resolution seismic imaging and the forward modeling framework presented here offers a pathway toward understanding how micro-scale kinetic processes govern mantle convection and shape Earth's interior seismic structure.

# Acknowledgements {.unnumbered}

This work was funded by the UKRI NERC Large Grant no. NE/V018477/1. All computations were undertaken on Barkla2, part of the High Performance Computing facilities at the University of Liverpool, who graciously provided expert support. We thank the Computational Infrastructure for Geodynamics ([https://geodynamics.org](https://geodynamics.org)) which is funded by the National Science Foundation under award EAR-0949446 and EAR-1550901 for supporting the development of ASPECT. We also extend our gratitude towards Sujoy Ghosh for his editorial handling, and two anonymous reviewers for their constructive feedback that improved the manuscript.

# Data Availability {.unnumbered #sec:data-availability}

All data, code, and relevant information for reproducing this work are archived on the OSF [@kerswell2026a] and Zenodo [@kerswell2026b] repositories. All code within these repositories is MIT Licensed and free for use and distribution (see license details). ASPECT version 3.0.0 [@aspect-doi-v3.0.0] was used for the computations in this study and is freely available under the GPL v2.0 or later license.

# Conflict of Interest {.unnumbered #sec:conflict-of-interest}

The authors declare there are no conflicts of interest for this manuscript.

\clearpage

# References {.unnumbered #sec:references}

::: {#refs}
:::
