<!--

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

Note that the deviatoric stress and strain rate tensors are symmetric such that $\sigma^{\prime}_{ij} = \sigma^{\prime}_{ji}$ and $\dot{\epsilon}^{\prime}_{ij} = \dot{\epsilon}^{\prime}_{ji}$, which by definition means that there is no rotational deformation in the fluid flow, only rigid-body rotation. Because of this symmetry, the full matrix form the deviatoric stress tensor can be written as:

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

-->

\clearpage

# Measuring 410 Displacement and Width {.unnumbered #sec:measuring-displacement-width}

Structure of the 410 was evaluated from the volume fraction field $X$ along vertical profiles that intersected either: 1) the widest phase transition zone or 2) the largest displacement of the olivine $\Leftrightarrow$ wadsleyite reaction found within model domain. The phase transition zone width was defined as the difference between the depths at $X$ = 0.9 and $X$ = 0.1, while the olivine $\Leftrightarrow$ wadsleyite reaction displacement was defined as the offset between the nominal equilibrium reaction depth and the depth at $X$ = 0.9. The exact vertical profile position was chosen heuristically by "best fit" criteria based on structural characteristics of the phase transition zone. The maximum reaction rate $\dot{X}$ and vertical velocity $\vec{u}_y$ were evaluated along the selected vertical profile, within the upper and lower bounds of the phase transition zone, after 100 Ma of evolution. Selected examples of vertical profile picking in slab simulations are shown in Figure \ref{fig:slab-composition-set1}.

\clearpage

![Slab simulations with intermediate strength contrasts ($B$ = 4) demonstrating ultra-sluggish (top row: $Z$ = 3.0e0 K s$^{-1}$), intermediate (middle row: $Z$ = 4.7e2 K s$^{-1}$), and quasi-equilibrium (bottom row: $Z$ = 7.0e7 K s$^{-1}$) kinetic regimes after 100 Ma evolution. Panels show dynamic temperature $\hat{T}$ (left column), reaction rate $\dot{X}$ (middle column), and volume fraction of wadsleyite $X$ (right column). Thin black lines indicate the vertical profile and depth bounds for measuring 410 structure (width and displacement).](../figs/simulation/compositions/slab-Z3.0e0-B4-Z4.7e2-B4-Z7.0e7-B4-set1-composition-0010.png){#fig:slab-composition-set1 width=100%}

\clearpage

![Plume simulations with intermediate strength contrasts ($B$ = 4) demonstrating ultra-sluggish (top row: $Z$ = 3.0e0 K s$^{-1}$), intermediate (middle row: $Z$ = 4.7e2 K s$^{-1}$), and quasi-equilibrium (bottom row: $Z$ = 7.0e7 K s$^{-1}$) kinetic regimes after 100 Ma evolution. Panels show dynamic temperature $\hat{T}$ (left column), reaction rate $\dot{X}$ (middle column), and volume fraction of wadsleyite $X$ (right column). Thin black lines indicate the vertical profile and depth bounds for measuring 410 structure (width and displacement).](../figs/simulation/compositions/plume-Z3.0e0-B4-Z4.7e2-B4-Z7.0e7-B4-set1-composition-0010.png){#fig:plume-composition-set1 width=100%}

\clearpage

| Simulation   | $Z$   |   $B$ |   Displacement |   Width |   Max $\vec{u}_y$ |   log$_{10}$ Max $\dot{X}$ |
|:-------------|------:|------:|---------------:|--------:|------------------:|---------------------------:|
| plume        | 3.0e0 |     2 |             20 |      19 |              1.71 |                      -0.01 |
| plume        | 3.0e0 |     4 |             20 |      19 |              1.78 |                      -0.01 |
| plume        | 3.0e0 |     6 |             19 |      20 |              1.85 |                      -0.01 |
| plume        | 3.0e0 |     8 |             19 |      21 |              1.91 |                      -0.01 |
| plume        | 3.0e0 |    10 |             18 |      21 |              1.96 |                      -0.01 |
| plume        | 7.0e0 |     2 |             21 |      14 |              2.06 |                       0.2  |
| plume        | 7.0e0 |     4 |             21 |      15 |              2.19 |                       0.21 |
| plume        | 7.0e0 |     6 |             20 |      15 |              2.32 |                       0.22 |
| plume        | 7.0e0 |     8 |             20 |      16 |              2.45 |                       0.22 |
| plume        | 7.0e0 |    10 |             19 |      17 |              2.58 |                       0.22 |
| plume        | 1.6e1 |     2 |             22 |      10 |              2.29 |                       0.39 |
| plume        | 1.6e1 |     4 |             22 |      11 |              2.46 |                       0.41 |
| plume        | 1.6e1 |     6 |             21 |      12 |              2.65 |                       0.42 |
| plume        | 1.6e1 |     8 |             21 |      12 |              2.84 |                       0.42 |
| plume        | 1.6e1 |    10 |             20 |      12 |              3.04 |                       0.44 |
| plume        | 3.7e1 |     2 |             23 |       8 |              2.42 |                       0.59 |
| plume        | 3.7e1 |     4 |             23 |       8 |              2.62 |                       0.61 |
| plume        | 3.7e1 |     6 |             22 |       8 |              2.85 |                       0.61 |
| plume        | 3.7e1 |     8 |             22 |       9 |              3.09 |                       0.61 |
| plume        | 3.7e1 |    10 |             21 |       9 |              3.34 |                       0.64 |
| plume        | 8.7e1 |     2 |             24 |       6 |              2.53 |                       0.76 |
| plume        | 8.7e1 |     4 |             24 |       6 |              2.76 |                       0.75 |
| plume        | 8.7e1 |     6 |             23 |       7 |              3    |                       0.8  |
| plume        | 8.7e1 |     8 |             23 |       7 |              3.29 |                       0.84 |
| plume        | 8.7e1 |    10 |             22 |       7 |              3.58 |                       0.86 |
| plume        | 2.0e2 |     2 |             24 |       4 |              2.57 |                       0.97 |
| plume        | 2.0e2 |     4 |             24 |       5 |              2.84 |                       0.94 |
| plume        | 2.0e2 |     6 |             24 |       6 |              3.07 |                       0.89 |
| plume        | 2.0e2 |     8 |             24 |       6 |              3.38 |                       0.98 |
| plume        | 2.0e2 |    10 |             23 |       6 |              3.72 |                       1.04 |
| plume        | 4.7e2 |     2 |             24 |       3 |              2.59 |                       1.19 |
| plume        | 4.7e2 |     4 |             24 |       4 |              2.85 |                       1.13 |
| plume        | 4.7e2 |     6 |             24 |       4 |              3.14 |                       1.07 |
| plume        | 4.7e2 |     8 |             24 |       4 |              3.48 |                       1.04 |
| plume        | 4.7e2 |    10 |             24 |       5 |              3.74 |                       1.1  |
| plume        | 1.1e3 |     2 |             24 |       3 |              2.65 |                       1.42 |
| plume        | 1.1e3 |     4 |             24 |       3 |              2.87 |                       1.3  |
| plume        | 1.1e3 |     6 |             24 |       3 |              3.15 |                       1.26 |
| plume        | 1.1e3 |     8 |             24 |       4 |              3.51 |                       1.19 |
| plume        | 1.1e3 |    10 |             24 |       4 |              3.86 |                       1.14 |
| plume        | 2.6e3 |     2 |             24 |       3 |              2.79 |                       1.6  |
| plume        | 2.6e3 |     4 |             24 |       3 |              2.91 |                       1.5  |
| plume        | 2.6e3 |     6 |             24 |       3 |              3.23 |                       1.44 |
| plume        | 2.6e3 |     8 |             24 |       3 |              3.46 |                       1.5  |
| plume        | 2.6e3 |    10 |             24 |       4 |              3.92 |                       1.27 |
| plume        | 6.0e3 |     2 |             24 |       2 |              2.5  |                       1.63 |
| plume        | 6.0e3 |     4 |             24 |       3 |              2.88 |                       1.93 |
| plume        | 6.0e3 |     6 |             24 |       3 |              3.29 |                       1.65 |
| plume        | 6.0e3 |     8 |             24 |       3 |              3.47 |                       1.86 |
| plume        | 6.0e3 |    10 |             24 |       4 |              3.88 |                       1.63 |
| plume        | 1.4e4 |     2 |             24 |       2 |              2.51 |                       1.7  |
| plume        | 1.4e4 |     4 |             24 |       3 |              2.82 |                       2.03 |
| plume        | 1.4e4 |     6 |             24 |       2 |              3.09 |                       1.68 |
| plume        | 1.4e4 |     8 |             24 |       3 |              3.55 |                       2.11 |
| plume        | 1.4e4 |    10 |             24 |       4 |              3.89 |                       2.02 |
| plume        | 3.3e4 |     2 |             24 |       2 |              2.42 |                       1.82 |
| plume        | 3.3e4 |     4 |             24 |       2 |              2.63 |                       1.66 |
| plume        | 3.3e4 |     6 |             24 |       3 |              2.99 |                       2.3  |
| plume        | 3.3e4 |     8 |             24 |       3 |              3.63 |                       2.36 |
| plume        | 3.3e4 |    10 |             24 |       3 |              4.03 |                       2.21 |
| plume        | 7.8e4 |     2 |             24 |       2 |              2.06 |                       1.64 |
| plume        | 7.8e4 |     4 |             24 |       3 |              2.9  |                       2.67 |
| plume        | 7.8e4 |     6 |             24 |       2 |              2.89 |                       2.06 |
| plume        | 7.8e4 |     8 |             24 |       3 |              3.53 |                       2.9  |
| plume        | 7.8e4 |    10 |             24 |       3 |              5.33 |                       2.44 |
| plume        | 1.8e5 |     2 |             24 |       2 |              2.06 |                       1.96 |
| plume        | 1.8e5 |     4 |             24 |       3 |              2.88 |                       3.08 |
| plume        | 1.8e5 |     6 |             24 |       2 |              2.93 |                       2.52 |
| plume        | 1.8e5 |     8 |             24 |       3 |              3.51 |                       3.28 |
| plume        | 1.8e5 |    10 |             24 |       2 |              4.26 |                       2.65 |
| plume        | 4.3e5 |     2 |             24 |       2 |              2.06 |                       2.34 |
| plume        | 4.3e5 |     4 |             24 |       3 |              2.87 |                       3.48 |
| plume        | 4.3e5 |     6 |             24 |       2 |              3.24 |                       2.93 |
| plume        | 4.3e5 |     8 |             24 |       3 |              3.51 |                       3.66 |
| plume        | 4.3e5 |    10 |             24 |       3 |              4.18 |                       3.13 |
| plume        | 1.0e6 |     2 |             24 |       2 |              2.06 |                       2.71 |
| plume        | 1.0e6 |     4 |             24 |       3 |              2.87 |                       3.85 |
| plume        | 1.0e6 |     6 |             24 |       2 |              3.23 |                       3.31 |
| plume        | 1.0e6 |     8 |             24 |       3 |              3.51 |                       4.03 |
| plume        | 1.0e6 |    10 |             24 |       3 |              4.19 |                       3.48 |
| plume        | 2.4e6 |     2 |             24 |       2 |              2.06 |                       3.09 |
| plume        | 2.4e6 |     4 |             24 |       3 |              2.87 |                       4.24 |
| plume        | 2.4e6 |     6 |             24 |       2 |              3.22 |                       3.7  |
| plume        | 2.4e6 |     8 |             24 |       3 |              3.51 |                       4.41 |
| plume        | 2.4e6 |    10 |             24 |       3 |              4.21 |                       3.83 |
| plume        | 5.6e6 |     2 |             24 |       2 |              2.06 |                       3.46 |
| plume        | 5.6e6 |     4 |             24 |       3 |              2.86 |                       4.61 |
| plume        | 5.6e6 |     6 |             24 |       2 |              3.22 |                       4.08 |
| plume        | 5.6e6 |     8 |             24 |       3 |              3.51 |                       4.77 |
| plume        | 5.6e6 |    10 |             24 |       2 |              4.24 |                       4.17 |
| plume        | 1.3e7 |     2 |             24 |       2 |              2.05 |                       3.83 |
| plume        | 1.3e7 |     4 |             24 |       3 |              2.86 |                       4.98 |
| plume        | 1.3e7 |     6 |             24 |       2 |              3.22 |                       4.45 |
| plume        | 1.3e7 |     8 |             24 |       3 |              3.51 |                       5.14 |
| plume        | 1.3e7 |    10 |             24 |       2 |              4.31 |                       4.45 |
| plume        | 3.0e7 |     2 |             24 |       2 |              2.05 |                       4.19 |
| plume        | 3.0e7 |     4 |             24 |       3 |              2.86 |                       5.34 |
| plume        | 3.0e7 |     6 |             24 |       2 |              3.21 |                       4.82 |
| plume        | 3.0e7 |     8 |             24 |       3 |              3.51 |                       5.5  |
| plume        | 3.0e7 |    10 |             24 |       2 |              3.71 |                       4.05 |
| plume        | 7.0e7 |     2 |             24 |       2 |              2.05 |                       4.56 |
| plume        | 7.0e7 |     4 |             24 |       3 |              2.86 |                       5.71 |
| plume        | 7.0e7 |     6 |             24 |       2 |              3.21 |                       5.18 |
| plume        | 7.0e7 |     8 |             24 |       3 |              3.51 |                       5.87 |
| plume        | 7.0e7 |    10 |             24 |       2 |              3.7  |                       4.4  |
| slab         | 3.0e0 |     2 |             83 |       7 |              0.09 |                      -2    |
| slab         | 3.0e0 |     4 |             85 |       8 |              0.08 |                      -1.92 |
| slab         | 3.0e0 |     6 |             88 |      10 |              0.1  |                      -1.85 |
| slab         | 3.0e0 |     8 |             91 |      14 |              0.13 |                      -1.8  |
| slab         | 3.0e0 |    10 |             95 |      20 |              0.16 |                      -1.74 |
| slab         | 7.0e0 |     2 |             66 |       9 |              0.09 |                      -1.74 |
| slab         | 7.0e0 |     4 |             69 |      11 |              0.11 |                      -1.71 |
| slab         | 7.0e0 |     6 |             73 |      14 |              0.13 |                      -1.66 |
| slab         | 7.0e0 |     8 |             78 |      20 |              0.15 |                      -1.6  |
| slab         | 7.0e0 |    10 |             84 |      32 |              0.23 |                      -1.58 |
| slab         | 1.6e1 |     2 |             50 |      11 |              0.11 |                      -1.56 |
| slab         | 1.6e1 |     4 |             54 |      13 |              0.12 |                      -1.53 |
| slab         | 1.6e1 |     6 |             60 |      18 |              0.16 |                      -1.49 |
| slab         | 1.6e1 |     8 |            144 |     150 |              0.72 |                      -1.53 |
| slab         | 1.6e1 |    10 |            167 |     172 |              0.84 |                      -1.39 |
| slab         | 3.7e1 |     2 |             38 |      13 |              0.13 |                      -1.37 |
| slab         | 3.7e1 |     4 |             45 |      18 |              0.18 |                      -1.37 |
| slab         | 3.7e1 |     6 |            144 |     155 |              1    |                      -1.38 |
| slab         | 3.7e1 |     8 |            169 |     179 |              1.09 |                      -1.23 |
| slab         | 3.7e1 |    10 |            145 |     154 |              0.99 |                      -1.17 |
| slab         | 8.7e1 |     2 |             48 |      32 |              0.29 |                      -1.5  |
| slab         | 8.7e1 |     4 |             60 |      49 |              0.46 |                      -1.46 |
| slab         | 8.7e1 |     6 |            169 |     183 |              1.34 |                      -1.11 |
| slab         | 8.7e1 |     8 |             88 |     101 |              0.91 |                      -0.99 |
| slab         | 8.7e1 |    10 |             57 |      70 |              0.72 |                      -0.97 |
| slab         | 2.0e2 |     2 |            158 |     184 |              1.11 |                      -1.69 |
| slab         | 2.0e2 |     4 |            139 |     158 |              1.09 |                      -1.17 |
| slab         | 2.0e2 |     6 |             48 |      64 |              0.68 |                      -0.83 |
| slab         | 2.0e2 |     8 |             32 |      49 |              0.56 |                      -0.93 |
| slab         | 2.0e2 |    10 |             27 |      43 |              0.57 |                      -0.9  |
| slab         | 4.7e2 |     2 |            143 |     167 |              1.35 |                      -1.21 |
| slab         | 4.7e2 |     4 |             89 |     106 |              1.15 |                      -0.52 |
| slab         | 4.7e2 |     6 |             25 |      45 |              0.54 |                      -0.78 |
| slab         | 4.7e2 |     8 |             14 |      35 |              0.52 |                      -0.81 |
| slab         | 4.7e2 |    10 |              9 |      29 |              0.54 |                      -0.75 |
| slab         | 1.1e3 |     2 |            102 |     124 |              1.54 |                      -0.51 |
| slab         | 1.1e3 |     4 |             56 |      75 |              1.17 |                      -0.32 |
| slab         | 1.1e3 |     6 |             16 |      40 |              0.65 |                      -0.69 |
| slab         | 1.1e3 |     8 |             -1 |      24 |              0.47 |                      -0.73 |
| slab         | 1.1e3 |    10 |             -4 |      19 |              0.5  |                      -0.61 |
| slab         | 2.6e3 |     2 |             64 |      87 |              1.75 |                      -0.32 |
| slab         | 2.6e3 |     4 |             33 |      56 |              1.35 |                      -0.29 |
| slab         | 2.6e3 |     6 |              7 |      33 |              0.92 |                      -0.48 |
| slab         | 2.6e3 |     8 |             -8 |      20 |              0.59 |                      -0.54 |
| slab         | 2.6e3 |    10 |            -15 |      11 |              0.37 |                      -0.45 |
| slab         | 6.0e3 |     2 |             34 |      61 |              2.01 |                      -0.26 |
| slab         | 6.0e3 |     4 |             14 |      41 |              1.59 |                      -0.23 |
| slab         | 6.0e3 |     6 |             -4 |      25 |              1.1  |                      -0.34 |
| slab         | 6.0e3 |     8 |            -18 |      11 |              0.53 |                      -0.27 |
| slab         | 6.0e3 |    10 |            -19 |       8 |              0.42 |                      -0.26 |
| slab         | 1.4e4 |     2 |             12 |      42 |              2.26 |                      -0.17 |
| slab         | 1.4e4 |     4 |             -2 |      29 |              1.76 |                      -0.17 |
| slab         | 1.4e4 |     6 |            -14 |      18 |              1.22 |                      -0.19 |
| slab         | 1.4e4 |     8 |            -21 |       9 |              0.62 |                      -0.11 |
| slab         | 1.4e4 |    10 |            -22 |       7 |              0.52 |                      -0.08 |
| slab         | 3.3e4 |     2 |             -5 |      29 |              2.46 |                      -0.05 |
| slab         | 3.3e4 |     4 |            -14 |      20 |              1.88 |                      -0.02 |
| slab         | 3.3e4 |     6 |            -22 |      11 |              1.17 |                       0.08 |
| slab         | 3.3e4 |     8 |            -25 |       7 |              0.73 |                       0.06 |
| slab         | 3.3e4 |    10 |            -24 |       6 |              0.57 |                       0.08 |
| slab         | 7.8e4 |     2 |            -17 |      20 |              2.59 |                       0.13 |
| slab         | 7.8e4 |     4 |            -24 |      11 |              1.87 |                       0.27 |
| slab         | 7.8e4 |     6 |            -27 |       9 |              1.25 |                       0.21 |
| slab         | 7.8e4 |     8 |            -27 |       5 |              0.82 |                       0.29 |
| slab         | 7.8e4 |    10 |            -27 |       4 |              0.6  |                       0.17 |
| slab         | 1.8e5 |     2 |            -27 |       8 |              2.51 |                       0.53 |
| slab         | 1.8e5 |     4 |            -28 |      10 |              1.93 |                       0.33 |
| slab         | 1.8e5 |     6 |            -30 |       6 |              1.32 |                       0.37 |
| slab         | 1.8e5 |     8 |            -29 |       8 |              0.82 |                       0.57 |
| slab         | 1.8e5 |    10 |            -28 |       3 |              0.62 |                       0.35 |
| slab         | 4.3e5 |     2 |            -31 |       8 |              2.6  |                       0.6  |
| slab         | 4.3e5 |     4 |            -32 |       7 |              1.96 |                       0.52 |
| slab         | 4.3e5 |     6 |            -32 |       4 |              1.35 |                       0.49 |
| slab         | 4.3e5 |     8 |            -30 |       7 |              0.82 |                       0.85 |
| slab         | 4.3e5 |    10 |            -29 |       2 |              0.63 |                       0.51 |
| slab         | 1.0e6 |     2 |            -34 |       6 |              2.66 |                       0.7  |
| slab         | 1.0e6 |     4 |            -35 |       5 |              1.97 |                       0.62 |
| slab         | 1.0e6 |     6 |            -34 |       3 |              1.34 |                       0.78 |
| slab         | 1.0e6 |     8 |            -30 |       6 |              0.88 |                       1.16 |
| slab         | 1.0e6 |    10 |            -29 |       2 |              0.63 |                       0.66 |
| slab         | 2.4e6 |     2 |            -37 |       4 |              2.69 |                       0.85 |
| slab         | 2.4e6 |     4 |            -37 |       3 |              1.98 |                       0.87 |
| slab         | 2.4e6 |     6 |            -34 |       7 |              1.34 |                       1.06 |
| slab         | 2.4e6 |     8 |            -31 |       5 |              0.8  |                       1.26 |
| slab         | 2.4e6 |    10 |            -29 |       2 |              0.62 |                       0.85 |
| slab         | 5.6e6 |     2 |            -38 |       4 |              2.7  |                       1.12 |
| slab         | 5.6e6 |     4 |            -37 |       3 |              1.99 |                       1.07 |
| slab         | 5.6e6 |     6 |            -35 |       3 |              1.31 |                       1.05 |
| slab         | 5.6e6 |     8 |            -32 |       5 |              0.78 |                       1.41 |
| slab         | 5.6e6 |    10 |            -29 |       2 |              0.69 |                       1.13 |
| slab         | 1.3e7 |     2 |            -39 |       7 |              2.71 |                       1.41 |
| slab         | 1.3e7 |     4 |            -38 |       2 |              2    |                       1.32 |
| slab         | 1.3e7 |     6 |            -36 |       1 |              1.31 |                       1.29 |
| slab         | 1.3e7 |     8 |            -32 |       5 |              0.78 |                       1.78 |
| slab         | 1.3e7 |    10 |            -29 |       2 |              0.71 |                       1.23 |
| slab         | 3.0e7 |     2 |            -39 |       7 |              2.72 |                       1.62 |
| slab         | 3.0e7 |     4 |            -38 |       2 |              2.01 |                       1.65 |
| slab         | 3.0e7 |     6 |            -36 |       5 |              1.32 |                       1.85 |
| slab         | 3.0e7 |     8 |            -32 |       5 |              0.79 |                       2.15 |
| slab         | 3.0e7 |    10 |            -29 |       2 |              0.71 |                       1.59 |
| slab         | 7.0e7 |     2 |            -39 |       7 |              2.67 |                       2.26 |
| slab         | 7.0e7 |     4 |            -38 |       2 |              2.01 |                       2.01 |
| slab         | 7.0e7 |     6 |            -36 |       5 |              1.32 |                       2.24 |
| slab         | 7.0e7 |     8 |            -32 |       5 |              0.79 |                       2.52 |
| slab         | 7.0e7 |    10 |            -29 |       1 |              0.56 |                       1.84 |

Table: Summary of the kinetic prefactor $Z$, rheological activation factor $B$, 410 structure (displacement and width), maximum vertical velocity $\vec{u}_y$, and maximum reaction rate $\dot{X}$ evaluated in plume and slab simulations after 100 Ma of evolution. Units are $Z$: K s$^{-1}$, $B$: none, displacement: km, width: km, $\vec{u}_y$: cm/yr, log$_{10}$ $\dot{X}$: log$_{10}$ Ma$^{-1}$. {#tbl:depth-profile-summary}

\clearpage

# Effects of Rheological Strength Contrasts on Flow Dynamics {.unnumbered #sec:rheological-strength-contrasts}

The following simulation snapshots demonstrate the effect of the rheological activation factor $B$ on 410 structure after 100 Ma of evolution.

\clearpage

![Slab simulations with low strength contrasts ($B$ = 2) demonstrating ultra-sluggish (top row: $Z$ = 3.0e0 K s$^{-1}$), intermediate (middle row: $Z$ = 4.7e2 K s$^{-1}$), and quasi-equilibrium (bottom row: $Z$ = 7.0e7 K s$^{-1}$) kinetic regimes after 100 Ma evolution. Panels show dynamic temperature $\hat{T}$ (left column), dynamic density $\hat{\rho}$ (middle column), and pressure-wave velocity $V_p$ (right column).](../figs/simulation/compositions/slab-Z3.0e0-B2-Z3.0e0-B6-Z3.0e0-B10-set2-composition-0010.png){#fig:slab-B2-composition-set2 width=100%}

\clearpage

![Slab simulations with intermediate strength contrasts ($B$ = 4) showing ultra-sluggish (top row: $Z$ = 3.0e0 K s$^{-1}$), intermediate (middle row: $Z$ = 4.7e2 K s$^{-1}$), and quasi-equilibrium (bottom row: $Z$ = 7.0e7 K s$^{-1}$) kinetic regimes after 100 Ma evolution. Panels show dynamic temperature $\hat{T}$ (left column), dynamic density $\hat{\rho}$ (middle column), and pressure-wave velocity $V_p$ (right column).](../figs/simulation/compositions/slab-Z3.0e0-B4-Z4.7e2-B4-Z7.0e7-B4-set2-composition-0010.png){#fig:slab-B4-composition-set2 width=100%}

\clearpage

![Slab simulations with intermediate strength contrasts ($B$ = 6) demonstrating ultra-sluggish (top row: $Z$ = 4.7e2 K s$^{-1}$), intermediate (middle row: $Z$ = 4.7e2 K s$^{-1}$), and quasi-equilibrium (bottom row: $Z$ = 7.0e7 K s$^{-1}$) kinetic regimes after 100 Ma evolution. Panels show dynamic temperature $\hat{T}$ (left column), dynamic density $\hat{\rho}$ (middle column), and pressure-wave velocity $V_p$ (right column).](../figs/simulation/compositions/slab-Z4.7e2-B2-Z4.7e2-B6-Z4.7e2-B10-set2-composition-0010.png){#fig:slab-B6-composition-set2 width=100%}

\clearpage

![Slab simulations with high strength contrasts ($B$ = 10) demonstrating ultra-sluggish (top row: $Z$ = 7.0e7 K s$^{-1}$), intermediate (middle row: $Z$ = 4.7e2 K s$^{-1}$), and quasi-equilibrium (bottom row: $Z$ = 7.0e7 K s$^{-1}$) kinetic regimes after 100 Ma evolution. Panels show dynamic temperature $\hat{T}$ (left column), dynamic density $\hat{\rho}$ (middle column), and pressure-wave velocity $V_p$ (right column).](../figs/simulation/compositions/slab-Z7.0e7-B2-Z7.0e7-B6-Z7.0e7-B10-set2-composition-0010.png){#fig:slab-B10-composition-set2 width=100%}

\clearpage

![Plume simulations with low strength contrasts ($B$ = 2) demonstrating ultra-sluggish (top row: $Z$ = 3.0e0 K s$^{-1}$), intermediate (middle row: $Z$ = 4.7e2 K s$^{-1}$), and quasi-equilibrium (bottom row: $Z$ = 7.0e7 K s$^{-1}$) kinetic regimes after 100 Ma evolution. Panels show dynamic temperature $\hat{T}$ (left column), dynamic density $\hat{\rho}$ (middle column), and pressure-wave velocity $V_p$ (right column).](../figs/simulation/compositions/plume-Z3.0e0-B2-Z3.0e0-B6-Z3.0e0-B10-set2-composition-0010.png){#fig:plume-B2-composition-set2 width=100%}

\clearpage

![Plume simulations with intermediate strength contrasts ($B$ = 4) showing ultra-sluggish (top row: $Z$ = 3.0e0 K s$^{-1}$), intermediate-sluggish (middle row: $Z$ = 4.7e2 K s$^{-1}$), and quasi-equilibrium (bottom row: $Z$ = 7.0e7 K s$^{-1}$) kinetic regimes after 100 Ma evolution. Panels show dynamic temperature $\hat{T}$ (left column), dynamic density $\hat{\rho}$ (middle column), and pressure-wave velocity $V_p$ (right column).](../figs/simulation/compositions/plume-Z3.0e0-B4-Z4.7e2-B4-Z7.0e7-B4-set2-composition-0010.png){#fig:plume-B4-composition-set2 width=100%}

\clearpage

![Plume simulations with intermediate strength contrasts ($B$ = 6) demonstrating ultra-sluggish (top row: $Z$ = 4.7e2 K s$^{-1}$), intermediate (middle row: $Z$ = 4.7e2 K s$^{-1}$), and quasi-equilibrium (bottom row: $Z$ = 7.0e7 K s$^{-1}$) kinetic regimes after 100 Ma evolution. Panels show dynamic temperature $\hat{T}$ (left column), dynamic density $\hat{\rho}$ (middle column), and pressure-wave velocity $V_p$ (right column).](../figs/simulation/compositions/plume-Z4.7e2-B2-Z4.7e2-B6-Z4.7e2-B10-set2-composition-0010.png){#fig:plume-B6-composition-set2 width=100%}

\clearpage

![Plume simulations with high strength contrasts ($B$ = 10) demonstrating ultra-sluggish (top row: $Z$ = 7.0e7 K s$^{-1}$), intermediate (middle row: $Z$ = 4.7e2 K s$^{-1}$), and quasi-equilibrium (bottom row: $Z$ = 7.0e7 K s$^{-1}$) kinetic regimes after 100 Ma evolution. Panels show dynamic temperature $\hat{T}$ (left column), dynamic density $\hat{\rho}$ (middle column), and pressure-wave velocity $V_p$ (right column).](../figs/simulation/compositions/plume-Z7.0e7-B2-Z7.0e7-B6-Z7.0e7-B10-set2-composition-0010.png){#fig:plume-B10-composition-set2 width=100%}

\clearpage

<!--
# References {.unnumbered #sec:references}

::: {#refs}
:::
-->
