# 1D Quasi-Steady Supersonic Nozzle Flow Simulation — MacCormack Method

---

## Table of Contents

1. [Problem Statement](#1-problem-statement)
2. [Nozzle Geometry](#2-nozzle-geometry)
3. [Governing Equations](#3-governing-equations)
4. [Non-Dimensional Form](#4-non-dimensional-form)
5. [Non-Conservative Formulation](#5-non-conservative-formulation)
6. [Conservative Formulation](#6-conservative-formulation)
7. [MacCormack Predictor–Corrector Scheme](#7-maccormack-predictorcorrector-scheme)
8. [Boundary Conditions](#8-boundary-conditions)
9. [Time Step — CFL Constraint](#9-time-step--cfl-constraint)
10. [Grid Independence Study](#10-grid-independence-study)
11. [File Structure](#11-file-structure)
12. [Parameters](#12-parameters)
13. [Results](#13-results)

---

## 1. Problem Statement

Simulate **quasi-1D steady isentropic flow** through a converging-diverging (de Laval) nozzle using the **MacCormack explicit finite difference method**, solved as a time-marching problem to steady state. Both **non-conservative** and **conservative** forms of the governing equations are implemented and compared.

The flow is subsonic in the converging section, reaches Mach 1 at the throat, and becomes supersonic in the diverging section — a classical benchmark for compressible flow solvers.

> All variables are **non-dimensionalised** with respect to reservoir (stagnation) conditions. The non-dimensional form removes unit dependencies and allows direct comparison with isentropic analytical solutions.

---

## 2. Nozzle Geometry

The nozzle cross-sectional area profile is parabolic, defined on $x \in [0, 3]$:

$$A(x) = 1 + 2.2\,(x - 1.5)^2$$

The throat is located at $x = 1.5$ where $A = 1$ (minimum area). The nozzle length is non-dimensional, with $x = 0$ as inlet and $x = 3$ as exit.

```
Area
 |    \               /
 |     \             /
 |      \           /
 |       \_________/    ← throat at x = 1.5, A = 1
 |__________________________ x
 0       1.5             3
```

---

## 3. Governing Equations

The quasi-1D Euler equations for compressible flow through a variable-area duct:

**Continuity:**

$$\frac{\partial \rho}{\partial t} + \frac{1}{A}\frac{\partial (\rho V A)}{\partial x} = 0$$

**Momentum:**

$$\frac{\partial V}{\partial t} + V\frac{\partial V}{\partial x} + \frac{1}{\gamma}\left(\frac{\partial T}{\partial x} + \frac{T}{\rho}\frac{\partial \rho}{\partial x}\right) = 0$$

**Energy:**

$$\frac{\partial T}{\partial t} + V\frac{\partial T}{\partial x} + (\gamma - 1)T\left(\frac{\partial V}{\partial x} + V\frac{\partial \ln A}{\partial x}\right) = 0$$

where $\rho$, $V$, $T$, $p = \rho T$ are non-dimensional density, velocity, temperature, and pressure respectively, and $\gamma = 1.4$.

---

## 4. Non-Dimensional Form

All flow variables are normalised by stagnation (reservoir) conditions:

| Variable | Non-dimensional form |
|---|---|
| Density | $\rho / \rho_0$ |
| Velocity | $V / a_0$ where $a_0 = \sqrt{\gamma R T_0}$ |
| Temperature | $T / T_0$ |
| Pressure | $p / p_0 = \rho T$ |
| Speed of sound | $a = \sqrt{T}$ (non-dim) |
| Mach number | $M = V / \sqrt{T}$ |
| Mass flow rate | $\dot{m} = \rho A V$ |

---

## 5. Non-Conservative Formulation

The governing equations are written directly in terms of primitive variables $(\rho, V, T)$. Expanding the divergence terms using the chain rule:

**Continuity:**

$$\frac{\partial \rho}{\partial t} = -\rho \frac{\partial V}{\partial x} - \rho V \frac{\partial \ln A}{\partial x} - V \frac{\partial \rho}{\partial x}$$

**Momentum:**

$$\frac{\partial V}{\partial t} = -V \frac{\partial V}{\partial x} - \frac{1}{\gamma}\left(\frac{\partial T}{\partial x} + \frac{T}{\rho}\frac{\partial \rho}{\partial x}\right)$$

**Energy:**

$$\frac{\partial T}{\partial t} = -V \frac{\partial T}{\partial x} - (\gamma - 1)\,T\left(\frac{\partial V}{\partial x} + V \frac{\partial \ln A}{\partial x}\right)$$

All spatial derivatives are evaluated using **forward differences** in the predictor step and **backward differences** in the corrector step.

**Implemented in:** `non_cons_nozzle_fresh.m`

---

## 6. Conservative Formulation

The equations are recast in terms of solution vectors $\mathbf{U}$ and flux vectors $\mathbf{F}$:

$$\frac{\partial \mathbf{U}}{\partial t} = -\frac{\partial \mathbf{F}}{\partial x} + \mathbf{J}$$

**Solution vectors:**

$$U_1 = \rho A, \quad U_2 = \rho A V, \quad U_3 = \rho A\left(\frac{T}{\gamma - 1} + \frac{\gamma}{2}V^2\right)$$

**Flux vectors:**

$$F_1 = U_2$$

$$F_2 = \frac{U_2^2}{U_1} + \frac{\gamma - 1}{\gamma}\left(U_3 - \frac{\gamma}{2}\frac{U_2^2}{U_1}\right)$$

$$F_3 = \frac{\gamma U_2 U_3}{U_1} - \frac{\gamma(\gamma-1)}{2}\frac{U_2^3}{U_1^2}$$

**Source term** (due to area variation — only in momentum equation):

$$J_2 = \frac{1}{\gamma}\,\rho\,T\,\frac{\partial A}{\partial x}$$

After each time step, primitive variables are recovered from the solution vectors:

$$\rho = \frac{U_1}{A}, \quad V = \frac{U_2}{U_1}, \quad T = (\gamma - 1)\left(\frac{U_3}{U_1} - \frac{\gamma}{2}V^2\right), \quad p = \rho T$$

> The conservative form conserves mass, momentum, and energy exactly in the discrete sense, making it more robust for capturing shocks. The non-conservative form can introduce errors at discontinuities due to the chain-rule expansion.

**Implemented in:** `cons_nozzle_fresh.m`

---

## 7. MacCormack Predictor–Corrector Scheme

The MacCormack method is a **2nd-order accurate explicit finite difference scheme** in both space and time. It uses a two-step predictor–corrector structure with alternating upwind/downwind differencing.

### Predictor Step (forward difference)

Compute predicted values $\bar{\phi}$ using forward spatial differences:

$$\bar{\phi}_j = \phi_j^n + \Delta t \left.\frac{\partial \phi}{\partial t}\right|_j^n$$

where spatial derivatives use $(\phi_{j+1} - \phi_j)/\Delta x$.

### Corrector Step (backward difference)

Compute corrector time derivatives using backward spatial differences on the predicted values $\bar{\phi}$:

$$\left.\frac{\partial \phi}{\partial t}\right|_j^{\text{corr}} = f\left(\bar{\phi}_j, \bar{\phi}_{j-1}\right)$$

using $(\bar{\phi}_j - \bar{\phi}_{j-1})/\Delta x$.

### Average and Update

$$\left.\frac{\partial \phi}{\partial t}\right|_j^{\text{avg}} = \frac{1}{2}\left(\left.\frac{\partial \phi}{\partial t}\right|_j^{\text{pred}} + \left.\frac{\partial \phi}{\partial t}\right|_j^{\text{corr}}\right)$$

$$\phi_j^{n+1} = \phi_j^n + \Delta t \left.\frac{\partial \phi}{\partial t}\right|_j^{\text{avg}}$$

This averaging gives 2nd-order accuracy in time and cancels the leading-order truncation error from the one-sided spatial differences.

**For the conservative form**, the same structure applies to $U_1$, $U_2$, $U_3$ — the flux vectors $F_1$, $F_2$, $F_3$ are updated after the predictor step before computing corrector fluxes.

---

## 8. Boundary Conditions

### Non-Conservative Form

| Boundary | $\rho$ | $V$ | $T$ |
|---|---|---|---|
| Inlet ($j=1$) | Fixed: $\rho = 1$ | Floating: $V_1 = 2V_2 - V_3$ | Fixed: $T = 1$ |
| Outlet ($j=n$) | Floating: $\rho_n = 2\rho_{n-1} - \rho_{n-2}$ | Floating: $V_n = 2V_{n-1} - V_{n-2}$ | Floating: $T_n = 2T_{n-1} - T_{n-2}$ |

### Conservative Form

| Boundary | $U_1$ | $U_2$ | $U_3$ |
|---|---|---|---|
| Inlet ($j=1$) | Fixed: $U_1 = \rho(1) \cdot A(1)$ | Floating: $U_2(1) = 2U_2(2) - U_2(3)$ | Recomputed from $T$, $V$ |
| Outlet ($j=n$) | Floating: $2U_{n-1} - U_{n-2}$ | Floating | Floating |

> **Floating** (linear extrapolation) means the boundary value is linearly extrapolated from the two nearest interior points — no physical condition is imposed, and the boundary evolves freely with the interior solution. **Fixed** values enforce stagnation conditions at the inlet reservoir.

---

## 9. Time Step — CFL Constraint

The time step is determined by the CFL (Courant–Friedrichs–Lewy) condition applied at every grid point, and the minimum is taken:

$$\Delta t = C \cdot \min_j \left(\frac{\Delta x}{a_j + V_j}\right)$$

where $a_j = \sqrt{T_j}$ is the local non-dimensional speed of sound and $C = 0.5$ is the Courant number.

> The denominator $a + V$ represents the fastest wave speed (acoustic + convective) at each node. Using the global minimum ensures stability across all grid points simultaneously.

---

## 10. Grid Independence Study

Three grid resolutions are tested to assess convergence:

| Grid | Nodes | $\Delta x$ |
|---|---|---|
| Coarse | $n = 31$ | $0.1$ |
| Medium | $n = 61$ | $0.05$ |
| Fine | $n = 91$ | $0.033$ |

For each grid, both conservative and non-conservative forms are run for 1400 time steps. Throat values of $\rho$, $V$, $T$, $p$, $M$, and $\dot{m}$ at the final time step are compared via bar plots across all three grids and both formulations.

**Key observation:** The conservative form shows better consistency of mass flow rate across the nozzle length (flat $\dot{m}$ vs. $x$ profile) compared to the non-conservative form, which exhibits slight non-uniformity due to the non-conservative discretisation of the area source term.

**Implemented in:** `challenge_7_final.m` (outer nodal loop), `br_plots.m` (bar plot comparison), `plot_compare.m` (side-by-side flow field comparison)

---

## 11. File Structure

| File | Role | Called by |
|---|---|---|
| `challenge_7_final.m` | Main driver: nodal loop, calls NC and C solvers, comparison plots | — |
| `non_cons_nozzle_fresh.m` | Non-conservative MacCormack solver (primitive variables) | `challenge_7_final.m` |
| `cons_nozzle_fresh.m` | Conservative MacCormack solver (flux vectors $U_1$, $U_2$, $U_3$) | `challenge_7_final.m` |
| `plot_compare.m` | Side-by-side flow field comparison (NC vs C) for all 3 grids | `challenge_7_final.m` |
| `plot_nc.m` | Individual flow field + throat history plots — non-conservative | standalone |
| `plot_c.m` | Individual flow field + throat history plots — conservative | standalone |
| `plot_compare_nc.m` | Standalone NC comparison plot function | standalone |
| `plot_compare_c.m` | Standalone C comparison plot function | standalone |
| `br_plots.m` | Bar plots: throat variables vs. grid size for NC and C | `challenge_7_final.m` |
| `challenge_7.m` | Early development script — non-conservative, animated | — |
| `CHALLENGE_7_practice.m` | Practice script — NC form, throat density history | — |
| `CHALLENGE_7_practice2.m` | Practice script — NC form, full plots + MFR | — |
| `CHALLENGE_7_practice3.m` | Practice script — Conservative form, inline (non-modular) | — |
| `CHALLENGE_7_practice21.m` | Practice script — NC function call wrapper | — |
| `CHALLENGE_7_practice_fresh.m` | Practice version of final driver | — |
| `non_cons_nozzle.m` | Early non-modular NC function (legacy) | — |

### Call Graph (Final Solver)

```
challenge_7_final.m
│
└── Nodal loop [n = 31, 61, 91]:
    ├── non_cons_nozzle_fresh.m    ← NC MacCormack (ρ, V, T)
    │   └── Predictor → BC → Corrector → Average → Update
    │
    └── cons_nozzle_fresh.m        ← C MacCormack (U1, U2, U3)
        └── Predictor → BC → Corrector → Average → Update → Recover primitives
│
├── plot_compare.m                 ← NC vs C comparison (all 3 grids)
└── br_plots.m                     ← Bar plots: throat values vs grid size
```

---

## 12. Parameters

| Parameter | Symbol | Value | Description |
|---|---|---|---|
| Specific heat ratio | $\gamma$ | `1.4` | Air (calorically perfect gas) |
| Nozzle length | — | $x \in [0, 3]$ | Non-dimensional |
| Throat location | $x_t$ | `1.5` | At $A = 1$ |
| Grid nodes | $n$ | `31, 61, 91` | Grid independence study |
| Time steps | $N_t$ | `1400` | Per simulation |
| CFL number | $C$ | `0.5` | Time step safety factor |
| Stagnation density | $\rho_0$ | `1` (non-dim) | Fixed inlet BC |
| Stagnation temperature | $T_0$ | `1` (non-dim) | Fixed inlet BC |

### Initial Conditions

**Non-conservative form** (smooth linear profile):

$$\rho(x,0) = 1 - 0.3146x, \quad T(x,0) = 1 - 0.2314x, \quad V(x,0) = (0.1 + 1.09x)\sqrt{T}$$

**Conservative form** (piecewise profile matching stagnation conditions):

$$\rho = \begin{cases} 1 & 0 \leq x \leq 0.5 \\ 1 - 0.366(x-0.5) & 0.5 \leq x \leq 1.5 \\ 0.634 - 0.3879(x-1.5) & x > 1.5 \end{cases}$$

Same piecewise structure applies to $T$; $V$ is initialised from $V = 0.59/(\rho A)$.

---

## 13. Results

### Mass Flow Rate Evolution — Non-Conservative Form

| n = 31 | n = 61 | n = 91 |
|:---:|:---:|:---:|
| ![MFR NC n=31](Results/1.png) | ![MFR NC n=61](Results/2.png) | ![MFR NC n=91](Results/3.png) |

### Mass Flow Rate Evolution — Conservative Form

| n = 31 | n = 61 | n = 91 |
|:---:|:---:|:---:|
| ![MFR C n=31](Results/4.png) | ![MFR C n=61](Results/5.png) | ![MFR C n=91](Results/6.png) |

### NC vs C — Flow Field Comparison (Steady State)

| n = 31 | n = 61 | n = 91 |
|:---:|:---:|:---:|
| ![Flow field n=31](Results/7.png) | ![Flow field n=61](Results/9.png) | ![Flow field n=91](Results/11.png) |

### NC vs C — Throat Variable History vs Time Steps

| n = 31 | n = 61 | n = 91 |
|:---:|:---:|:---:|
| ![Throat history n=31](Results/8.png) | ![Throat history n=61](Results/10.png) | ![Throat history n=91](Results/12.png) |

### NC vs C — Normalised Mass Flow Rate (Final Time Step)

| n = 31 | n = 61 | n = 91 |
|:---:|:---:|:---:|
| ![MFR normalised n=31](Results/13.png) | ![MFR normalised n=61](Results/14.png) | ![MFR normalised n=91](Results/15.png) |

### Grid Independence — Bar Plot Comparison at Throat

| Pressure | Velocity |
|:---:|:---:|
| ![Pressure at throat](Results/16.png) | ![Velocity at throat](Results/17.png) |

| Temperature | Mach Number |
|:---:|:---:|
| ![Temperature at throat](Results/18.png) | ![Mach number at throat](Results/19.png) |

| Mass Flow Rate | Density |
|:---:|:---:|
| ![MFR at throat](Results/20.png) | ![Density at throat](Results/21.png) |

> The conservative form produces a flatter mass flow rate profile across the nozzle — a direct consequence of the flux-conservative discretisation which preserves mass conservation in the discrete sense. The non-conservative form exhibits slight variation in $\dot{m}(x)$ due to the source term treatment. Both forms converge to the same steady-state throat Mach number $M \approx 1$ as expected for fully choked flow.

---

### Reference

> Anderson, J. D. (1995). *Computational Fluid Dynamics: The Basics with Applications.* McGraw-Hill. (Chapter 7 — Quasi-1D Nozzle Flow)

---

*Solver: 1D Quasi-Steady Supersonic Nozzle — MacCormack Explicit Predictor–Corrector, Conservative and Non-Conservative Forms, Grid Independence Study. Implemented in MATLAB.*
