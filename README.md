# HighOrderRoeWENOSolver

A MATLAB implementation of high-order finite-volume solvers for the one-dimensional compressible Euler equations. This project supports multiple Riemann solvers (Roe, Roe-Pike, HLL, HLLC, Osher), high-order WENO/MUSCL reconstruction, and strong stability preserving Runge-Kutta (SSPRK) time integration.

## Features

- Multiple **Riemann solvers**: Roe, Roe-Pike, HLL, HLLC, Osher
- High-order **reconstruction schemes**: WENO3, WENO5, MUSCL
- **SSPRK** time integration (e.g., SSPRK5-6)
- Configurable **flux limiters** (for MUSCL schemes)
- Optional plotting of additional primitive variables (Mach, entropy, etc.)
- Built-in verification driver using the method of manufactured solutions (MMS)

## Project Structure

```
/ (project root)
  main.m                - entrypoint: configure and run solver
  README.md             - this documentation
  src/                  - implementation sources
    runSolver.m         - core loop + configuration driver
    solvers/            - Riemann solvers + flux evaluation
    reconstruction/     - WENO / MUSCL reconstruction routines
    timeintegrators/    - time integrators (RK, SSPRK)
    utils/              - helpers (naming, paths, time step, boundaries)
    plotting/           - plot configuration and plotting helpers
    verification/       - MMS / error norm verification tools
```

## Requirements

- MATLAB (R2018b or newer is recommended)

## Quick Start

From the project root, run:

```matlab
main
```

This will run the solver with the default configuration in `main.m` and print final L2 and Linf error norms.

## Configuration Overview

All solver options are configured in `main.m` by setting values on the `config` struct and then calling `runSolver(config)`.

### Test Cases (`test_num`)

The following test cases are supported. The configuration values map to known shock tube test problems (see `src/utils/testCase.m`):

1. **Toro Test 1 (Modified Sod shock tube)**
   - Reference: Toro, "Riemann Solvers and Numerical Methods for Fluid Dynamics", Chapter 4.
2. **Toro Test 2 (123 Problem)**
   - Reference: Toro, Chapter 4.
3. **Toro Test 3 (Left Blast Wave)**
   - Reference: Toro, Chapter 4.
4. **Toro Test 4 (Right Blast Wave)**
   - Reference: Toro, Chapter 4.
5. **Toro Test 5 (Resulting shock from Toro Test 4 & 5)**
   - Reference: Toro, Chapter 4.
6. **Sod's Shock Tube**
   - Reference: G. A. Sod, JCP 27 (1978).
7. **Left expansion + right strong shock**
8. **Right expansion + left strong shock**
9. **Double shock**
10. **Double expansion**
11. **Cavitation**
12. **Sod shock tube (JCP 27:1, 1978)**
13. **Lax test case**
    - Reference: M. Arora and P. L. Roe, JCP 132:3-11 (1997).
14. **Mach 3 test case**
    - Reference: M. Arora and P. L. Roe, JCP 132:3-11 (1997).
15. **Shock tube with supersonic zone**
16. **Stationary shock**
17. **Right side of 2D Riemann case (case 14)**
18. **Right side of 2D Riemann case (case 15)**
19. **User-specified test**
20. **Manufactured Solution (MMS)**

### Riemann Solvers (`riemann_solver`)

Select one or more of:
- `'Roe'`
- `'Roe-Pike'`
- `'HLL'`
- `'HLLC'`
- `'Osher'`

### Reconstruction (`recon`)

Supported reconstruction schemes:
- `'WENO5'` (5th order WENO)
- `'WENO3'` (3rd order WENO)
- `'MUSCL'` (requires additional MUSCL options)

### MUSCL Options (when `recon = {'MUSCL'}`)

- `flux_limiter` - numeric code for limiter (see mapping below)
- `lambda` - scaling factor for MUSCL limiter
- `kappa` - MUSCL kappa parameter

#### Flux limiter mapping

The limiter is selected via a numeric ID passed in `flux_limiter`. The code uses `src/reconstruction/limitFlux.m`.

- **1** - Minmod
- **2** - van Leer
- **3** - Barth-Jesperson
- **4** - Superbee
- **5** - van Albada 2 (not strictly 2nd-order TVD)
- **6** - van Albada 1
- **7** - CHARM (not 2nd-order TVD)
- **8** - HCUS (not 2nd-order TVD)
- **9** - HQUICK (not 2nd-order TVD)
- **10** - Koren
- **11** - Monotonized central (MC)
- **12** - Osher
- **13** - OSPRE
- **14** - SMART
- **15** - Sweby
- **16** - UMIST
- **17** - Generalized Minmod (beta = 2)
- **18** - No limiter (sets limiter value to zero)

### Time Integration (`time_int_method`)

Common options:
- `'SSPRK5-6'` (6-stage, 5th-order SSP Runge-Kutta)
- `'SSPRK6-5'` (5-stage, 4th-order SSPRK)
- `'RK2'`, `'RK3'`, `'RK4'`, `'RK45'`, `'forwardEuler'`

### Spatial Grid

- `cells` - number of grid cells (spatial resolution)
- `xstart`, `xend` - domain boundaries
- `cfl` - CFL number (stability parameter)

## Central Differencing (CD Term)

The solver includes an optional *central differencing* dissipation term to improve stability and shock resolution. This is controlled by `config.CD_Term_Order`:

- `1` — 2nd-order central differencing (minimal additional dissipation)
- `2` — 4th-order central differencing
- `4` — 6th-order central differencing
- `6` — 8th-order central differencing (strongest dissipation)

The central-differencing term is applied inside the Riemann flux evaluation (see `src/solvers/roeSolver.m`).

### Plotting Options

- `extra_prim_var_plot` - `true`/`false` to plot additional primitive variables such as Mach number and entropy.
- `plot_in_time` - `true`/`false` to plot at every time step (slower).

## Verification / Error Norms

At the end of a run, the solver computes error norms for the last simulation configuration:

- L2 norm: $||u - u_{exact}||_2$
- Linf norm: $||u - u_{exact}||_{\infty}$

These values are printed in the MATLAB command window.

## Example: Custom Script (Standalone)

Create a new file `runExample.m` next to `main.m` containing:

```matlab
% runExample.m

config = struct();
config.test_num = 3;               % Left Blast Wave
config.riemann_solver = {'Roe'};   % Roe solver
config.recon = {'WENO5'};         % 5th order WENO
config.flux_limiter = 1;          % (only needed for MUSCL)
config.lambda = 1;                % (only needed for MUSCL)
config.kappa = 1;                 % (only needed for MUSCL)
config.time_int_method = {'SSPRK5-6'};
config.cells = 200;
config.xstart = 0;
config.xend = 1;
config.cfl = 0.8;
config.extra_prim_var_plot = true;
config.plot_in_time = false;
config.CD_Term_Order = [1];

results = runSolver(config);

last = results(end);
L2 = L2Norm(last.U(:,last.xidx), last.Uexact, config.cells);
Linf = LinfNorm(last.U(:,last.xidx), last.Uexact);
fprintf('L2 = %g, Linf = %g\n', L2, Linf);
```

Run it from MATLAB:

```matlab
runExample
```

## Notes / Tips

- If you add new solvers or reconstruction methods, follow the existing folder structure under `src/solvers/` and `src/reconstruction/`.
- Use `src/verification/l2VerificationMain.m` to run built-in MMS verification cases.

---