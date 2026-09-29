# Numerical-Modelling-of-MSMA-Nanoparticles
This repository contains the numerical modelling for MSMA nanoparticle prediction through micro-EDM.
# Micro-EDM Nanoparticle Generation Model

A MATLAB/Octave implementation of a coupled thermal → nucleation → coagulation
→ cooling-rate model for predicting nanoparticle formation during
micro-electro-discharge machining (micro-EDM).

The model framework follows:

> Mitra, S., Muralidhara, Vasa, N. J., & Singaperumal, M. (2010).
> *Investigation on particle generation by micro-electro discharge
> machining.* Proc. SPIE 7590, 75900J.

Given a discharge current, voltage, frequency, and duty cycle, the model
simulates a single discharge pulse's vapor cloud from the moment of
vaporization through cooling, nucleation of stable particle embryos, and
their subsequent growth by Brownian coagulation — returning the resulting
average particle diameter and critical nucleation diameter as functions of
time.

---

## Repository contents

| File | Purpose |
|---|---|
| `microEDM_nanoparticle_model.m` | The model itself, packaged as a single MATLAB/Octave **function**. Contains the thermal model, nucleation-coagulation solver, and cooling-rate model, plus all helper routines as local subfunctions. |
| `test_cases.m` | Example driver script: runs the model for several parameter sets (e.g. a sweep over discharge frequency) and overlays all results on one pair of comparison plots. |
| `compare_avg.png` | Example output — average particle diameter vs. time, overlaid for multiple discharge frequencies. |
| `compare_crit.png` | Example output — critical nucleation diameter vs. time, overlaid for multiple discharge frequencies. |
| `README.md` | This file. |

There are no other dependencies — everything the model needs is contained in
`microEDM_nanoparticle_model.m`.

---

## Requirements

- MATLAB (no additional toolboxes required — only base functions such as
  `erfc` and `accumarray` are used), **or**
- GNU Octave (tested on Octave 8.x)

Both files must be in the same folder (or otherwise both on the MATLAB/Octave
path) for `test_cases.m` to find the model function.

---

## Quick start

```matlab
% Run once with a specific set of process parameters:
%   microEDM_nanoparticle_model(I_curr, V, freq, duty)
%     I_curr : average discharge current   [A]
%     V      : discharge voltage           [V]
%     freq   : discharge frequency         [Hz]
%     duty   : duty cycle, as a FRACTION (0-1), e.g. 0.3 = 30%

results = microEDM_nanoparticle_model(3, 40, 6000, 0.3);

plot(results.t*1e6, results.d_avg);
xlabel('Time (\mus)'); ylabel('Average particle diameter (nm)');
```

Running `test_cases.m` as-is reproduces `compare_avg.png` and
`compare_crit.png`: a sweep over three discharge frequencies (2.5 kHz, 6 kHz,
10 kHz) at fixed current/voltage/duty cycle, all overlaid on a single pair of
plots for direct comparison.

---

## Model overview

The simulation proceeds in four coupled stages for each discharge pulse:

**1. Thermal model.** An empirical correlation (based on discharge current
and pulse-on duration) gives the heat-flux radius at the discharge site.
Solving the transient conduction equation at the centre of this zone gives
the peak temperature reached during the pulse; the average of this peak
temperature and the workpiece material's boiling point is taken as the
initial vapor temperature that seeds everything downstream.

**2. Nucleation model.** As the vapor cools, its actual pressure exceeds the
equilibrium (saturation) vapor pressure at the new, lower temperature — this
excess (the saturation ratio, *S*) drives the vapor to nucleate stable
particle embryos once *S* > 1. Classical nucleation theory gives both the
critical embryo size (below which clusters simply re-evaporate) and the rate
at which stable nuclei are produced.

**3. Coagulation model.** Newly nucleated particles grow further by
colliding and merging with each other under Brownian motion. Particle sizes
are grouped into a geometrically-spaced set of discrete size "nodes"
(sectional / nodal method), and a fixed-pivot mass-conserving scheme tracks
how the population redistributes across nodes as collisions occur.

**4. Cooling-rate model.** Each particle loses heat to the surrounding
dielectric fluid by convection; the resulting cooling rate scales inversely
with particle size, which is why very small, newly-nucleated particles pull
the whole system's temperature down extremely quickly, while coagulation
continues to grow the population's average size over a longer timescale
after the temperature has largely settled.

These four stages are solved together in one explicit time march (an
adaptively-ramped time step: very fine at the start to resolve the fast
initial cooling/nucleation transient, coarsening later once temperature has
stabilized and only slower coagulation remains).

---

## Function reference

### Required inputs

```matlab
results = microEDM_nanoparticle_model(I_curr, V, freq, duty)
```

| Argument | Description | Units |
|---|---|---|
| `I_curr` | Average discharge current | A |
| `V` | Discharge voltage | V |
| `freq` | Discharge frequency | Hz |
| `duty` | Duty cycle, as a **fraction** (0–1), e.g. `0.3` for 30% | – |

The pulse-on duration is derived internally as `Ton = duty / freq`.

### Optional Name-Value inputs

All of the following have defaults already validated for a steel-workpiece
/ deionized-water system; pass any of them to override just that one
setting, e.g. `microEDM_nanoparticle_model(3, 40, 6000, 0.3, 'd_max_nm', 500)`.

**Workpiece / particle material properties**

| Name | Description | Units | Default |
|---|---|---|---|
| `rho_p` | Density | kg/m³ | 8000 |
| `K_th` | Thermal conductivity | W/m·K | 12.5 |
| `alpha_th` | Thermal diffusivity | m²/s | 3.5×10⁻⁶ |
| `Tb` | Boiling point | K | 3250 |
| `Hlv` | Latent heat of vaporization | J/kg | 5.5×10⁶ |
| `Mmolar` | Molar mass | kg/mol | 0.06617 |
| `sigma` | Surface tension of the liquid metal | N/m | 1.8 |
| `Cp_p` | Specific heat capacity | J/kg·K | 550 |
| `Cfrac` | Fraction of discharge energy transferred to workpiece | – | 0.39 |

**Dielectric fluid properties**

| Name | Description | Units | Default |
|---|---|---|---|
| `Tf` | Bulk fluid temperature | K | 300 |
| `mu_f` | Dynamic viscosity | Pa·s | 8.9×10⁻⁴ |
| `kf` | Thermal conductivity | W/m·K | 0.61 |
| `Cp_f` | Specific heat capacity | J/kg·K | 4180 |
| `rho_f` | Density | kg/m³ | 997 |
| `lambda_mfp` | Fluid mean free path | m | 4×10⁻¹⁰ |

**Nodal (sectional) grid**

| Name | Description | Units | Default |
|---|---|---|---|
| `d_min_nm` | Smallest node particle diameter | nm | 1 |
| `d_max_nm` | Largest node particle diameter | nm | 200 |
| `q_ratio` | Node-to-node volume ratio (Hounslow grid) | – | 2 |

> `d_max_nm` should comfortably exceed the largest particle diameter you
> expect the run to reach — check `results.d_avg(end)` against it and raise
> it if the two are close.

**Simulation control**

| Name | Description | Default |
|---|---|---|
| `total_sim_time` | Total simulated time, in seconds | `1/freq` (the inter-pulse period) if left empty |
| `n_coarse_steps` | Number of steps after the time-step has ramped up to its coarse value | 3000 |
| `n_ramp_steps` | Number of steps used to ramp the time step from fine to coarse | 500 |
| `Verbose` | Print thermal-model and summary text to the console | `true` |
| `MakePlot` | Create this run's own pair of figures | `false` |

> `total_sim_time` defaults to one full inter-pulse period (`1/freq`) because
> the model only simulates a single discharge event and does not account for
> interaction between successive pulses — there is no physical basis for
> simulating further than the point where the next pulse would fire. If you
> have a specific reason to look at a different window, override it
> explicitly; otherwise, run once and check whether `results.d_avg` has
> levelled off by the end of the window before trusting a "final" diameter.

### Output

`results` is a struct with the following fields:

| Field | Description |
|---|---|
| `I_curr`, `V`, `freq`, `duty`, `Ton` | The process parameters used for this run |
| `R_um` | Heat-flux radius (μm) |
| `T_peak` | Peak discharge temperature (K) |
| `Ti` | Initial vapor temperature (K) |
| `t` | Time vector (s) |
| `T` | Temperature history (K) |
| `d_avg` | Average particle diameter history (nm) |
| `d_star` | Critical nucleation diameter history (nm); `NaN` once the saturation ratio drops below 1 and nucleation has stopped |
| `label` | A ready-made legend string summarizing this run's parameters |

---

## Running multiple cases and comparing them

`test_cases.m` demonstrates the intended pattern for parameter sweeps: define
a matrix of `[I_curr, V, freq, duty]` rows, loop over them calling the model
with `'MakePlot', false` (so it doesn't pop up its own figures), collect the
returned `results` structs, then plot all of them together using each
struct's `.t`, `.d_avg` / `.d_star`, and `.label` fields. Edit the `cases`
matrix at the top of the script to sweep whatever parameter(s) you're
interested in — frequency, voltage, current, duty cycle, or any combination.

---

## Modelling notes — corrections to the published equations

The source paper is a scanned/OCR'd PDF, and a handful of its printed
equations are not directly usable as typeset (a dimensionally-inconsistent
exponent, an ambiguous repeated symbol, an unstated unit convention). Each
correction — and the reasoning behind it — is documented in full at the top
of `microEDM_nanoparticle_model.m`. In summary:

1. The Clausius–Clapeyron relation (saturation vapor pressure) is implemented
   in its standard, dimensionally-consistent form rather than the paper's
   printed (dimensionally inconsistent, overflow-prone) version.
2. The empirical heat-flux-radius correlation's numeric output is in
   micrometers even though its inputs (current, pulse-on time) are plugged
   in directly in SI units — this matches realistic micro-EDM crater sizes
   and is handled internally.
3. The non-dimensional surface-tension term in the nucleation rate uses the
   standard classical-nucleation-theory form (based on monomer surface
   area), resolving an ambiguous repeated symbol in the original text.
4. The dielectric fluid's mean free path (needed for the cooling-rate model)
   is not derived in the paper and is exposed as an editable parameter
   (`lambda_mfp`) with a reasonable default for water.
5. The Reynolds number used in the Nusselt-number correlation is computed
   with its standard definition, which the paper does not spell out.
6. The finest time step used at the start of the simulation is fixed in
   absolute terms (~1 ps) rather than scaled from `total_sim_time`, to avoid
   numerical overshoot (the temperature jumping straight to the fluid
   temperature in a single step) when a much longer `total_sim_time` is
   requested.

None of these corrections change the modelling framework — only the specific
formulas needed fixing to produce physically meaningful, numerically stable
results.

---

## Choosing physically sensible inputs

A few parameters are easy to set to values that are numerically valid but
physically meaningless, which silently produces misleading results rather
than an error. Worth double-checking before trusting a run:

- **`sigma` (surface tension)** must be a liquid-*metal* value (order
  1–2 N/m for most metals), not a value for water or another dielectric
  fluid (~0.07 N/m). Too low a value collapses the critical nucleus size
  below the size of a single atom, which is not physically meaningful.
- **`d_max_nm`** should comfortably exceed the diameter the run actually
  reaches (check `results.d_avg(end)`), or growth will appear to plateau
  artificially at the grid's ceiling.
- **`total_sim_time`** should be checked against whether `results.d_avg` has
  actually levelled off by the end of the window — Brownian coagulation does
  not have a sharp mathematical plateau, so "run it until it looks flat" is
  not always a well-defined stopping point. The inter-pulse period (`1/freq`,
  the default) is the most physically defensible choice, since the model
  does not track interaction between successive pulses.

---

## Citation

If you use this model, please cite the original paper:

> Mitra, S., Muralidhara, Vasa, N. J., & Singaperumal, M. (2010).
> Investigation on particle generation by micro-electro discharge machining.
> *Proceedings of SPIE*, 7590, 75900J. https://doi.org/10.1117/12.845709
