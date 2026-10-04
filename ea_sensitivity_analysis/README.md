# Standalone ROE parameter-sensitivity analysis

This module implements the parameter-sensitivity workflow for the
`MC_EA.m` region-of-expanse (ROE) model without changing or calling the
existing simulation files.

## Method

- Paired Monte Carlo scenarios are shared by all beam types, divergence
  levels, and factor levels.
- The ROE scan cycle contains the lower boundary circle, an equidistant
  spherical spiral, the upper boundary circle, and the reverse spiral.
- Line and annular beams use the manuscript baseline width `wd = 0.8 mrad`.
- The physical model includes finite-target sampling, target photons,
  atmospheric backscatter, sky background, dark counts, and optimal
  multi-pulse SNR accumulation.
- The primary response is
  `kappa = minimum annular wD / minimum line wD` at a common time-limited
  detection-probability requirement. `kappa < 1` indicates a smaller
  off-axis field requirement for the annular beam.
- Endpoint screening covers nine factors. The two largest absolute effects
  are refined at five levels and evaluated on a 3-by-3 interaction grid.

Both screening and refinement default to `N = 3000`. No ten-thousand-level
run is defined or started by this module.

## Run

```matlab
addpath("D:\lzx\MatlabCode\codexagent_MC\ea_sensitivity_analysis")

% Complete screening, five-level refinement, interaction analysis, and plot
outputs = run_ea_sensitivity_analysis;

% Endpoint screening only
outputs = run_ea_sensitivity_analysis("Mode", "screen");

% Very small end-to-end verification run
outputs = run_ea_sensitivity_analysis( ...
    "Mode", "smoke", "UseParallel", false, "Force", true);

% Rebuild the figure from existing result tables
run_ea_sensitivity_analysis("Mode", "plotonly");
```

The workflow writes resumable raw checkpoints and CSV/MAT summaries under
`results`, and writes the paper-oriented PNG/FIG output under `figures`.
Set `Force=true` to replace compatible checkpoints.

## Main outputs

- `ea_sensitivity_screen_summary.csv`
- `ea_sensitivity_screen_timely.csv`
- `ea_sensitivity_screen_iso_performance.csv`
- `ea_sensitivity_screen_effects.csv`
- `ea_sensitivity_refined_iso_performance.csv`
- `ea_sensitivity_interaction_iso_performance.csv`
- `selected_factors.csv`
- `figures/ea_sensitivity_main.png`
- `figures/ea_sensitivity_main.fig`

## Tests

```matlab
results = runtests(fullfile( ...
    "ea_sensitivity_analysis", "tests", "eaSensitivityTest.m"));
```
