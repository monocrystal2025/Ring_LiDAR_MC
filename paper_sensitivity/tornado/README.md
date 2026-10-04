# EA-mode tornado sensitivity analysis

This module performs a paired, one-at-a-time sensitivity analysis for the
EA/ROE moving-target model. It compares line and annular beams through

```text
kappa = minimum annular wD / minimum line wD
```

at the common requirement `P_D(15 s) >= 60%`. Smaller `kappa` means a
stronger annular angular-width advantage. The module is independent of the
older `ea_sensitivity_analysis` directory and never overwrites its outputs.

## Scientific design

- Pulse energy, reflectivity, target speed, and atmospheric extinction are
  varied one at a time at `0.8:0.1:1.2` of their baselines.
- The 17 unique configurations share the same targets, velocity directions,
  and scan phases. The baseline is simulated only once.
- The formal run uses `N=3000` and 1000 paired bootstrap resamples.
- Minimum width is the first raw probability crossing, interpolated only
  between adjacent simulated widths. No smoothing, cumulative maximum, or
  extrapolation is used.
- The initial grid is `5:5:150 mrad`. Cases are extended to 250 and then
  400 mrad only when the point estimate is unattained or paired-bootstrap
  reachability is below 97.5%.
- Extinction changes keep `alpha/beta = 50 sr`.

## Usage

```matlab
addpath("D:/lzx/MatlabCode/codexagent_MC/paper_sensitivity/tornado")

% Complete resumable formal workflow (long-running)
outputs = run_ea_tornado_sensitivity;

% Individual stages
run_ea_tornado_sensitivity(Mode="simulate");
run_ea_tornado_sensitivity(Mode="analyze");
run_ea_tornado_sensitivity(Mode="plot", Visible="on");

% Small non-scientific end-to-end verification
outputs = run_ea_tornado_sensitivity( ...
    Mode="smoke", OutputDirectory=string(tempname), ...
    UseParallel=false);
```

Compatible checkpoints resume automatically. A configuration or source-code
signature mismatch is rejected; use a new output directory rather than
mixing experiments. Completed files are committed atomically after both
beams finish for one configuration and width.

## Outputs

- `ea_tornado_level_results.csv`: all four factors and five levels.
- `ea_tornado_effects.csv`: ranked endpoint tornado effects.
- `ea_tornado_probability_curves.csv`: traceable raw empirical curves.
- `ea_tornado_analysis.mat`: self-contained analysis and provenance.
- `figures/ea_tornado_sensitivity.png/.pdf/.fig`.
- `methodology_and_results.md`: bilingual method, caption, and summary.
- `verification.json`: completeness and cross-file audit.

The formal workflow requires MATLAB R2025a, `readNPY`, Parallel Computing
Toolbox when `UseParallel=true`, and the EA trajectory files under
`G:/BeamVEC_NEW` by default.

## Tests

```matlab
results = runtests("paper_sensitivity/tornado/tests");
```

The smoke test is tagged `Integration` and `Slow`. Smoke outputs use two
targets and a 0.2 s horizon and are only for software verification; they
must not be used as scientific results.
