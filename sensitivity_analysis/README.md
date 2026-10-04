# Standalone ROI sensitivity analysis

This folder contains a self-contained, paired Monte Carlo sensitivity study
for the spot, line, and annular LiDAR models. It does not modify or overwrite
the existing simulation code or result files in the parent folder.

## Run

Add this folder to the MATLAB path and run:

```matlab
addpath('D:\lzx\MatlabCode\codexagent_MC\sensitivity_analysis')
outputs = run_sensitivity_analysis;
```

Useful development modes:

```matlab
run_sensitivity_analysis('Mode','smoke','UseParallel',false,'Force',true)
run_sensitivity_analysis('Mode','screen')
run_sensitivity_analysis('Mode','plotonly')
```

The full run uses 3000 scenarios per screening condition and 10000 scenarios
per refined condition. Checkpoints are stored in `results`; rerunning without
`Force=true` resumes compatible completed cases.

## Outputs

- `results/sensitivity_screen_summary.csv`
- `results/sensitivity_screen_comparisons.csv`
- `results/sensitivity_screen_effects.csv`
- `results/sensitivity_refined_summary.csv`
- `results/baseline_validation.csv`
- `results/baseline_validation_refined_N10000.csv`
- `results/effect_direction_validation.csv`
- `results/selected_factors.csv`
- `figures/sensitivity_analysis_main.png`
- `figures/sensitivity_analysis_main.fig`

Run the class-based tests with:

```matlab
results = runtests(fullfile('sensitivity_analysis','tests'));
```
