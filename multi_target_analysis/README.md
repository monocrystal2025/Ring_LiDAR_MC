# Multi-UAV annular-beam simulation

This standalone MATLAB module extends the repository's single-target LiDAR
model to paired scene-level multi-UAV simulations. It does not modify the
existing `MC_*` implementation.

## Entry point

```matlab
addpath('multi_target_analysis')

% Fast functional run
outputs = run_multiuav_analysis('Mode','smoke','UseParallel',false,'Force',true);

% Paper-parameter preview (300 scenes by default)
outputs = run_multiuav_analysis('Mode','preview','Force',true);

% Larger grids
run_multiuav_analysis('Mode','screen');
run_multiuav_analysis('Mode','formal','N',10000, ...
    'UseParallel',true,'FormationSpacing',15);
run_multiuav_analysis('Mode','formationopt','N',300, ...
    'TargetCounts',[4 8],'UseParallel',true);
run_multiuav_analysis('Mode','plotonly');
```

The output contains per-target detections and times, scene-level summaries,
equal-width annular-versus-line comparisons, CSV tables, MAT checkpoints,
and a preview figure.

Use `multiuav_validate_single_target(300)` to compare the M=1 limit with
the current single-target implementation. Preview numerical conclusions
are recorded in `PREVIEW_REPORT.md`.

`formationopt` independently searches line and ring beams over the
manuscript grid `wD=5:5:500 mrad`, reports discrete maxima, and performs a
paired five-fold out-of-sample comparison. Its latest interpretation is in
`FORMATIONOPT_REPORT.md`. Use `WDGrid_mrad` for a confirmatory sub-grid;
the default remains the complete manuscript range.

Each target count is simulated independently so that range-cell conflicts
from targets outside the requested M cannot contaminate smaller-M metrics.
The default range resolution is `c*50 ns/2 = 7.4948 m`. Checkpoints are
saved after every beam/path/scenario/M case; non-default formation spacing
is included in the checkpoint name.
