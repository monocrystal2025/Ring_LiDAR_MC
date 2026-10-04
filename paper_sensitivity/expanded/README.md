# Expanded symmetric parameter robustness

This module extends the previous three-row analysis. It keeps the same target bank (moving N=512, static N=2048), 15 s / 60% angular requirement, and 125 mrad reference width. Identical old trials are reused only after source checks. All additional values are physically simulated; no curve smoothing, pulse scaling or replacement of unattained outcomes is used.

Candidates: energy, reflectivity and scan angular speed +/-50%; extinction +/-90%, in 10% increments. The selected range for each parameter is continuous on this discrete grid and symmetric. All interior levels remain. Candidate reference results, excluded points and rejection reasons are retained. Selection uses positive finite estimates and adequate bootstrap reachability, not statistical significance. This is exploratory operating-domain identification, not preregistered validation.

The optional scan-speed column must support at least +/-30% before it enters the expanded main figure. A narrower qualified interval is retained in analysis.qualifiedHalfRange, and all scan-speed data remain in the candidate table. Static blind-zone reversal at -20% already prevents a broad symmetric scan-speed column. Target speed is not substituted only to fill a column, because it does not affect the static experiment.

Scan angular speed regenerates every pulse direction, boundary circles, initial-phase support, scan period and accumulation window at fixed pulse repetition frequency and pulse energy. Baseline regenerated trajectories at 5, 60 and 125 mrad match the old NPY files exactly. Static blind fractions use one complete cycle at each speed, with different elapsed durations.

## MATLAB reproduction

```matlab
addpath('D:/lzx/MatlabCode/codexagent_MC')
addpath('D:/lzx/MatlabCode/codexagent_MC/paper_sensitivity')
addpath('D:/lzx/MatlabCode/codexagent_MC/paper_sensitivity/expanded')
run_expanded_campaign
finish_expanded_robustness
```

Do not start another writer while the campaign is already running. The first command saves each completed task, screens all candidates at 125 mrad, then completes 5:5:150 mrad searches for candidates with positive reference benefits. The second extends candidates whose initial maximum detection fraction is below 65% for either beam to 400 mrad, then analyzes, plots and verifies. Changes in numerical sources invalidate cache reuse. The original trajectory disk and readNPY are required for simulation.

After simulation, plotting can be repeated with plot_expanded_robustness only. Main results and all candidates are in results/all_candidates.csv; range selection is in results/range_decisions.csv; common-success pulse controls are in results/pulse_controls.csv. analysis.mat contains paired bootstrap estimates and the data provenance. verification.json is written only after numerical and artifact checks pass. Figures are exported in PNG, vector PDF and editable FIG formats.

The main figure gives pointwise 95% paired-bootstrap intervals (1000 draws). Simultaneous intervals are also exported. These do not include selection, grid or model uncertainty. No inference is made about every continuous parameter value or arbitrary joint perturbations. Successful-window pulse counts are conditional on each beam's successful reports; the same-success target control is retained.
Angular searches also extend to 400 mrad when either beam has not reached a 65% empirical fraction in the initial grid. This five-percentage-point margin checks near-threshold censoring before the 60% criterion and bootstrap reachability are evaluated; the reported detection requirement remains 60%.

Final verified ranges are E +/-40%, reflectivity +/-40%, extinction +/-90%. The optional scan-speed parameter only supports +/-10% under the stated selection rule and is omitted from the main figure. See manuscript_zh.md for Chinese manuscript paragraphs, physical ranges, boundary evidence and interpretation limits, and completion_audit.md for the final checks. The English manuscript.md is regenerated from analysis.mat; the Chinese discussion is a separately reviewed document. These results use N=512 moving and N=2048 static targets, not the earlier N=3000 batch.
