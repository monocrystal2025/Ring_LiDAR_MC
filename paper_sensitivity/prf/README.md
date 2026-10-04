# PRF sensitivity, 1–10 kHz

User-confirmed constraint: fixed single-pulse energy 150 microjoules. Average transmitted power therefore varies from 0.15 to 1.50 W. Scan angular speed stays at 2*pi rad/s. Other physics, target positions/velocities and normalized initial phases match the expanded experiment (512 moving, 2048 static targets). This is not the older N=3000 batch.

Every PRF regenerates the spiral and boundary circles on its actual pulse time grid. Dynamics use dt=1/f, the 15-second limit is converted to pulses using the current f, and accumulation length is ceil(f*width/omega). Static coverage uses a complete cycle. At 5 kHz the already validated raw data are reused after numerical source and bank signature checks. No pulse counts or outcomes are algebraically scaled from another PRF.

Run `run_prf_campaign` after adding this directory to the MATLAB path. It completes all 10 reference cases at 125 mrad, then all 5:5:150 mrad angular grids without filtering on benefit sign. If either beam's empirical maximum is below 65%, the search extends to 400 mrad, retaining a margin above the unchanged 60% requirement. Searches use raw first crossings and adjacent linear interpolation, not monotonic smoothing. Unattained widths and mathematically undefined ratios remain NaN with status labels; they are not imputed.

`analyze_prf_sensitivity` reports paired 95% bootstrap intervals (1000 resamples), absolute blind-fraction differences, relative gains, conditional pulse counts and common-success controls. Intervals are conditional on finite values when >=97.5% of resamples are valid. These are pointwise Monte Carlo sampling intervals, not simultaneous, selection-adjusted or physical-model error bounds. Changing PRF resamples the scan; exact target-by-target monotonicity is not imposed because the sampled pulse grids are not generally nested.

`plot_prf_sensitivity` exports PNG, vector PDF and editable FIG. A positive blind-fraction difference means fewer blind targets with the annular beam. This difference remains meaningful if the line-beam blind fraction is zero, whereas the relative reduction is undefined. Pulse counts are conditional on successful detection windows and do not directly measure probability or detection speed.

All data are independent of the prior expanded figure and do not overwrite its results. Use only one writer at a time for this directory.
