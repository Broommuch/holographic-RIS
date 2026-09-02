# Chapter V-B: High-Monte-Carlo DOA RMSE and CRLB

This folder contains the replacement simulation for the DOA RMSE/CRLB figure
in Chapter V-B.

The driver is derived from
`../ch5_channel_estimation_2026/generate_ch5_channel_estimation_results.m`. It
retains the same channel, reference, post-detection Gaussian-noise model,
Gauss--Newton channel estimator, DOA manifold fitting, and nuisance-parameter
CRLB. The only statistical change is that the SNR sweep uses 1000 independent
Monte Carlo realizations instead of 30.

The CSV file also records a delta-method 95% confidence interval for the RMSE
and the RMS magnitude of the empirical DOA bias. The plotted error bars are
these confidence intervals.

Run `generate_ch5_doa_crlb_high_mc_results.m` in MATLAB. The `results` folder
contains the numerical CSV/MAT data and the figure in EPS, PNG, and editable
MATLAB FIG formats.
