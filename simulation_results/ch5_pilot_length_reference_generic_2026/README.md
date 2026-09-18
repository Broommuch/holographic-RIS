# Pilot-count and reference-diversity experiment

Run `generate_ch5_pilot_length_reference_generic_results.m` in MATLAB to
regenerate the figure and data in `results/`.

The experiment uses the per-element two-user unstructured training model
`z = |X_p h + b|^2 + w`, which is the atomic block of the full RIS model.
All pilot lengths are prefixes of one QPSK sequence. The no-reference,
constant, random-phase, and four-phase designs use the same channel and noise
draws in each of 10000 paired Monte Carlo trials. Panel (a) reports the empirical
noiseless recovery rate of a multistart damped Gauss--Newton solver; panel (b)
reports channel NMSE at 12 dB. The lines at `L`, `2L`, and `3L` mark the field,
regular-local, and global-unstructured counting thresholds, respectively.

The empirical recovery rate is an algorithmic statistic for random instances;
it is not itself a proof of uniform global injectivity.
