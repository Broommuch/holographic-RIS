# Reference-design NMSE-versus-SNR comparison

This folder contains the replacement simulation for the Chapter V comparison
of training-reference designs.

The experiment retains the original two-user channel, 8-by-8 RIS, 32 QPSK
pilots, reference amplitude, reference codebook, Gauss--Newton estimator, and
intensity-domain SNR definition. It compares no reference, constant phase,
random phase, uniformly spaced phase, and codebook-optimized references over
SNRs from 0 to 24 dB.

Each SNR and reference design uses 300 Monte Carlo trials. Common
standard-normal noise and matched randomized solver restarts are used across
the curves within each trial to reduce comparison variance without changing
the marginal distribution at any operating point. The script also calculates
95% Monte Carlo confidence intervals and saves them in the CSV and MAT files,
although the publication figure shows only the mean curves for clarity.

The optimized codebook search selects the uniformly spaced reference in this
configuration. These identical curves are therefore combined in the plot
legend rather than drawn twice.

Running generate_ch5_reference_design_snr_results.m writes EPS, PNG, FIG, CSV,
and MAT outputs to the results folder. Random seeds are fixed for
reproducibility.
