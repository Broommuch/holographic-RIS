# Spatial-reference modulation-order comparison

This experiment evaluates perfect-CSI exact-ML detection for normalized QPSK,
16-QAM, 64-QAM, and 256-QAM under the fixed simultaneous spatial four-phase
reference used in Chapter V. It uses the same two-user 8-by-8 surface, channel,
reference amplitude, and intensity-domain SNR definition as the current paper.

The codeword-difference subspace projection is lossless for Euclidean ML
detection and is used only to make the high-order codebooks computationally
manageable. Each SNR point contains 10,000 symbol-vector trials, or 20,000
user-symbol decisions. Empirical zero-error points are omitted from the plot.

Run `generate_spatial_modulation_order_results.m` to generate the EPS, PNG, FIG,
CSV, and MAT files in `results`. Use
`replot_spatial_modulation_order_results.m` to regenerate only the figure from
the saved MAT file.
