# Extended reference and RIS measurement-diversity simulations

This folder contains the revised simulations for the final subsection of
Chapter V. The driver is based on
`../ch5_reference_ris_diversity_2026/generate_ch5_reference_ris_diversity_results.m`
and retains its two-user QPSK channel and intensity-domain SNR definition.

The RIS-size sweep now uses 13 square arrays ranging from 2-by-2 to 32-by-32
(4 to 1024 element-level observations). Each point uses 5000 Monte Carlo
symbol-vector trials, or 10000 user-symbol decisions, and is evaluated in
chunks to control memory usage.

The reference-state sweep now covers 1 through 16 uniformly spaced phase
states. Each point uses 10000 symbol-vector trials. All two-metric figures,
including the unchanged angular-separation experiment, use a horizontal
two-panel layout suitable for double-column placement.

For the equal-overhead reference-only and joint reference--RIS comparison,
the distance-only bar panel is omitted. The retained SER figure reports the
proposed reference-assisted GS detector for both measurement designs, with
exact ML shown as a benchmark on the same paired observations. Each SNR point
uses 10000 symbol-vector trials.

Running `generate_ch5_reference_ris_diversity_extended_results.m` writes EPS,
PNG, FIG, CSV, and MAT outputs to the `results` folder. Random seeds are fixed
for reproducibility.
