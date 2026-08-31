# Chapter V-D: End-to-End Detection with Estimated CSI

This folder contains the complete simulation source and generated results used
by the subsection **End-to-End Symbol Detection with Estimated CSI** in
`main.tex`.

## Provenance

The archive contains successful early joint estimation/detection programs
(`test6_4_1.m`, `test7.m`, and the modular `test7_1.m`), but their noise and
reference models predate the manuscript's current formulation. The driver in
this folder therefore combines the two current, publication-oriented baselines:

- `../ch5_channel_estimation_2026/generate_ch5_channel_estimation_results.m`;
- `../ch5_perfect_csi_detection_2026/generate_ch5_perfect_csi_results.m`.

It reuses their damped intensity-domain GN channel estimator and ML/projected
GN/reference-assisted GS data detectors. Training and data measurements use the
same real Gaussian post-detection noise model and intensity-domain SNR definition
as `main.tex`.

## Reproduction

Run `generate_ch5_end_to_end_results.m` in MATLAB. The script is self-contained
and writes all outputs to `results/`:

- perfect- versus estimated-CSI SER for ML, projected GN, and
  reference-assisted GS;
- channel NMSE and end-to-end ML SER versus pilot length;
- ML SER versus an exactly controlled channel NMSE;
- end-to-end performance for five training-reference designs;
- end-to-end performance for five data-reference designs;
- CSV tables and one MAT archive containing all numerical results.

EPS files are referenced by LaTeX; PNG files support quick inspection and FIG
files remain editable in MATLAB. Zero-error Monte Carlo points remain zero in
CSV/MAT and are displayed only at the half-error plotting floor on logarithmic
axes.
