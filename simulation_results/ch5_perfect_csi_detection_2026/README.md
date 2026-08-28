# Chapter V-C: Perfect-CSI Symbol Detection

This folder contains the complete simulation source and generated results used
by the subsection **Symbol Detection Performance with Perfect CSI** in
`main.tex`.

## Provenance

The starting point was the newest complete multiuser symbol-detection program
in the archive:

`../01_核心_全息收发联合仿真/02_v2波形域与符号级_2026/tx_and_rx_test4_1.m`

That program already provides verified ML/GN multiuser QPSK results.  The
publication driver here retains its 8-by-8 RIS channel construction, exhaustive
finite-alphabet ML detector, and reference-aided initialization.  It supplements
the missing reference-design, energy-distance, and modulation-order sweeps, and
uses the real Gaussian post-detection noise model stated in `main.tex` throughout
so that all curves in the subsection share one SNR definition.

## Reproduction

Run `generate_ch5_perfect_csi_results.m` in MATLAB.  The script has no external
function dependencies and writes all outputs to `results/`:

- detector SER versus SNR for ML, projected GN, and reference-assisted GS;
- ML SER and normalized energy-domain minimum distance for five reference
  designs;
- ML SER versus the normalized minimum distance over a 32-entry equal-power
  phase/spatial-support reference codebook;
- QPSK and 16-QAM ML SER under the same optimized reference;
- the underlying CSV tables and one MAT archive.

EPS files are the versions referenced by LaTeX.  PNG files are provided for
quick inspection and FIG files retain editable MATLAB figures.

Zero-error Monte Carlo points remain zero in the CSV/MAT data.  Only their plot
markers are displayed at the conventional half-error floor
`0.5/(number of detected user symbols)` so that they remain visible on a
logarithmic axis.
