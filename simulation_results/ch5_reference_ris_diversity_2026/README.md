# Chapter V-E: Reference and RIS Measurement Diversity

This folder is the self-contained simulation package used by the fifth
subsection of Chapter V in `main.tex`.

## Provenance

The archived simulations do not contain one recent driver covering every sweep
requested by this subsection. The present driver therefore consolidates the
verified intensity-domain ML detector and channel estimator from:

- `../ch5_perfect_csi_detection_2026/generate_ch5_perfect_csi_results.m`;
- `../ch5_end_to_end_estimated_csi_2026/generate_ch5_end_to_end_results.m`.

Those drivers were themselves based on the newest complete archived programs
`tx_and_rx_test4_1.m` and `tx_and_rx_test5_2.m`. The observation and noise model
remain consistent with the simulation setup in `main.tex`:

`z = |G s + b|^2 + w`, with real post-detection Gaussian noise.

## Experiments

Run `generate_ch5_reference_ris_diversity_results.m` in MATLAB. It produces:

1. normalized energy-domain minimum distance and ML SER versus the number of
   RIS element-level observations;
2. the same metrics versus the number of phase-shifted reference states;
3. array-manifold conditioning and perfect/estimated-CSI SER versus user
   angular separation;
4. an equal-overhead comparison between reference-only diversity and joint
   reference-RIS measurement diversity.

For the equal-overhead comparison, both designs use four states and sixteen
element observations per state. The reference-only design repeats a central
4-by-4 subarray. The joint design uses four state-dependent 16-element masks
over the full 8-by-8 aperture and selects them from a 4096-entry equal-overhead
codebook by maximizing the normalized minimum energy distance.

## Outputs

The `results` folder contains the complete MAT workspace, numerical CSV tables,
and every publication figure in EPS, PNG, and editable MATLAB FIG formats.
