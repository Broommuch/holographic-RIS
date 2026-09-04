# End-to-end reference-design SER curves

This folder contains the replacement simulation for the reference-design
figure in the estimated-CSI subsection of Chapter V. It is adapted from
`../ch5_end_to_end_estimated_csi_2026/generate_ch5_end_to_end_results.m` and
retains the same two-user channel, 8-by-8 RIS, 32-symbol pilot, four data
reference states, QPSK alphabet, Gauss--Newton channel estimator, and
mismatched-CSI ML detector.

The first panel sweeps the training SNR while fixing the data SNR at -6 dB,
thereby isolating the effect of the training-reference design. The second
panel uses one common estimated-channel ensemble obtained at a training SNR
of 12 dB and sweeps the data SNR, thereby isolating the effect of the
data-reference design.

Each point averages 20 independently estimated channels and 250 transmitted
symbol vectors per channel, corresponding to 10000 user-symbol decisions.
Pilot-noise and data-noise realizations are paired across reference designs
to reduce comparison variance. The color, line, and marker mapping matches
the earlier reference-design figures. Because the uniform and optimized data
references coincide in this configuration, they are represented by one curve
in the second panel; the distinct optimized training reference uses a purple
triangle curve in the first panel.

Running `generate_ch5_end_to_end_reference_ser_snr_results.m` writes EPS, PNG,
FIG, CSV, and MAT outputs to the `results` folder. Random seeds are fixed for
reproducibility.
