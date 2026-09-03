# Perfect-CSI reference-design SER comparison

This folder contains the replacement simulation for the reference-design
comparison in the perfect-CSI subsection of Chapter V.

The experiment retains the original two-user QPSK channel, 8-by-8 RIS, four
reference states, reference amplitude, ML detector, reference codebook, and
intensity-domain SNR definition. It compares no reference, constant phase,
random phase, uniformly spaced phase, and codebook-optimized references over
SNRs from -15 to 3 dB.

Each SNR uses 20000 transmitted symbol-vector trials, corresponding to 40000
user-symbol decisions. The same transmitted symbols and standardized noise are
used across reference designs at each SNR. Detection is evaluated in chunks to
control memory usage.

The publication figure contains only the SER curves. It uses the same color,
line, and marker mapping as the training-reference NMSE figure in the preceding
subsection. The energy-distance panel is removed; the normalized distances are
retained in the CSV, MAT, and manuscript discussion.

Empirical zero-error points are displayed at the conventional half-event level
0.5 divided by the number of symbol decisions, so they remain visible on the
logarithmic axis. The CSV records the error counts, empirical SER, plotting
value, zero-error indicator, and Wilson 95% confidence interval.

Running generate_ch5_perfect_csi_reference_ser_results.m writes EPS, PNG, FIG,
CSV, and MAT outputs to the results folder. Random seeds are fixed for
reproducibility.
