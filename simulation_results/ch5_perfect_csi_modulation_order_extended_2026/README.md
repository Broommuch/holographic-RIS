# Extended modulation-order comparison under perfect CSI

This folder contains the replacement simulation for the modulation-order
figure in the perfect-CSI subsection of Chapter V.

The experiment retains the original two-user channel, 8-by-8 RIS, four
reference states, optimized data-phase reference, ML detector, and
intensity-domain SNR definition. The modulation comparison is extended from
QPSK and 16-QAM to QPSK, 16-QAM, 64-QAM, and 256-QAM.

For two users, 256-QAM produces 65536 joint symbol candidates. The script
performs exact ML detection after projecting the intensity codebook onto the
real affine subspace spanned by all codeword differences. Orthogonal noise
components add the same term to every Euclidean ML metric, so this projection
does not change the ML decision. A residual check verifies that every codeword
difference lies in the reduced subspace to numerical precision.

Each SNR uses 5000 transmitted symbol-vector trials, corresponding to 10000
user-symbol decisions. The same transmitted candidate indices and standardized
noise realizations are reused across SNR values for each modulation order.
Empirical zero-error points are displayed at the half-event level of 5e-5.
The CSV records the error counts, empirical and plotted SER values, Wilson 95%
confidence intervals, reduced-subspace ranks, and projection residuals.

Running generate_ch5_modulation_order_extended_results.m writes EPS, PNG, FIG,
CSV, and MAT outputs to the results folder. Random seeds are fixed for
reproducibility.
