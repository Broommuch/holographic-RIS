# Theory-guided algorithm simulations

This folder contains the new Chapter V experiments that connect the analysis
in Section III to the receiver algorithms in Section IV.

Run generate_ch5_theory_guided_algorithm_results.m in MATLAB R2023b or later.
The script produces three publication figures and their CSV/MAT data in the
results directory:

- stability-aware parametric Gauss--Newton convergence and recovery success;
- distance-certified finite-alphabet GS accuracy and iteration count;
- end-to-end detection accuracy and nonlinear workload.

The random seeds and all simulation parameters are fixed in the script.
