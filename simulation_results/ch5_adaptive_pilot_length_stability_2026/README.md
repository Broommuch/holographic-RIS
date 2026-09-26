# Stability-guided adaptive pilot length

This experiment evaluates the stability analysis as a minimum-training-length
selector rather than as an additional damping term.

The receiver starts from the shortest full-rank pilot block and increases the
number of four-phase pilot cycles one at a time.  For every candidate length it
runs the data-driven initialization in Algorithm 1 and constructs the local
FIM.  The parameter covariance bound is propagated through the channel
Jacobian to predict the channel CRLB-NMSE; the gain-nuisance Schur complement
also provides a DOA CRLB.  The first length satisfying the channel-NMSE target,
the DOA target, and the Jacobian-conditioning threshold is selected.
Conventional LM is then run from the corresponding Algorithm-1 initial point.

Run `generate_adaptive_pilot_length_results.m` in MATLAB R2023b or later.
All figures and numerical results are saved under `results`.  The script does
not modify `main.tex`.
