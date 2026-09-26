# Result summary

The production run uses 300 paired trials per angular separation for the GN
experiment and 100 channel trials with 100 data vectors per channel for the
end-to-end experiment.

## What changed

The previous script initialized every GN run by perturbing the true DOAs and
path gains.  It also assigned identical pilot sequences to the two users.  The
new script never accesses the true parameters during initialization.  It uses
four phase-shifted intensity measurements to recover each complex field
snapshot, estimates the two unstructured user channels by least squares, and
then performs the angular-grid and gain-projection steps in Algorithm 1.

Because identical user pilots give a rank-one pilot matrix and cannot support
that per-user LS initialization, the new QPSK pilot matrix is correlated but
full rank.  It has rank 2 and a Gram-matrix condition number of 8.55.  Both the
old and new experiments retain 16 total intensity snapshots.

## Main observations

- At 12 dB training SNR, the data-driven initializer has approximately
  1.01--1.05 degree joint DOA RMSE across the tested angular separations.  GN
  reduces the final RMSE to approximately 0.51--0.56 degree.
- The three GN variants now have the same final channel NMSE and recovery
  probability to plotting precision.  Their success probabilities range from
  0.95 to 0.993.  Stability-aware damping needs more iterations at the closest
  separations (about 7.5 versus about 4.9 for undamped GN and conventional LM)
  because the practical initializer already places the iterate in a reliable
  basin.
- Relative to the old oracle/rank-one-pilot result, the stability-aware channel
  NMSE improves by roughly 6.7--8.1 dB for separations up to 8 degrees.  This is
  primarily a pilot-identifiability and initialization effect, not a clean
  damping-only comparison.
- In the end-to-end experiment, conventional LM and stability-aware GN yield
  essentially identical channel NMSE and SER.  The SER is 9.5e-4 at 0 dB,
  1.0e-4 at 4 dB, and zero observed errors from 8 dB onward.  At 0 dB the final
  channel NMSE is -13.56 dB, compared with -9.79 dB immediately after the
  Algorithm-1 initialization.

## Interpretation

The old figure is not a valid demonstration of the complete receiver because
its rank-one pilot matrix prevents Algorithm 1's stated per-user channel
initialization, while its truth-perturbed start bypasses that obstacle.  The
new result is internally consistent and shows that the proposed initialization
is effective.  It also shows that, under this setting, stability-aware damping
does not improve the converged solution once the initializer is already inside
the attraction basin.  A separate damping stress test should therefore vary a
full-rank pilot correlation, training SNR, or initialization quality while
always starting from Algorithm 1; otherwise pilot identifiability and optimizer
robustness remain confounded.
