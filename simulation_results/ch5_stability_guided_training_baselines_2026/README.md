# Stability-guided adaptive training with practical baselines

This folder contains a self-contained simulation for comparing adaptive
pilot-length policies under the practical initialization in Algorithm 1.  No
truth-perturbed or oracle initialization is used by any estimator.

## Compared policies

All policies use the same nested QPSK pilots, four-phase field recovery,
unstructured channel least squares, angular-grid projection, and parametric LM
refinement.  They differ only in the number of acquired pilot measurements.

- **Fixed short:** always uses 8 intensity measurements.
- **Residual-adaptive:** uses the noise-aware reduced-residual discrepancy
  rule with threshold 1.10.
- **Stability-guided:** selects the first length satisfying the predicted
  channel CRLB-NMSE target of -25 dB, the 1-degree DOA-CRLB target, and the
  normal-matrix condition-number limit of 5000.
- **Oracle stopping:** uses the true channel NMSE only to provide an
  unattainable lower-overhead reference.
- **Fixed long:** always uses 40 intensity measurements.

The experiment uses an 8-by-8 RIS, two users separated by 2 degrees, candidate
lengths from 8 to 40 measurements in increments of four, SNRs from 0 to 20 dB,
and 1000 paired Monte Carlo trials per SNR.  The same noise realization is used
by all five policies in each trial.

## Reproduction

Run `generate_stability_guided_training_baselines.m` in MATLAB R2023b or
later.  Set the environment variable `HOLO_QUICK_CHECK=1` for a 20-trial code
check.  The full run uses a four-worker process pool when the Parallel
Computing Toolbox is available and remains reproducible through per-trial
random seeds.

The `results` folder contains the vector EPS figure, a 300-dpi PNG preview, the
MATLAB figure, complete numerical tables, selection distributions, and the MAT
workspace.
