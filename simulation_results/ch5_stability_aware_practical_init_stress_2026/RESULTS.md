# Result summary

All reported results use 300 paired Monte Carlo trials, an 8-by-8 RIS, two
users separated by 2 degrees, a post-detection SNR of 8 dB, and the complete
four-phase initialization in Algorithm 1.  No true parameter is used by the
initializer.

## Direct stability-aware damping

When the pilot Gram-matrix condition number is at most 79, conventional LM and
stability-aware GN converge to essentially the same solution.  At condition
199, stability-aware GN changes the success probability only from 0.887 to
0.890 and the channel NMSE from -11.71 to -11.91 dB.  Thus, once Algorithm 1
provides a useful initial point, adaptive LM already absorbs most of the local
conditioning problem.  The Jacobian floor is a conservative safeguard, not a
general acceleration mechanism.

## Stability-guided adaptive training

The adaptive receiver first collects 16 intensity snapshots.  It acquires 16
complementary snapshots only if either the initial pilot Gram condition exceeds
100 or the Algorithm-1 Jacobian normal-matrix condition exceeds 1000.  The
additional pilot block cancels the initial cross-correlation, so the combined
pilot matrix is orthogonal.

- With pilot Gram condition at most 19, the trigger never fires; the adaptive
  receiver uses 16 snapshots and exactly matches the fixed-16 baseline.
- At condition 39, the trigger probability is 0.103.  The average training
  length is 17.65 snapshots and the channel NMSE improves from -19.45 to
  -20.30 dB.
- At condition 79, the trigger probability is 0.570.  With an average of 25.12
  snapshots, success improves from 0.830 to 0.887 and NMSE improves from
  -16.42 to -19.18 dB.
- At conditions 199 and 399, the pilot-condition test always triggers.  The
  adaptive receiver then matches the fixed-32 benchmark: success is 0.990 and
  NMSE is -28.36 dB.  The fixed-16 success probabilities are only 0.720 and
  0.560, with NMSEs of -9.25 and -4.54 dB, respectively.

These results indicate that the strongest algorithmic value of the stability
analysis is diagnostic and adaptive: it identifies when the initial training
block is insufficient and activates extra diversity only when necessary.  It
does not need to alter a well-conditioned GN trajectory.
