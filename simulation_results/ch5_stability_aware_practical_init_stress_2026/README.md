# Stability-aware GN with Algorithm-1 initialization

This folder contains two experiments that retain the complete data-driven
initialization in Algorithm 1.

The two user-pilot columns have equal energy and an exactly controlled
correlation coefficient `rho_p < 1`.  Hence every tested pilot matrix remains
full rank, but its Gram-matrix condition number grows as
`(1 + rho_p)/(1 - rho_p)`.  For every Monte Carlo trial, the code performs
four-phase field recovery, unstructured channel least squares, coarse angular
search, gain projection, and then either conventional LM or stability-aware
GN.  The two solvers receive the same data-driven initial point and the same
noisy observation.

Run `generate_stability_aware_practical_init_stress.m` in MATLAB R2023b or
later.  EPS, PNG, FIG, CSV, and MAT files are written to `results`.

The first experiment applies the Jacobian-conditioned damping floor directly.
It shows that adaptive LM already handles most moderately conditioned cases;
the extra floor mainly reduces rejected updates near the ill-conditioned
boundary and should be viewed as a safeguard rather than a universal speedup.

The second experiment uses the stability analysis as an online training
decision.  It starts with 16 intensity snapshots and checks both the known
pilot Gram-matrix condition and the Jacobian condition at the Algorithm-1
initial point.  If either exceeds its threshold, four complementary pilot
cycles are acquired.  Their cross-correlation cancels that of the initial
block, making the combined 32-snapshot pilot matrix orthogonal.  This adaptive
experiment is the clearer demonstration of how the stability theory can guide
the receiver algorithm.
