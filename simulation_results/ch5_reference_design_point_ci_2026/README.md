# Reference-design point-estimate comparison

This folder contains the replacement simulation for the Chapter V comparison
of training-phase reference designs.

The experiment retains the signal model and baseline parameters of
generate_ch5_channel_estimation_results.m. The statistical and presentation
changes are:

- 300 Monte Carlo trials per reference design instead of 30;
- common random numbers and matched randomized GN restarts across designs;
- storage of every trial-level NMSE value;
- 95% Monte Carlo confidence intervals for the mean channel NMSE;
- a two-panel point-estimate plot instead of bars and connected categorical
  points.

Panel (a) shows the mean channel NMSE with its 95% confidence interval. Panel
(b) shows the minimum singular value of the training Jacobian. The singular
value below 1e-12 in the no-reference case is plotted as numerical zero, which
keeps the reference-assisted values readable on a linear scale. Since the
reference designs are unordered categories, no line is drawn between them.

Running generate_ch5_reference_design_point_ci_results.m writes the EPS, PNG,
FIG, CSV, and MAT outputs to the results folder. Random seeds are fixed for
reproducibility.
