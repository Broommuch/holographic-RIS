# Chapter V-B channel-estimation simulations

This folder contains the paper-oriented simulations used by Chapter V,
Subsection B of `main.tex`.

The entry point is `generate_ch5_channel_estimation_results.m`. It was derived
from the latest verified channel/DOA result script
`tx_and_rx_test5_2.m`, while changing the observation noise to the additive
real Gaussian intensity-domain model assumed in the paper.

Running the script creates `results/` with:

- channel NMSE versus SNR;
- DOA RMSE and the CRLB versus SNR;
- channel NMSE versus pilot length;
- channel NMSE versus reference amplitude;
- channel NMSE and Jacobian conditioning for five reference designs;
- EPS, PNG, FIG, CSV, and MAT versions of the results.

The baseline uses two users, an 8-by-8 RIS, 32 QPSK pilot observations, and 30
Monte Carlo trials per operating point. Random seeds are fixed in the script.
The angular-separation sweep from the original manuscript placeholder is not
included because the latest full-channel estimator separates user channels by
their pilots; such a sweep should be added only after the joint parametric
DOA/gain estimator is used consistently.
