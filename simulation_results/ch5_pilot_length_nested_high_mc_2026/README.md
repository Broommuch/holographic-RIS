# Nested high-Monte-Carlo pilot-length simulation

This folder contains the replacement simulation for the Chapter V
channel-NMSE-versus-pilot-length figure.

The original channel-estimation script already used prefixes of one
maximum-length QPSK pilot sequence. However, it independently redesigned the
reference phase matrix at each pilot length, used independent noise and solver
randomness across operating points, and averaged only 30 Monte Carlo trials.
Consequently, both the number and the type of measurements changed along the
horizontal axis.

The file generate_ch5_pilot_length_nested_results.m isolates the pilot-length
effect by:

- using prefixes of one 48-symbol QPSK pilot sequence;
- using prefixes of one jointly selected 48-row reference matrix;
- pairing the standard-normal noise and randomized GN restarts across pilot
  lengths within each trial;
- increasing the number of Monte Carlo trials from 30 to 300;
- reporting the linear-domain Monte Carlo mean and its approximate 95% interval
  after conversion to dB.

The script writes EPS, PNG, FIG, CSV, and MAT outputs to the results folder.
Random seeds are fixed for reproducibility.
