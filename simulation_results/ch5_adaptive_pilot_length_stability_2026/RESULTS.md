# Result summary

The experiment uses an 8-by-8 RIS, two users separated by 2 degrees, four
reference phases per pilot cycle, and 300 paired Monte Carlo trials per SNR.
The shortest and longest training lengths are 8 and 40 intensity measurements.
The adaptive receiver targets a predicted channel CRLB-NMSE of -25 dB, a DOA
CRLB below 1 degree, and a Jacobian normal-matrix condition number below 5000.

## Main result

At 0 and 2 dB, even the maximum candidate length cannot reach the -25 dB
target, so the selector correctly uses all 40 measurements.  At 4 dB it uses
38.91 measurements on average and nearly matches fixed long training.

In the transition region, the selector substantially reduces overhead while
maintaining approximately the requested accuracy:

- at 6 dB, 25.67 measurements give -25.27 dB NMSE, compared with -20.09 dB
  for fixed short training and -27.02 dB for fixed long training;
- at 8 dB, 16.69 measurements give -25.60 dB NMSE, compared with -22.25 dB
  and -29.40 dB for the short and long baselines;
- at 10 dB, 12 measurements give -26.12 dB NMSE, compared with -24.33 dB and
  -31.17 dB for the two baselines.

From 12 dB onward, the shortest 8-measurement block already satisfies the
target, so the adaptive receiver adds no training.  Its NMSE and success
probability then coincide exactly with the fixed-short baseline.

## Agreement with the theory

The predicted channel CRLB-NMSEs at 6, 8, and 10 dB are -25.37, -25.44, and
-26.00 dB.  The corresponding empirical NMSEs are -25.27, -25.60, and
-26.12 dB.  This close agreement is the key result: the stability/FIM analysis
predicts a small sufficient training length without using the true channel or
the empirical estimation error.

The experiment therefore demonstrates a clearer algorithmic role for the
analysis than additional GN damping.  It turns the theoretical conditioning
and CRLB quantities into an online stopping rule that trades training overhead
against a prescribed channel-accuracy target.
