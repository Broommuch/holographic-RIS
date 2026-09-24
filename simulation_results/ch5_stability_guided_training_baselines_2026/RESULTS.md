# Result summary

The two-panel result compares channel-estimation accuracy and training overhead
for five policies that all use the practical Algorithm 1 initialization and
the same LM refinement.  Each point is based on 1000 paired Monte Carlo trials.
The 95% confidence half-width of every plotted channel-NMSE point is below
0.16 dB.

## Main observations

The reduced-residual rule nearly always stops at the shortest block: its
average length is only 8.3--8.4 measurements.  This behavior is expected
because a residual close to the noise floor does not certify that the physical
parameters are stably estimated.  At 6 dB, for example, it uses 8.42
measurements and obtains -20.37 dB NMSE.

The stability-guided policy responds to this ambiguity by increasing the
training length only where needed:

- at 4 dB, 39.11 measurements give -25.13 dB NMSE;
- at 6 dB, 25.58 measurements give -25.28 dB NMSE;
- at 8 dB, 16.50 measurements give -25.54 dB NMSE;
- at 10 dB, 12 measurements give -26.15 dB NMSE;
- from 12 dB onward, the shortest 8-measurement block is already sufficient.

At 0 and 2 dB, the selected length reaches the 40-measurement budget, but the
-25 dB target remains unattainable at that maximum length.  Fixed-long training
continues to improve with SNR and therefore provides the highest accuracy, but
always pays five times the minimum training overhead.  The oracle curve shows
that the proposed selector follows the correct overhead transition without
using the true channel; the small remaining gap is the price of making the
decision from estimated FIM/CRLB quantities.

## Interpretation

The result isolates the role of the theory more clearly than a damping-only
comparison.  Algorithm 1 supplies a realizable basin-reaching initialization,
LM performs the local refinement, and the uniqueness/stability analysis decides
whether the current measurements contain enough information.  The experiment
therefore demonstrates a measurable accuracy--overhead benefit rather than
attributing a fundamental estimation gain to additional numerical damping.
