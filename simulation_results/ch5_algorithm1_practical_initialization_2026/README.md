# Algorithm-1 practical-initialization simulations

This folder re-evaluates the Chapter V channel-estimation and end-to-end
experiments without truth-perturbed (oracle) initialization.

Run `generate_algorithm1_practical_initialization_results.m` in MATLAB R2023b
or later.  The total number of intensity snapshots remains `T_p = 16`.  They
are organized as four pilot vectors, each observed with the four reference
phases `0`, `pi/2`, `pi`, and `3pi/2`.  The two user-pilot columns are
correlated but full rank.  Consequently, the script can execute the actual
Algorithm-1 initialization:

1. recover the complex field by opposite-phase energy differences;
2. estimate the unstructured user channels by least squares;
3. obtain coarse DOAs from a two-dimensional angular grid;
4. estimate path gains by steering-vector projection; and
5. run the three Gauss--Newton variants from the same data-driven initial point.

The script writes EPS, PNG, FIG, CSV, and MAT outputs to `results`.  It also
loads the published oracle-initialized CSV files from the preceding simulation
folder and writes diagnostic old-versus-new comparisons.  These comparisons
must be interpreted with the pilot-design distinction in mind: the old
rank-one pilot matrix cannot support Algorithm 1's per-user unstructured
initialization, whereas the new matrix is deliberately full rank.
