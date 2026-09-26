# Spatial-reference-consistent Chapter V simulations

This folder regenerates the Chapter V results after replacing repeated
four-phase measurements of one field coordinate with the intended hardware
model: one fixed reference phase per receiving unit, with different units
operating simultaneously.

- `generate_spatial_doa_crlb_results.m`: DOA RMSE and CRLB with a fixed spatial
  four-state phase map.
- `replot_valid_pilot_identifiability.m`: retains the two hardware-consistent
  columns of the existing 10,000-trial unstructured pilot-count experiment.
- `generate_spatial_stability_guided_training.m`: linearized spatial-code
  initialization, parametric LM refinement, adaptive pilot selection, and
  end-to-end detection.
- `generate_spatial_detection_results.m`: perfect-CSI ML, linearized, and GS
  detection with exact-distance stopping.
- `generate_spatial_pairwise_distance_results.m`: direct Monte Carlo validation
  of the pairwise-error formula over fixed spatial phase maps.
- `generate_spatial_diversity_energy_results.m`: receiving-unit diversity,
  angular separation, and reference-to-signal energy allocation.

All publication figures and numerical summaries are stored in `results/`.
