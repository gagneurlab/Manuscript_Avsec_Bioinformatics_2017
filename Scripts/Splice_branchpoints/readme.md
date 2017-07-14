# Obtain training data for Signal et al

See the following repository: 
- https://i12g-gagneurweb.informatik.tu-muenchen.de/gitlab/avsec/splice_branchpoints

# Run the following scripts for further pre-processing

1. 0_append_id.R
2. 1_preproc_data.R
3. 2_simple_model.R

# Start hyper-parameter tuning with hyperopt

Start workers and run train.py

# Post-process:

Run:
- python postprocess.py
- python interpret_model.py
- Rscript 3_bootstrap_predictions.R

# Plot

Using `Scripts/Figures/Fig3.R` and `Scripts/Figures/Fig4.R`.
