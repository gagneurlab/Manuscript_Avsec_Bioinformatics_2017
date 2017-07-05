# TODO 

- [ ] re-shuffle the provided train, valid and test csv files
  - [ ] Confirm you are getting exactly the same thing as before...
- [ ] re-plot everything with new data


# Running the whole eCLIP pipeline

Final result of the pipeline are the following files:

- .csv files for model training and testing:
  - `"data/eclip/processed/design_matrix/<data split>/<RBP>_extended.csv"`
- Tidy tables of peak locations:
  - `"data/eclip/processed/protein_peak_overlaps.rds"`
  - `"data/eclip/processed/peak_center-gene_mapping.rds"`
- Model test-set predictions: `"data/eclip/processed/predictions/<RBP>/<METHOD>.csv"`


## Pre-process eCLIP data

### Download all the encode data

From the repository root run:

```{bash}
mkdir -p data/eclip/raw
cp Scripts/RBP/Eclip/files.txt data/eclip/raw/
cd data/eclip/raw
xargs -n 1 curl -O -L < files.txt
gunzip *.gz
```

### Pre-process

From the `Scripts/RBP/Eclip` directory run:

```{bash}
snakemake
```

See Snakefile and [snakemake docs](https://snakemake.readthedocs.io/) for more information. You can run things in parallel with `--cores=<ncores>` flag.

## Train the models

Move to the `predictive_models` dir.

```{bash}
cd predictive_models
```

### Start workers

Setup the hyperopt mongo database (see [hyperopt docs](https://github.com/hyperopt/hyperopt/wiki/Parallelizing-Evaluations-During-Search-via-MongoDB))

Edit the `MONGO_IP` and `THIS_DIR` in `hyperopt_worker.bash` and start many instances of it (the following example is in [SLURM](https://slurm.schedmd.com/)):

```{bash}
cd Eclip/predictive_models
for i in {1..20}; do sbatch <path>/hyperopt_worker.bash; done
```

### Start hyper-param optimization main script

Only for 6 RBP's and 3 different methods (maybe a few hours with 20 workers in parallel):
```{bash}
./train_subset_start-end.py
```

For all 112 RBP's (this might take a few days running 20 workers in parallel):

```{bash}
./train_all.py
```

After everything is done, you can check how many trials were actually done for a particular experiment and RBP (sanity check if anything crashed):

```{python}
from concise.hyopt import CMongoTrials
trials = CMongoTrials("RBP__Eclip", "DeepNN_2_UPF1, ip='localhost', kill_timeout=30 * 60)
len(trials)
```

See [concise documentation](https://i12g-gagneurweb.in.tum.de/project/concise/hyopt/) for more info on CMongoTrials object.

### Get test-set predictions

Retrieve the test-set predictions for the best model from each experiment:

```{bash}
./get_predictions.py
```

Note: you might have to edit the `HOST` variable.

## Plot

Now, you can run `Scripts/Figures/Fig2.R` to reproduce Figure 2.
