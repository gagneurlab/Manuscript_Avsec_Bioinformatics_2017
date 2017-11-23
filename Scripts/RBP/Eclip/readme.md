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

All the experiments utilizing DNN's were done through hyperopt. To start the master, run

```{bash}
./train_all.py
```

Note that the this might take a few days with 20 workers in parallel on 4 CPU threads each.:

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

This will save the csv file of predictions with columns `y_true,y_pred` to `data/eclip/processed/predictions/{rbp}/{method}.csv`.

Note: you might have to edit the `HOST` variable.

### Score new sequences

The best models were saved to: `data/eclip/models`. To score your own sequences for binding affinities, use: TODO 

## Plot

Now, you can run `Scripts/Figures/Fig2.R` to reproduce Figure 2.

## TODO

- [ ] add the best models to `data/eclip/models`
- [ ] run `dump_best_models.py` - to save all the best models
- [ ] add scoring your own sequences
