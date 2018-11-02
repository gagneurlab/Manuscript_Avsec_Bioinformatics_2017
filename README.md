# Manuscript_Avsec_Bioinformatics_2017

Code for [Avsec et al, Bioinformatics 2017](https://doi.org/10.1093/bioinformatics/btx727) (Bioarxiv preprint [Avsec et al, Bioarxiv 2017](https://doi.org/10.1101/165183)).

## Installing requirements

### Python

Note, use python>3.5 from anaconda and use virtual environments to not interfere with your default python environment.

```{python}
pip install -r python3_requirements.txt
```

### R

Start R using: `R --vanilla`.

Install CRAN R packages:

```{r}
install.packages(readLines("r_packages.txt"))
```

Install Bioconductor packages:

```{r}
source("https://bioconductor.org/biocLite.R")
sapply(readLines("r_bioc_packages.txt"), biocLite)
```

## General repository organization notes

In each experiment folder training predictive models, the main files are:

- readme.md - contains further instructions
- data.py - contains a `data()` function returning a tuple of train and test-set arrays
- model.py - contains a `model()` function returning a Keras model
- train.py - runs model training and hyper-parameter optimization

R should be started from the repository root.

All data are located either in `Data` (smaller clip data) or in `data/` (everything else).

## Download the data

To download the rest of the data not contained in the repository, run:

```bash
wget https://i12g-gagneurweb.in.tum.de/public/paper/Avsec_Bioinformatics_2017/data.tar.gz
tar xvfz data.tar.gz
```

`data/` includes intermediary results, trained models, as well as model training/test datasets.

## Support

Let me know if you have any problems by creating an issue or sending me an email to avsec-at-in.tum.de.
