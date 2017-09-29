## Spline transformation paper

**Disclaimer:** I'm in the process of porting the code from my research repository and haven't tested the whole pipeline end-to-end. I will also provide the intermediary raw files soon.

[Avsec *et al* 2017 Bioarxiv](https://doi.org/10.1101/165183)

## Installing requirements

### Python

Note, use python>3.5 from anaconda and use virtual environments to not interfere with your default python environment.

```{python}
pip install -r python3_requirements.txt
```

### R

Normal R packages

```{r}
install.packages(readLines("r_packages.txt"))
```

Bioconductor packages

```{r}
source("https://bioconductor.org/biocLite.R")
sapply(readLines("r_bioc_packages.txt"), biocLite)
```


## General repository organization notes

In each experiment folder, the main files are:

- readme.md - contains further instructions
- data.py - contains a `data()` function returning a tuple of train and test-set arrays
- model.py - contains a `model()` function returning a Keras model
- train.py - runs model training and hyper-parameter optimization

R should be started from the repository root.

## Code for producing the plots

Is located at `Scripts/Figures/`.
