## Spline transformation paper

Disclaimer: I haven't tested the whole pipeline end-to-end. Hence you might encounter some dependency bugs or so.

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
