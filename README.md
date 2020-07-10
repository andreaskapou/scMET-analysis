# scMET-analysis
This repository contains scripts to reproduce the results of the manuscript:
__scMET: Bayesian modelling of DNA methylation heterogeneity at single-cell resolution__.

The github repository for the scMET package is [https://github.com/andreaskapou/scMET](https://github.com/andreaskapou/scMET)

## Installation
```R
# install.packages("devtools")
devtools::install_github("andreaskapou/scMET")
```

# Structure
The structure of the project is the following:

* `/ecker2017/` folder: Contains scripts for the analysis of the mouse frontal cortex dataset by [Luo et al. (2017)](https://science.sciencemag.org/content/357/6351/600.abstract).
* `/gastrulation/` folder: Contains  scripts for the analysis of the scNMT-seq gastrulation dataset by [Argelaguet et al. (2019)](https://www.nature.com/articles/s41586-019-1825-8).
* `/synthetic/` folder: Contains scripts for analysis on synthetic data.
