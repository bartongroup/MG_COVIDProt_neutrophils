# Large-scale proteomics of COVID-19 neutrophils

Maintainer: Marek Gierlinski (M.Gierlinski@dundee.ac.uk)

Software to accompany the manuscript "Investigating the impact of dipeptidyl peptidase-1 inhibition in humans using multi-omics: results from a double-blind, randomized, placebo-controlled trial in patients with COVID-19" by Long et al. (2024), in preparation.

## Usage

We suggest using RStudio. Start in the top project directory. The first step is to create environment using 'renv':

```
install.packages("renv")
renv::restore()
```

This will install all necessary packages. Run the `targets` pipeline.

```
targets::tar_make()
```

This will analyse data, create figures (in directory `manuscript_fig`) and a report (in directory `doc`). 
