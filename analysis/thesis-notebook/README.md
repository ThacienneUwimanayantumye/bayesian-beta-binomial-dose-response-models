# Thesis notebook (unchanged)

This folder holds the original R Markdown analysis from the master’s thesis.

The BUGS model, priors, seed (`123`), and MCMC settings (`3` chains, `10 000` iterations, `1 000` burn-in, thin `2`) are what they were.

To knit it, set the R working directory to the **repository root** and copy `data/raw/salmonella.csv` and `inst/bugs/hierarchical_dose_response_model_sigmapriors.txt` beside the notebook, or run the pipeline instead:

```r
Rscript scripts/run_pipeline.R
```
