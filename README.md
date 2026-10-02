# Hierarchical beta-binomial dose-response for Salmonella

Master’s thesis analysis: how the **probability of infection** rises with dose, and how that curve **differs by Salmonella strain**.

The model is a hierarchical beta-binomial. Strains share a common mean, each strain has its own intercept, and a **new, unseen strain** can be predicted from the same hierarchy. Inference is Bayesian (OpenBUGS, three chains).

This repository is that thesis notebook, turned into a pipeline you can rerun. The model file, priors, seed, and MCMC settings are unchanged.

Thesis figures (PDF): [infection curves and ED50](docs/figures/HBB_Pinf_with_ED50_segments.pdf) · [log10(α)–log10(β) contours](docs/figures/HBB_cont.pdf) · [u histograms](docs/figures/HBB_U_with_CI_in_Title.pdf) · [MCMC traces](docs/figures/HBB_trace.pdf)

## What it answers

Given outbreak counts \(Y\) infected out of \(N\) exposed at a recorded dose, estimate \(P(\text{infection} \mid \text{dose}, \text{strain})\) for:

- four serovars in the fitted subset (*S. enteritidis*, *S. heidelberg*, *S. oranienburg*, *S. typhimurium*)
- a pooled “overall” curve
- a new strain drawn from the same population of strains

Host status in the published fit is **Normal** only (susceptible hosts are excluded, as in the thesis).

## Data

`data/raw/salmonella.csv` — outbreak rows with:

| Column | Meaning |
|---|---|
| `t` | Serovar |
| `S` | Host status (`Normal` / `Susceptible`) |
| `log10dose` | Log10 of the ingested dose |
| `N` | Number exposed |
| `Y` | Number of cases |

## Model

BUGS code: [`inst/bugs/hierarchical_dose_response_model_sigmapriors.txt`](inst/bugs/hierarchical_dose_response_model_sigmapriors.txt)

- Infection probability is the beta-Poisson / beta-binomial form used in the thesis (`p_inf` from `alpha`, `beta`, and dose).
- Strain-level `w[k]`, `z[k]` sit under global `w_0`, `z_0` with half-normal `sigma_w`, `sigma_z`.
- `wnew`, `znew` are the predictive draw for a new strain.

MCMC: seed `123`, `n.chains = 3`, `n.iter = 10000`, `n.burnin = 1000`, `n.thin = 2`.

## How to rerun

1. Install [OpenBUGS](https://www.mrc-bsu.cam.ac.uk/software/bugs/openbugs/) and the R packages:

```r
install.packages(c("R2OpenBUGS", "ggplot2", "MASS", "tidyr", "gridExtra", "coda", "here", "testthat"))
```

2. From the repository root:

```r
testthat::test_dir("tests/testthat")
```

```bash
Rscript scripts/run_pipeline.R
```

The first fit writes `data/derived/bugs_fit.rds`. Later runs reuse it. Pass `--refit` to sample again with the same settings.

Figures from a new run go to `output/`. The PDFs in `docs/figures/` are the thesis plots.

## Layout

| Path | Role |
|---|---|
| `R/` | Data prep, OpenBUGS fit, figures |
| `scripts/run_pipeline.R` | One-command reproduce |
| `inst/bugs/` | Model (source of truth) |
| `data/raw/` | Outbreak table |
| `docs/figures/` | Thesis figures |
| `analysis/thesis-notebook/` | Original R Markdown, unmodified |

## License

MIT. Outbreak figures are from the thesis data file in `data/raw/`.
