# Supplementary material

**When can a dataset tell Lindley-type count models apart? Sample-size limits and the behaviour of AIC**

Ahmad A. Hanandeh. Code and data accompanying the manuscript submitted to *Statistical Papers*.

All code is base R (it uses only the `parallel` package that ships with R), plus Python 3 scripts: `matched_kl_series.py` uses only the standard library, and the scripts in `python_theory/` (symbolic and numerical checks of Theorems 1–5) need `sympy`, `mpmath`, `numpy`, `scipy`, `matplotlib` and `rdata`.

## R package

`countfrontier_0.1.0.tar.gz` — install with `R CMD INSTALL countfrontier_0.1.0.tar.gz`. Main functions:

| Function | Purpose |
| --- | --- |
| `cf_families()` | available families (Poisson, Poisson–Gold, Poisson–Lindley, –Shanker, –Akash, –Sujatha, –Amarendra, –Aradhana, –Ishita, –Suja, negative binomial, zero-inflated Poisson) |
| `cf_fit(freq)` | maximum likelihood fits and AIC |
| `cf_frontier(true, mean, rival)` | frontier n* and universal floor |
| `cf_lookup()` | table of frontiers for many pairs and means |
| `cf_npmle(freq)` | NPMLE of the Poisson mixing law with optimality certificate |
| `cf_gap_test(freq, family)` | test of a family against the NPMLE |
| `cf_mcs(freq)` | model confidence set |
| `cf_audit(freq)` | all of the above in one call |

Example: `print(countfrontier::cf_audit(c(187, 77, 40, 17, 6, 2, 1)))`.

## Scripts reproducing the paper

| File | Purpose |
| --- | --- |
| `frontier_core.R` | core functions (shared by all scripts and by the package) |
| `datasets.R` | the ten benchmark frequency tables with sources |
| `01_frontier.R` | Table 2 (frontier and floor), Table 1 (crossings), Figure 2 curves |
| `02_simulation.R` | Vuong power (Table 4, Figure 3) and model confidence set (Table 5) |
| `02b_simulation_large.R` | AIC selection with seven candidates (Figure 4) and NPMLE gap test (Table 6) |
| `03_audit.R`, `03b_pgd_frontier.R` | audit of the ten benchmark datasets (Tables 8–9) |
| `review.R` | review of 42 published claims (Table 7, Figure 5) |
| `general_theory_check.R` | numerical verification of Theorem 3 and Proposition 3 for 21 pairs |
| `06_lookup.R` | lookup table (Appendix B, Table B1) and the coefficients kappa |
| `07_lifetime.R` | lifetime laws: frontier curves and unit dependence (Table 10, Figure 6) |
| `04_figures.R`, `04b_figures_new.R` | Figures 1–6 |
| `05_tables.R`, `05b_tables_new.R` | LaTeX tables from the saved results |
| `matched_kl_series.py` | exact rational series for the Poisson–Gold family |
| `python_theory/` | Data S3: symbolic and numerical checks of Theorems 1–5 (see its README) |

## Order of running

```
Rscript 01_frontier.R
Rscript 02_simulation.R S1; Rscript 02_simulation.R S3
Rscript 02b_simulation_large.R all
Rscript 03_audit.R; Rscript 03b_pgd_frontier.R
Rscript review.R; Rscript general_theory_check.R
Rscript 06_lookup.R; Rscript 07_lifetime.R
Rscript 04_figures.R; Rscript 04b_figures_new.R
Rscript 05_tables.R; Rscript 05b_tables_new.R
python3 matched_kl_series.py
```

Results go to `results/` (`lookup_table.csv` holds the full lookup table), figures to `figures/` and tables to `tables/`. On a two-core machine, `02b_simulation_large.R` takes several hours (mostly the gap test) and `03_audit.R` about one hour; the other scripts take minutes. Seeds are fixed in every script.
