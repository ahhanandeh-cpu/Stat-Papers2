## =====================================================================
##  Audit of ten benchmark datasets
##  For each dataset: fits of six families, AIC, NPMLE gap test, model
##  confidence set, and the frontier n* / floor for the AIC winner against
##  the runner-up (at the winner's fitted parameters).
## =====================================================================
source("frontier_core.R"); source("datasets.R")
library(parallel)
CORES <- max(1, detectCores())
dir.create("results", showWarnings = FALSE)
B_GAP <- 199; B_MCS <- 2000

audit_one <- function(b) {
  f <- b$f; d <- tab_of(f); n <- sum(f); x <- rep(d$x, d$w)
  m <- mean(x); D <- var(x) * (n - 1) / n / m
  fits <- fit_all(f)
  ll <- sapply(fits, `[[`, "loglik")
  aic <- -2 * ll + 2 * sapply(names(fam), function(G) fam[[G]]$np)
  ord <- order(aic); best <- names(fam)[ord[1]]; second <- names(fam)[ord[2]]
  gaps <- t(sapply(names(fam), function(G) gap_test(f, G, B = B_GAP, seed = 1)))
  np <- npmle(f); pk <- npmle_peaks(np)
  ms <- mcs(f, B = B_MCS, seed = 1, fits = fits)
  pbest <- exp(fam[[best]]$lp(fits[[best]]$par, xs))
  fr <- frontier_pmf(pbest, second)
  ## spread of AIC among the four Poisson-Gold members
  pg <- c("PGD1", "PGD2", "PGD3", "PGD4")
  list(
    summary = data.frame(id = b$id, name = b$name, n = n, mean = m, disp = D,
      best = best, second = second, dAIC_second = sort(aic)[2] - min(aic),
      range_AIC_PGD = diff(range(aic[pg])), npmle_ll = np$loglik,
      gap_best = gaps[best, "gap"], gap_p_best = gaps[best, "p"],
      n_gap_reject = sum(gaps[, "p"] <= 0.05),
      mcs_size = sum(ms$in_set), n_star_best_vs_second = fr$n_star,
      floor_best_vs_second = fr$n_floor, claimed = b$claimed,
      npmle_peaks = paste(sprintf("%.2f (%.2f)", pk$lambda, pk$weight), collapse = "; ")),
    detail = data.frame(id = b$id, family = names(fam), loglik = ll, AIC = aic,
      dAIC = aic - min(aic), gap = gaps[, "gap"], gap_p = gaps[, "p"],
      MCS_p = ms$MCS_p, in_mcs = ms$in_set, row.names = NULL))
}

res <- mclapply(benchmarks, audit_one, mc.cores = CORES)
summ <- do.call(rbind, lapply(res, `[[`, "summary"))
det <- do.call(rbind, lapply(res, `[[`, "detail"))
saveRDS(list(summary = summ, detail = det), "results/T_audit.rds")
print(summ[, setdiff(names(summ), c("name", "claimed", "npmle_peaks"))], row.names = FALSE, digits = 4)
print(summ[, c("id", "claimed", "npmle_peaks")], row.names = FALSE)
print(det, row.names = FALSE, digits = 4)
