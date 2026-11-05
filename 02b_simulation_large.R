## =====================================================================
##  Larger simulation study (supersedes S2 and S4 of 02_simulation.R)
##   S2b  AIC selection among seven candidates, five true families,
##        1,000 replicates per cell
##   S4b  NPMLE gap test, 200 replicates x 99 bootstrap samples
##  Usage: Rscript 02b_simulation_large.R [S2b|S4b|all]
## =====================================================================
source("frontier_core.R")
library(parallel)
CORES <- max(1, detectCores())
stage <- commandArgs(TRUE); if (!length(stage)) stage <- "all"
dir.create("results", showWarnings = FALSE)
RNGkind("L'Ecuyer-CMRG")

pmf_family <- function(G, mu) {
  if (G == "PGD1") p <- exp(ldpgd(xs, theta_for_mean(mu, 1), 1))
  else p <- exp(ldgm(gm_spec[[G]], xs, exp(gm_rate_for_mean(gm_spec[[G]], mu))))
  p / sum(p)
}
draw <- function(p, n) tabulate(sample(xs, n, TRUE, p) + 1)
cand <- c("PGD1", "PL", "PSh", "PAk", "PS", "PA", "NB")

if (stage %in% c("S2b", "all")) {
  set.seed(2202)
  grid <- expand.grid(truth = c("PL", "PSh", "PAk", "PS", "PA"), mu = c(0.5, 1, 2, 5),
                      n = c(100, 300, 1000, 3000, 10000), stringsAsFactors = FALSE)
  R <- 1000
  res <- mclapply(seq_len(nrow(grid)), function(g) {
    p <- pmf_family(grid$truth[g], grid$mu[g]); n <- grid$n[g]; tr <- grid$truth[g]
    sel <- character(R); close <- logical(R)
    for (r in seq_len(R)) {
      fits <- fit_all(draw(p, n), cand)
      aic <- sapply(cand, function(G) -2 * fits[[G]]$loglik + 2 * fam_all[[G]]$np)
      sel[r] <- cand[which.min(aic)]
      close[r] <- min(aic[cand != tr]) - aic[tr] < 2
    }
    data.frame(truth = tr, mu = grid$mu[g], n = n, p_correct = mean(sel == tr),
               p_rival_within2 = mean(close), t(sapply(cand, function(G) mean(sel == G))))
  }, mc.cores = CORES, mc.set.seed = TRUE)
  res <- do.call(rbind, res)
  saveRDS(res, "results/S2b_aic.rds")
  print(res[, 1:5], row.names = FALSE)
}

two_point <- function() {
  lam <- c(0.16, 1.61); w <- c(0.56, 0.44)
  p <- w[1] * dpois(xs, lam[1]) + w[2] * dpois(xs, lam[2]); p / sum(p)
}
if (stage %in% c("S4b", "all")) {
  set.seed(4404)
  grid <- expand.grid(truth = c("PGD2", "TWO"), n = c(300, 1000), stringsAsFactors = FALSE)
  R4 <- 200; B4 <- 99
  res <- list()
  for (g in seq_len(nrow(grid))) {
    p <- if (grid$truth[g] == "TWO") two_point() else { q <- exp(ldpgd(xs, theta_for_mean(0.8, 2), 2)); q / sum(q) }
    n <- grid$n[g]
    pv <- unlist(mclapply(seq_len(R4), function(r)
      gap_test(draw(p, n), "PGD2", B = B4, seed = r, tol = 1e-3)[["p"]],
      mc.cores = CORES, mc.set.seed = TRUE))
    res[[g]] <- data.frame(truth = grid$truth[g], n = n, reps = sum(!is.na(pv)),
                           reject_05 = mean(pv <= 0.05, na.rm = TRUE),
                           reject_10 = mean(pv <= 0.10, na.rm = TRUE))
    cat(sprintf("S4b truth=%s n=%d: rejection at 5%% = %.3f\n", grid$truth[g], n, res[[g]]$reject_05))
    saveRDS(do.call(rbind, res), "results/S4b_gap.rds")
  }
}
cat("done", stage, "\n")
