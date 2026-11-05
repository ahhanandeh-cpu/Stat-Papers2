## =====================================================================
##  Simulation study
##   S1  Vuong power: predicted (frontier formula) vs empirical, n = c * n*
##   S2  AIC selection: probability of choosing the true family vs n
##   S3  Model confidence set: coverage and size
##   S4  NPMLE gap test: size (truth inside the families) and power
##       (truth = two-point Poisson mixture, outside all families)
##  Usage: Rscript 02_simulation.R [S1|S2|S3|S4|all]
## =====================================================================
source("frontier_core.R")
library(parallel)
CORES <- max(1, detectCores())
stage <- commandArgs(TRUE); if (!length(stage)) stage <- "all"
dir.create("results", showWarnings = FALSE)
RNGkind("L'Ecuyer-CMRG")

pmf_pgd <- function(k, mu) { p <- exp(ldpgd(xs, theta_for_mean(mu, k), k)); p / sum(p) }
draw <- function(p, n) tabulate(sample(xs, n, TRUE, p) + 1)

## ---------------- S1: power of the Vuong test ----------------------------
if (stage %in% c("S1", "all")) {
  set.seed(101)
  cases <- data.frame(F = c("PGD3", "PGD4", "PGD2", "PGD2"),
                      G = c("PGD1", "PGD2", "PGD3", "NB"),
                      mu = c(2, 5, 2, 5))
  mult <- c(0.25, 0.5, 1, 2); R1 <- 1000
  out <- list()
  for (i in seq_len(nrow(cases))) {
    fr <- frontier(cases$F[i], cases$G[i], cases$mu[i])
    k <- as.integer(sub("PGD", "", cases$F[i])); p <- pmf_pgd(k, cases$mu[i])
    for (m in mult) {
      n <- ceiling(m * fr$n_star)
      z <- unlist(mclapply(seq_len(R1), function(r) {
        f <- draw(p, n); d <- tab_of(f); x <- rep(d$x, d$w)
        fits <- list(project(d$w, cases$F[i], d$x), project(d$w, cases$G[i], d$x))
        l <- fits[[1]]$lp[x + 1] - fits[[2]]$lp[x + 1]
        sqrt(n) * mean(l) / sd(l)
      }, mc.cores = CORES, mc.set.seed = TRUE))
      out[[length(out) + 1]] <- data.frame(
        F = cases$F[i], G = cases$G[i], mu = cases$mu[i], n_star = fr$n_star,
        mult = m, n = n,
        power_pred = pnorm(sqrt(n) * fr$KL / fr$omega - z_a),
        power_emp = mean(z > z_a), mc_se = sqrt(mean(z > z_a) * (1 - mean(z > z_a)) / R1))
      cat(sprintf("S1 %s vs %s mu=%g n=%d: predicted %.3f empirical %.3f\n", cases$F[i],
                  cases$G[i], cases$mu[i], n, tail(out, 1)[[1]]$power_pred, tail(out, 1)[[1]]$power_emp))
    }
  }
  saveRDS(do.call(rbind, out), "results/S1_power.rds")
}

## ---------------- S2: AIC selection ----------------------------------------
if (stage %in% c("S2", "all")) {
  set.seed(202)
  cand <- c("PGD1", "PGD2", "PGD3", "PGD4", "NB")
  grid <- expand.grid(k = 2:4, mu = c(0.5, 1, 2, 5), n = c(100, 300, 1000, 3000, 10000))
  R2 <- 500
  res <- mclapply(seq_len(nrow(grid)), function(g) {
    k <- grid$k[g]; mu <- grid$mu[g]; n <- grid$n[g]; p <- pmf_pgd(k, mu)
    truth <- paste0("PGD", k)
    sel <- character(R2); close <- logical(R2)
    for (r in seq_len(R2)) {
      f <- draw(p, n); fits <- fit_all(f, cand)
      aic <- sapply(cand, function(G) -2 * fits[[G]]$loglik + 2 * fam[[G]]$np)
      sel[r] <- cand[which.min(aic)]
      close[r] <- min(aic[cand != truth]) - aic[truth] < 2   # a rival within 2 AIC units
    }
    data.frame(k = k, mu = mu, n = n, p_correct = mean(sel == truth),
               p_rival_within2 = mean(close),
               t(sapply(cand, function(G) mean(sel == G))))
  }, mc.cores = CORES, mc.set.seed = TRUE)
  res <- do.call(rbind, res)
  saveRDS(res, "results/S2_aic.rds")
  print(res[, 1:5], row.names = FALSE)
}

## ---------------- S3: model confidence set -------------------------------------
two_point <- function() {    # NPMLE of the summer thunderstorm data (two peaks)
  lam <- c(0.16, 1.61); w <- c(0.56, 0.44)
  p <- w[1] * dpois(xs, lam[1]) + w[2] * dpois(xs, lam[2]); p / sum(p)
}
if (stage %in% c("S3", "all")) {
  set.seed(303)
  grid <- expand.grid(truth = c("PGD2", "TWO"), mu = c(0.75, 2), n = c(300, 1000, 5000),
                      stringsAsFactors = FALSE)
  grid <- grid[!(grid$truth == "TWO" & grid$mu == 2), ]
  R3 <- 300
  res <- mclapply(seq_len(nrow(grid)), function(g) {
    p <- if (grid$truth[g] == "TWO") two_point() else pmf_pgd(2, grid$mu[g])
    n <- grid$n[g]
    sz <- cover <- numeric(R3)
    for (r in seq_len(R3)) {
      m <- mcs(draw(p, n), B = 500, seed = r)
      sz[r] <- sum(m$in_set)
      cover[r] <- if (grid$truth[g] == "TWO") NA else m$in_set[m$family == "PGD2"]
    }
    data.frame(truth = grid$truth[g], mu = if (grid$truth[g] == "TWO") sum(xs * p) else grid$mu[g],
               n = n, coverage = mean(cover), mean_size = mean(sz), p_all6 = mean(sz == 6))
  }, mc.cores = CORES, mc.set.seed = TRUE)
  res <- do.call(rbind, res)
  saveRDS(res, "results/S3_mcs.rds")
  print(res, row.names = FALSE)
}

## ---------------- S4: NPMLE gap test --------------------------------------------
if (stage %in% c("S4", "all")) {
  set.seed(404)
  grid <- expand.grid(truth = c("PGD2", "TWO"), n = c(300, 1000), stringsAsFactors = FALSE)
  R4 <- 100; B4 <- 59
  res <- list()
  for (g in seq_len(nrow(grid))) {
    p <- if (grid$truth[g] == "TWO") two_point() else pmf_pgd(2, 0.8)
    n <- grid$n[g]
    pv <- unlist(mclapply(seq_len(R4), function(r)
      gap_test(draw(p, n), "PGD2", B = B4, seed = r, tol = 1e-3)[["p"]],
      mc.cores = CORES, mc.set.seed = TRUE))
    res[[g]] <- data.frame(truth = grid$truth[g], n = n, reject_05 = mean(pv <= 0.05),
                           reject_10 = mean(pv <= 0.10))
    cat(sprintf("S4 truth=%s n=%d: rejection rate at 5%% = %.3f\n", grid$truth[g], n, res[[g]]$reject_05))
    saveRDS(do.call(rbind, res), "results/S4_gap.rds")
  }
  saveRDS(do.call(rbind, res), "results/S4_gap.rds")
}
cat("done", stage, "\n")
