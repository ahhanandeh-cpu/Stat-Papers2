## =====================================================================
##  Frontier tables and curves
##   - frontier n* and universal floor for every pair, means 0.5, 1, 2, 5
##   - dispersion-curve crossings
##   - n* as a function of the mean (curves for Figure 2)
##   - matched-mean KL vs the KL-closest rival (shifting the mean barely helps)
## =====================================================================
source("frontier_core.R")
dir.create("results", showWarnings = FALSE)

## ---- 1. frontier table ----
means <- c(0.5, 1, 2, 5)
pairs <- rbind(c("PGD2", "PGD1"), c("PGD3", "PGD2"), c("PGD2", "PGD3"), c("PGD4", "PGD3"),
               c("PGD3", "PGD4"), c("PGD2", "PGD4"), c("PGD4", "PGD2"),
               c("PGD2", "NB"), c("PGD3", "NB"), c("PGD4", "NB"))
tab1 <- do.call(rbind, lapply(seq_len(nrow(pairs)), function(i)
  do.call(rbind, lapply(means, function(m) frontier(pairs[i, 1], pairs[i, 2], m)))))
stopifnot(max(tab1$tail) < 1e-10)
saveRDS(tab1, "results/T1_frontier.rds")
print(tab1[, c("true", "rival", "mean", "D_true", "KL", "n_star", "n_floor")], row.names = FALSE)

## ---- 2. crossings of dispersion curves ----
Dk <- function(m, k) pgd_disp(theta_for_mean(m, k), k)
cross <- do.call(rbind, lapply(list(c(1, 2), c(1, 3), c(1, 4), c(2, 3), c(2, 4), c(3, 4)), function(pr) {
  f <- function(m) Dk(m, pr[1]) - Dk(m, pr[2])
  g <- seq(0.02, 20, by = 0.01); v <- sapply(g, f); s <- which(diff(sign(v)) != 0)
  root <- if (length(s)) uniroot(f, g[c(s[1], s[1] + 1)], tol = 1e-12)$root else NA
  th <- if (is.na(root)) c(NA, NA) else c(theta_for_mean(root, pr[1]), theta_for_mean(root, pr[2]))
  data.frame(k1 = pr[1], k2 = pr[2], mu_cross = root, D_cross = if (is.na(root)) NA else Dk(root, pr[1]),
             theta1 = th[1], theta2 = th[2], n_roots = length(s))
}))
saveRDS(cross, "results/T2_crossings.rds"); print(cross, row.names = FALSE)

## ---- 3. n* and floor as functions of the mean ----
mg <- exp(seq(log(0.1), log(5), length.out = 120))
curves <- do.call(rbind, lapply(list(c("PGD2", "PGD3"), c("PGD2", "PGD4"), c("PGD3", "PGD4"),
                                     c("PGD2", "NB")), function(pr)
  do.call(rbind, lapply(mg, function(m) frontier(pr[1], pr[2], m)))))
saveRDS(curves, "results/F2_curves.rds")

## ---- 4. matched-mean KL vs KL-closest rival ----
klm <- do.call(rbind, lapply(list(c(2, 3), c(2, 4), c(3, 4), c(1, 2)), function(pr)
  do.call(rbind, lapply(c(0.3, 0.5, 1, 2, 5), function(m) {
    pa <- exp(ldpgd(xs, theta_for_mean(m, pr[1]), pr[1]))
    pb <- exp(ldpgd(xs, theta_for_mean(m, pr[2]), pr[2]))
    ok <- pa > 1e-300
    data.frame(k1 = pr[1], k2 = pr[2], mean = m,
               KL_matched = sum(pa[ok] * (log(pa[ok]) - log(pb[ok]))),
               KL_closest = frontier(paste0("PGD", pr[1]), paste0("PGD", pr[2]), m)$KL,
               dD2_over_8 = (Dk(m, pr[1]) - Dk(m, pr[2]))^2 / 8)
  }))))
klm$ratio <- klm$KL_closest / klm$KL_matched
saveRDS(klm, "results/T3_matched.rds"); print(klm, row.names = FALSE)
