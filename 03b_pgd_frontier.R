## =====================================================================
##  For each benchmark dataset: frontier n* and floor for the best-fitting
##  one-parameter Poisson-Gold member against (a) the runner-up member and
##  (b) the negative binomial, at the fitted parameters of the best member.
## =====================================================================
source("frontier_core.R"); source("datasets.R")
pg <- c("PGD1", "PGD2", "PGD3", "PGD4")
out <- do.call(rbind, lapply(benchmarks, function(b) {
  fits <- fit_all(b$f, c(pg, "NB"))
  ll <- sapply(fits[pg], `[[`, "loglik"); o <- order(-ll)
  best <- pg[o[1]]; second <- pg[o[2]]
  pb <- exp(fam[[best]]$lp(fits[[best]]$par, xs))
  a <- frontier_pmf(pb, second); nb <- frontier_pmf(pb, "NB")
  data.frame(id = b$id, n = sum(b$f), best_pgd = best, second_pgd = second,
             dll = ll[o[1]] - ll[o[2]],
             n_star_pgd = a$n_star, floor_pgd = a$n_floor,
             n_star_nb = nb$n_star, floor_nb = nb$n_floor)
}))
saveRDS(out, "results/T_pgd_frontier.rds")
print(out, row.names = FALSE, digits = 4)
cat("n below floor (PGD pair):", sum(out$n < out$floor_pgd), "of", nrow(out), "\n")
cat("n below n* (PGD pair):", sum(out$n < out$n_star_pgd), "of", nrow(out), "\n")
