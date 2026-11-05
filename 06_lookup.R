## =====================================================================
##  Lookup table for authors and referees:
##  frontier n* (Vuong, alpha = 0.05, power 0.80) and universal floor for
##  every ordered pair of seven Lindley-type families, at eight means.
##  Also: the closed-form small-mean coefficient kappa of each family.
## =====================================================================
source("frontier_core.R")
dir.create("results", showWarnings = FALSE)
truths <- c("PGD1", "PL", "PSh", "PAk", "PS", "PA")
rivals <- c("PGD1", "PL", "PSh", "PAk", "PS", "PA", "NB")
means <- c(0.25, 0.5, 0.75, 1, 1.5, 2, 3, 5)
pmf_family <- function(G, mu) {
  if (G == "PGD1") p <- exp(ldpgd(xs, theta_for_mean(mu, 1), 1))
  else p <- exp(ldgm(gm_spec[[G]], xs, exp(gm_rate_for_mean(gm_spec[[G]], mu))))
  p / sum(p)
}
rows <- list()
for (F in truths) for (mu in means) {
  pF <- pmf_family(F, mu)
  for (G in setdiff(rivals, F)) {
    fr <- frontier_pmf(pF, G)
    rows[[length(rows) + 1]] <- data.frame(true = F, rival = G, mean = mu,
                                           n_star = fr$n_star, floor = fr$n_floor, KL = fr$KL)
  }
}
look <- do.call(rbind, rows)
saveRDS(look, "results/T_lookup.rds")
write.csv(look, "results/lookup_table.csv", row.names = FALSE)

## kappa = 4 c2/c0 - 2 (c1/c0)^2 from the weights of the Gamma(2), Gamma(3) parts:
## kappa = 2 * omega3 - 2 * omega2^2, omega2 = lim w2 / s, omega3 = lim w3 / s^2
kappa_num <- function(G) {
  if (G == "PGD1") return(0)
  sp <- gm_spec[[G]]; t <- 1e5; s <- 1 / t
  w <- gm_weights(sp, t); w2 <- sum(w[sp$a == 2]); w3 <- sum(w[sp$a == 3])
  2 * w3 / s^2 - 2 * (w2 / s)^2
}
kap <- sapply(truths, kappa_num)
saveRDS(kap, "results/T_kappa.rds"); print(round(kap, 3))
print(head(look[order(look$true, look$rival, look$mean), ], 20), row.names = FALSE)
