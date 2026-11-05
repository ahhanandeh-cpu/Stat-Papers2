## Numerical checks of the general small-mean theorem for Lindley-type families
## (i)  KL(P||Q) / {(D_F - D_G)^2 / 8} -> 1 as mu -> 0 at equal means
## (ii) H^2 / KL -> 1/4, omega^2 / KL -> 2
## (iii) moment-gap regularity: max_r |delta_r| / (r! mu^(r-2) |delta_2|), r = 3..12
source("frontier_core.R")
fams <- names(gm_spec)
pmf_gm <- function(nm, mu) { p <- exp(ldgm(gm_spec[[nm]], xs, exp(gm_rate_for_mean(gm_spec[[nm]], mu)))); p / sum(p) }
mom <- function(nm, mu, r) gm_mixmom(gm_spec[[nm]], exp(gm_rate_for_mean(gm_spec[[nm]], mu)), r)
disp <- function(nm, mu) gm_disp(gm_spec[[nm]], exp(gm_rate_for_mean(gm_spec[[nm]], mu)))
out <- list()
for (i in 1:(length(fams) - 1)) for (j in (i + 1):length(fams)) for (mu in c(0.005, 0.01, 0.02, 0.05, 0.1)) {
  a <- fams[i]; b <- fams[j]
  P <- pmf_gm(a, mu); Q <- pmf_gm(b, mu); ok <- P > 1e-300 & Q > 1e-300
  ## cancellation-free forms: KL = sum P {eps - log1p(eps)}, eps = (Q - P)/P,
  ## H^2 = (1/2) sum (P - Q)^2 / (sqrt P + sqrt Q)^2
  eps <- (Q[ok] - P[ok]) / P[ok]
  KL <- sum(P[ok] * ifelse(abs(eps) < 1e-4, eps^2 / 2 - eps^3 / 3 + eps^4 / 4, eps - log1p(eps)))
  H2 <- 0.5 * sum((P[ok] - Q[ok])^2 / (sqrt(P[ok]) + sqrt(Q[ok]))^2)
  L <- -log1p(eps); om2 <- sum(P[ok] * (L - KL)^2)
  dD <- disp(a, mu) - disp(b, mu)
  d2 <- mom(a, mu, 2) - mom(b, mu, 2)
  reg <- max(sapply(3:12, function(r) abs(mom(a, mu, r) - mom(b, mu, r)) / (factorial(r) * mu^(r - 2) * abs(d2))))
  out[[length(out) + 1]] <- data.frame(F = a, G = b, mu = mu, dD = dD, KL = KL,
    ratio_KL = KL / (dD^2 / 8), H2_over_KL = H2 / KL, om2_over_KL = om2 / KL, reg = reg)
}
res <- do.call(rbind, out)
saveRDS(res, "results/T_general_theory.rds")
## order and constant of the dispersion gap, from mu = 0.005 and 0.01
summ <- do.call(rbind, lapply(split(res, paste(res$F, res$G)), function(d) {
  a <- d[d$mu == 0.005, ]; b <- d[d$mu == 0.01, ]
  s <- log(abs(b$dD) / abs(a$dD)) / log(2)
  data.frame(F = a$F, G = a$G, order = round(s), const = a$dD / 0.005^round(s),
             ratio_KL_005 = a$ratio_KL, ratio_KL_05 = d$ratio_KL[d$mu == 0.05],
             H2_over_KL = a$H2_over_KL, om2_over_KL = a$om2_over_KL, reg_max = max(d$reg))
}))
saveRDS(summ, "results/T_general_summary.rds")
print(summ, row.names = FALSE, digits = 3)
