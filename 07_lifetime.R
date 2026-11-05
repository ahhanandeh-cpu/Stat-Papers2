## =====================================================================
##  Lindley-type LIFETIME distributions
##  Density f(x; theta) proportional to sum_j b_j(theta) x^(a_j - 1) exp(-theta x)
##  (the mixing laws of frontier_core.R used directly as lifetime laws):
##  exponential, Lindley, Shanker, Akash, Sujatha, Amarendra; rivals gamma and
##  Weibull. Frontier and floor by numerical integration; effect of the unit
##  of measurement on two classical datasets.
## =====================================================================
source("frontier_core.R")
dir.create("results", showWarnings = FALSE)
lt_spec <- c(list(EXP = list(name = "Exponential", a = 1, b = function(t) 1)),
             gm_spec[c("PL", "PSh", "PAk", "PS", "PA")])
lt_name <- c(EXP = "Exponential", PL = "Lindley", PSh = "Shanker", PAk = "Akash",
             PS = "Sujatha", PA = "Amarendra", GAM = "Gamma", WEI = "Weibull")
ld_lt <- function(G, x, par) {
  if (G == "GAM") return(dgamma(x, shape = exp(par[1]), rate = exp(par[2]), log = TRUE))
  if (G == "WEI") return(dweibull(x, shape = exp(par[1]), scale = exp(par[2]), log = TRUE))
  sp <- lt_spec[[G]]; t <- exp(par); w <- gm_weights(sp, t)
  A <- sapply(seq_along(w), function(j) log(w[j]) + dgamma(x, shape = sp$a[j], rate = t, log = TRUE))
  if (is.null(dim(A))) A <- matrix(A, nrow = length(x))
  m <- apply(A, 1, max); m + log(rowSums(exp(A - m)))
}
np_lt <- function(G) if (G %in% c("GAM", "WEI")) 2 else 1
lt_mean <- function(G, t) if (G == "EXP") 1 / t else gm_mean(lt_spec[[G]], t)
lt_rate_for_mean <- function(G, mu) uniroot(function(lt) lt_mean(G, exp(lt)) - mu, c(-25, 25), tol = 1e-13)$root

## maximum likelihood for data x
fit_lt <- function(x, G) {
  obj <- function(p) -sum(ld_lt(G, x, p))
  if (np_lt(G) == 1) { o <- optimize(obj, c(-20, 20), tol = 1e-12); list(par = o$minimum, loglik = -o$objective) }
  else {
    m <- mean(x); v <- var(x)
    p0 <- if (G == "GAM") log(c(m^2 / v, m / v)) else c(log(1.2), log(m))
    o <- optim(p0, obj, method = "BFGS", control = list(reltol = 1e-14, maxit = 5000)); list(par = o$par, loglik = -o$value)
  }
}

## divergences by numerical integration on a fine grid (log-spaced)
grid_for <- function(mu) { u <- exp(seq(log(mu * 1e-6), log(mu * 80), length.out = 6000)); u }
integ <- function(u, y) sum(diff(u) * (head(y, -1) + tail(y, -1)) / 2)
frontier_lt <- function(F, parF, G, mu) {
  u <- grid_for(mu); lf <- ld_lt(F, u, parF); f <- exp(lf)
  obj <- function(p) -integ(u, f * ld_lt(G, u, p))
  if (np_lt(G) == 1) p <- optimize(obj, c(-20, 20), tol = 1e-12)$minimum
  else {
    v <- integ(u, f * u^2) - mu^2
    p0 <- if (G == "GAM") log(c(mu^2 / v, mu / v)) else c(log(1.2), log(mu))
    p <- optim(p0, obj, method = "BFGS", control = list(reltol = 1e-14, maxit = 5000))$par
  }
  L <- lf - ld_lt(G, u, p)
  KL <- integ(u, f * L); om2 <- integ(u, f * (L - KL)^2)
  hobj <- function(q) 0.5 * integ(u, (sqrt(f) - exp(ld_lt(G, u, q) / 2))^2)
  H2 <- if (np_lt(G) == 1) optimize(hobj, p + c(-1, 1), tol = 1e-12)$objective
        else optim(p, hobj, method = "BFGS", control = list(reltol = 1e-14))$value
  data.frame(KL = KL, omega2 = om2, H2 = H2,
             n_star = ceiling((z_a + z_b)^2 * om2 / KL^2),
             floor = ceiling((1 - alpha - beta)^2 / (2 * H2)))
}

## ---- 1. frontier against the mean (unit-dependence) ----
pairs <- list(c("PL", "PS"), c("PL", "PSh"), c("PL", "PAk"), c("PL", "EXP"), c("PL", "GAM"), c("PL", "WEI"))
mus <- c(0.1, 0.3, 1, 3, 10, 30, 100)
cur <- do.call(rbind, lapply(pairs, function(pr) do.call(rbind, lapply(mus, function(mu) {
  t <- lt_rate_for_mean(pr[1], mu)
  cbind(data.frame(true = pr[1], rival = pr[2], mean = mu), frontier_lt(pr[1], t, pr[2], mu))
}))))
saveRDS(cur, "results/T_lifetime_frontier.rds")
print(cur[, c("true", "rival", "mean", "KL", "n_star", "floor")], row.names = FALSE, digits = 4)

## ---- 2. two datasets in three units ----
bank <- c(0.8, 0.8, 1.3, 1.5, 1.8, 1.9, 1.9, 2.1, 2.6, 2.7, 2.9, 3.1, 3.2, 3.3, 3.5, 3.6, 4, 4.1, 4.2, 4.2,
          4.3, 4.3, 4.4, 4.4, 4.6, 4.7, 4.7, 4.8, 4.9, 4.9, 5.0, 5.3, 5.5, 5.7, 5.7, 6.1, 6.2, 6.2, 6.2, 6.3,
          6.7, 6.9, 7.1, 7.1, 7.1, 7.1, 7.4, 7.6, 7.7, 8, 8.2, 8.6, 8.6, 8.6, 8.8, 8.8, 8.9, 8.9, 9.5, 9.6,
          9.7, 9.8, 10.7, 10.9, 11.0, 11.0, 11.1, 11.2, 11.2, 11.5, 11.9, 12.4, 12.5, 12.9, 13.0, 13.1, 13.3,
          13.6, 13.7, 13.9, 14.1, 15.4, 15.4, 17.3, 17.3, 18.1, 18.2, 18.4, 18.9, 19.0, 19.9, 20.6, 21.3,
          21.4, 21.9, 23, 27, 31.6, 33.1, 38.5)
bladder <- c(0.08, 2.09, 3.48, 4.87, 6.94, 8.66, 13.11, 23.63, 0.2, 2.23, 0.26, 0.31, 0.73,
             0.52, 4.98, 6.97, 9.02, 13.29, 0.4, 2.26, 3.57, 5.06, 7.09, 11.98, 4.51, 2.07,
             0.22, 13.8, 25.74, 0.5, 2.46, 3.64, 5.09, 7.26, 9.47, 14.24, 19.13, 6.54, 3.36,
             0.82, 0.51, 2.54, 3.7, 5.17, 7.28, 9.74, 14.76, 26.31, 0.81, 1.76, 8.53, 6.93,
             0.62, 3.82, 5.32, 7.32, 10.06, 14.77, 32.15, 2.64, 3.88, 5.32, 3.25, 12.03, 8.65,
             0.39, 10.34, 14.83, 34.26, 0.9, 2.69, 4.18, 5.34, 7.59, 10.66, 4.5, 20.28, 12.63,
             0.96, 36.66, 1.05, 2.69, 4.23, 5.41, 7.62, 10.75, 16.62, 43.01, 6.25, 2.02, 22.69,
             0.19, 2.75, 4.26, 5.41, 7.63, 17.12, 46.12, 1.26, 2.83, 4.33, 8.37, 3.36, 5.49,
             0.66, 11.25, 17.14, 79.05, 1.35, 2.87, 5.62, 7.87, 11.64, 17.36, 12.02, 6.76,
             0.4, 3.02, 4.34, 5.71, 7.93, 11.79, 18.1, 1.46, 4.4, 5.85, 2.02, 12.07)
stopifnot(length(bank) == 100, length(bladder) == 128)
cands <- c("EXP", "PL", "PSh", "PAk", "PS", "PA", "GAM", "WEI")
one_param <- c("EXP", "PL", "PSh", "PAk", "PS", "PA")
units <- list(bank = list(x = bank, units = c(minutes = 1, hours = 1 / 60, seconds = 60)),
              bladder = list(x = bladder, units = c(months = 1, years = 1 / 12, days = 30.4375)))
rows <- list()
for (ds in names(units)) for (un in names(units[[ds]]$units)) {
  x <- units[[ds]]$x * units[[ds]]$units[[un]]
  fits <- lapply(setNames(cands, cands), function(G) fit_lt(x, G))
  aic <- sapply(cands, function(G) -2 * fits[[G]]$loglik + 2 * np_lt(G))
  ## frontier between the best and second-best one-parameter Lindley-type laws
  o1 <- one_param[order(aic[one_param])]
  fr <- frontier_lt(o1[1], fits[[o1[1]]]$par, o1[2], mean(x))
  rows[[length(rows) + 1]] <- data.frame(data = ds, unit = un, n = length(x), mean = mean(x),
    best = names(which.min(aic)), best_1p = o1[1], second_1p = o1[2],
    dAIC_1p = aic[o1[2]] - aic[o1[1]], range_AIC_1p = diff(range(aic[one_param])),
    n_star = fr$n_star, floor = fr$floor, t(round(aic - min(aic), 2)))
}
uni <- do.call(rbind, rows)
saveRDS(uni, "results/T_lifetime_units.rds")
print(uni, row.names = FALSE, digits = 4)
cat("Lindley log-likelihood, bank data in minutes:", fit_lt(bank, "PL")$loglik, "\n")

## ---- 3. curves of n* against the mean, for Figure 6 ----
mg <- exp(seq(log(0.05), log(1000), length.out = 40))
curves_lt <- do.call(rbind, lapply(pairs, function(pr) do.call(rbind, lapply(mg, function(mu) {
  t <- lt_rate_for_mean(pr[1], mu)
  cbind(data.frame(true = pr[1], rival = pr[2], mean = mu), frontier_lt(pr[1], t, pr[2], mu))
}))))
saveRDS(curves_lt, "results/F6_lifetime_curves.rds")
cat("curves done\n")
