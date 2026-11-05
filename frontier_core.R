## =====================================================================
##  How many count distributions can a dataset tell apart?
##  Core functions (base R only)
##
##  Families: Poisson-Gold PGD(k, theta), k = 1..4, negative binomial (NB)
##  and zero-inflated Poisson (ZIP).
##    PGD(k = 1) = geometric (Poisson-exponential)
##    PGD(k = 2) = Poisson-Lindley      (Lindley parameter = 1/theta)
##    PGD(k = 3) = Poisson-Sujatha      (Sujatha parameter = 1/theta)
##    PGD(k = 4) = Poisson-Amarendra
##  (Hanandeh & Al-Nasser, 2021, Austrian Journal of Statistics, Table 1)
##
##  Contents
##   ldpgd, pgd_mean, pgd_disp, theta_for_mean   Poisson-Gold building blocks
##   fam, lab                                    family registry
##   project(w, G, x)       maximum likelihood / KL projection onto family G
##   project_hell(p, G)     Hellinger projection onto family G
##   frontier(F, G, mu)     KL, omega, Vuong frontier n*, universal floor
##   npmle(tab)             nonparametric MLE of the mixing law (EM on a grid)
##   gap_test(tab, G)       Proposition: NPMLE gap test (parametric bootstrap)
##   mcs(tab)               model confidence set (Hansen, Lunde & Nason, 2011)
##   vuong_z(tab, F, G)     Vuong statistic
## =====================================================================

options(width = 160)
alpha <- 0.05; beta <- 0.20            # level and 1 - power used for n*
z_a <- qnorm(1 - alpha / 2); z_b <- qnorm(1 - beta)
XMAX <- 400                            # support truncation for theoretical pmfs
xs <- 0:XMAX

## ---------- Poisson-Gold pmf (log scale, numerically stable) ------------
## P(X = x) = sum_j Gamma(x+j) q^(x+j) / ( Gamma(x+1) sum_j Gamma(j) theta^j ),
## q = theta / (1 + theta)
ldpgd <- function(x, theta, k) {
  q <- theta / (1 + theta)
  j <- seq_len(k)
  A <- outer(x, j, function(x, j) lgamma(x + j) + (x + j) * log(q))
  m <- apply(A, 1, max)
  num <- m + log(rowSums(exp(A - m)))
  B <- lgamma(j) + j * log(theta)
  den <- max(B) + log(sum(exp(B - max(B))))
  num - lgamma(x + 1) - den
}
pgd_mean <- function(theta, k) {
  j <- seq_len(k)
  sum(gamma(j + 1) * theta^j) / sum(gamma(j) * theta^(j - 1))
}
pgd_mixmom <- function(theta, k, r) {  # E[lambda^r] = r-th factorial moment of X
  j <- seq_len(k)
  sum(gamma(j + r) * theta^(j + r - 1)) / sum(gamma(j) * theta^(j - 1))
}
pgd_disp <- function(theta, k) {       # index of dispersion Var/mean
  m1 <- pgd_mixmom(theta, k, 1); m2 <- pgd_mixmom(theta, k, 2)
  1 + (m2 - m1^2) / m1
}
theta_for_mean <- function(mu, k)
  uniroot(function(t) pgd_mean(t, k) - mu, c(1e-8, 1e4), tol = 1e-12)$root

## ---------- family registry --------------------------------------------
## lp(p, x): log-pmf at counts x for unconstrained parameter vector p
fam <- list(
  PGD1 = list(np = 1, lp = function(p, x = xs) ldpgd(x, exp(p), 1)),
  PGD2 = list(np = 1, lp = function(p, x = xs) ldpgd(x, exp(p), 2)),
  PGD3 = list(np = 1, lp = function(p, x = xs) ldpgd(x, exp(p), 3)),
  PGD4 = list(np = 1, lp = function(p, x = xs) ldpgd(x, exp(p), 4)),
  NB   = list(np = 2, lp = function(p, x = xs)
                dnbinom(x, size = exp(p[1]), mu = exp(p[2]), log = TRUE)),
  ZIP  = list(np = 2, lp = function(p, x = xs) {
                pi0 <- plogis(p[1]); lam <- exp(p[2])
                out <- log1p(-pi0) + dpois(x, lam, log = TRUE)
                z <- x == 0
                out[z] <- log(pi0 + (1 - pi0) * exp(-lam)); out })
)
lab <- c(PGD1 = "Geometric (PGD k=1)", PGD2 = "Poisson-Lindley (PGD k=2)",
         PGD3 = "Poisson-Sujatha (PGD k=3)", PGD4 = "Poisson-Amarendra (PGD k=4)",
         NB = "Negative binomial", ZIP = "Zero-inflated Poisson")
short <- c(PGD1 = "Geometric", PGD2 = "P-Lindley", PGD3 = "P-Sujatha",
           PGD4 = "P-Amarendra", NB = "Neg. binomial", ZIP = "ZIP")

## ---------- Lindley-type (finite gamma-mixture) families ---------------
## Mixing density f(lambda; theta) proportional to
##   sum_j b_j(theta) lambda^(a_j - 1) exp(-theta lambda),   theta = rate,
## so the Poisson mixture is sum_j w_j(theta) NB(size a_j, prob theta/(1+theta)),
##   w_j proportional to b_j(theta) Gamma(a_j) theta^(-a_j).
gm_spec <- list(
  PL  = list(name = "Poisson-Lindley",   a = c(1, 2),       b = function(t) c(1, 1)),
  PSh = list(name = "Poisson-Shanker",   a = c(1, 2),       b = function(t) c(t, 1)),
  PAk = list(name = "Poisson-Akash",     a = c(1, 3),       b = function(t) c(1, 1)),
  PS  = list(name = "Poisson-Sujatha",   a = c(1, 2, 3),    b = function(t) c(1, 1, 1)),
  PA  = list(name = "Poisson-Amarendra", a = c(1, 2, 3, 4), b = function(t) c(1, 1, 1, 1)),
  PAr = list(name = "Poisson-Aradhana",  a = c(1, 2, 3),    b = function(t) c(1, 2, 1)),
  PIs = list(name = "Poisson-Ishita",    a = c(1, 3),       b = function(t) c(t, 1)),
  PSu = list(name = "Poisson-Suja",      a = c(1, 5),       b = function(t) c(1, 1))
)
gm_weights <- function(spec, t) {
  lw <- log(spec$b(t)) + lgamma(spec$a) - spec$a * log(t)
  w <- exp(lw - max(lw)); w / sum(w)
}
gm_mixmom <- function(spec, t, r)          # E[lambda^r]
  sum(gm_weights(spec, t) * exp(lgamma(spec$a + r) - lgamma(spec$a)) / t^r)
gm_mean <- function(spec, t) gm_mixmom(spec, t, 1)
gm_disp <- function(spec, t) { m1 <- gm_mixmom(spec, t, 1); 1 + (gm_mixmom(spec, t, 2) - m1^2) / m1 }
gm_rate_for_mean <- function(spec, mu)
  uniroot(function(lt) gm_mean(spec, exp(lt)) - mu, c(-20, 20), tol = 1e-13)$root
ldgm <- function(spec, x, t) {
  w <- gm_weights(spec, t); pr <- t / (1 + t)
  A <- sapply(seq_along(w), function(j) log(w[j]) + dnbinom(x, size = spec$a[j], prob = pr, log = TRUE))
  if (is.null(dim(A))) A <- matrix(A, nrow = 1)
  m <- apply(A, 1, max); m + log(rowSums(exp(A - m)))
}
gm_fam <- lapply(gm_spec, function(sp) list(np = 1, spec = sp,
  lp = local({ s <- sp; function(p, x = xs) ldgm(s, x, exp(p)) })))
fam_all <- c(fam, gm_fam,
  list(POIS = list(np = 1, lp = function(p, x = xs) dpois(x, exp(p), log = TRUE))))
lab_all <- c(lab, sapply(gm_spec, `[[`, "name"), POIS = "Poisson")

start_par <- function(G, w, x) {
  m <- sum(w * x) / sum(w); v <- sum(w * x^2) / sum(w) - m^2
  if (G == "NB") return(c(if (v > m) log(m^2 / (v - m)) else 5, log(max(m, 1e-3))))
  if (G == "ZIP") {
    lam <- max(v / m + m - 1, 1.01 * m)
    return(c(qlogis(min(max(1 - m / lam, 0.01), 0.99)), log(lam)))
  }
  NULL
}

## maximum likelihood (= KL projection when w is a pmf) of family G.
## w: weights (frequencies or probabilities) at counts x.
project <- function(w, G, x = xs) {
  f <- fam_all[[G]]
  keep <- w > 0; ww <- w[keep]; xx <- x[keep]
  obj <- function(p) -sum(ww * f$lp(p, xx))
  if (f$np == 1) {
    par <- optimize(obj, c(-14, 8), tol = 1e-10)$minimum
  } else {
    o <- optim(start_par(G, ww, xx), obj, method = "BFGS",
               control = list(reltol = 1e-12, maxit = 2000))
    par <- o$par
  }
  list(par = par, lp = f$lp(par, x), loglik = -obj(par))
}

## Hellinger projection: minimise H^2 = 1 - sum sqrt(p q)
project_hell <- function(p, G) {
  f <- fam_all[[G]]
  obj <- function(par) 1 - sum(sqrt(p * exp(f$lp(par))))
  start <- project(p, G)$par
  if (f$np == 1) {
    o <- optimize(obj, start + c(-1, 1), tol = 1e-12); list(par = o$minimum, H2 = o$objective)
  } else {
    o <- optim(start, obj, method = "BFGS", control = list(reltol = 1e-14, maxit = 2000))
    list(par = o$par, H2 = o$value)
  }
}

## ---------- frontier for one (true family F, rival family G, mean mu) ----
frontier <- function(F, G, mu) {
  stopifnot(grepl("^PGD", F))
  k <- as.integer(sub("PGD", "", F)); th <- theta_for_mean(mu, k)
  lpF <- ldpgd(xs, th, k); pF <- exp(lpF)
  pr <- project(pF, G)
  L <- lpF - pr$lp
  Delta <- sum(pF * L)
  omega <- sqrt(sum(pF * (L - Delta)^2))
  H2 <- project_hell(pF, G)$H2
  data.frame(true = F, rival = G, mean = mu, D_true = pgd_disp(th, k),
             KL = Delta, omega = omega, H2 = H2,
             n_star = ceiling(((z_a + z_b) * omega / Delta)^2),
             n_floor = ceiling((1 - alpha - beta)^2 / (2 * H2)),
             tail = 1 - sum(pF))
}

## same quantities for an arbitrary true pmf pF on xs (e.g. a fitted model)
frontier_pmf <- function(pF, G, alpha = 0.05, beta = 0.20) {
  z_a <- qnorm(1 - alpha / 2); z_b <- qnorm(1 - beta)
  pF <- pF / sum(pF); lpF <- log(pF)
  pr <- project(pF, G)
  ok <- pF > 0
  L <- lpF[ok] - pr$lp[ok]
  Delta <- sum(pF[ok] * L)
  omega <- sqrt(sum(pF[ok] * (L - Delta)^2))
  H2 <- project_hell(pF, G)$H2
  data.frame(KL = Delta, omega = omega, H2 = H2,
             n_star = ceiling(((z_a + z_b) * omega / Delta)^2),
             n_floor = ceiling((1 - alpha - beta)^2 / (2 * H2)))
}

## ---------- data helpers --------------------------------------------------
## a dataset is a frequency vector f, f[i] = number of observations equal to i-1
tab_of <- function(f) list(x = seq_along(f) - 1, w = f)
fit_all <- function(f, fams = names(fam)) {
  d <- tab_of(f)
  lapply(setNames(fams, fams), function(G) project(d$w, G, d$x))
}

## ---------- NPMLE of the mixing law (Lindsay, 1983) -----------------------
## Vertex-direction method with free atoms: alternately (a) optimise the
## locations and weights of a finite set of atoms, and (b) add an atom where
## Lindsay's gradient function
##     D(lambda) = (1/n) sum_x f_x Pois(x; lambda) / p(x) - 1
## is largest. By Lindsay (1983), the fit is the NPMLE iff D <= 0 for all
## lambda >= 0, and by concavity loglik(NPMLE) - loglik(fit) <= n * max D.
## The returned 'bound' is that certificate.
npmle <- function(f, tol = 1e-5, maxit = 60) {
  xv <- which(f > 0) - 1; cnt <- f[f > 0]; n <- sum(cnt)
  top <- 2 * max(xv) + 5
  lam_grid <- c(0, exp(seq(log(1e-4), log(top), length.out = 600)))
  Kg <- outer(xv, lam_grid, dpois)
  ## atoms: an optional fixed atom at 0 (weight w0) plus free atoms lam > 0
  pmf_of <- function(lam, w, w0) as.vector(outer(xv, lam, dpois) %*% w) + w0 * (xv == 0)
  fit_atoms <- function(lam, w, w0) {
    J <- length(lam); z <- w0 > 0
    wt <- c(if (z) w0 else NULL, w)
    par0 <- c(log(lam), log(wt[-1] / wt[1]))
    par0[!is.finite(par0)] <- -60; par0 <- pmin(pmax(par0, -60), 60)
    unpack <- function(par) {
      l <- exp(pmin(pmax(par[1:J], -30), 8)); a <- c(0, pmin(pmax(par[-(1:J)], -60), 60))
      ww <- exp(a - max(a)); ww <- ww / sum(ww)
      if (z) list(lam = l, w = ww[-1], w0 = ww[1]) else list(lam = l, w = ww, w0 = 0)
    }
    obj <- function(par) { u <- unpack(par); p <- pmf_of(u$lam, u$w, u$w0)
      -sum(cnt * log(pmax(p, 1e-300))) }
    o <- tryCatch(optim(par0, obj, method = "BFGS", control = list(reltol = 1e-12, maxit = 2000)),
                  error = function(e) optim(par0, obj, method = "Nelder-Mead",
                                            control = list(reltol = 1e-12, maxit = 20000)))
    unpack(o$par)
  }
  lam <- sum(xv * cnt) / n; w <- 1; w0 <- 0
  for (it in seq_len(maxit)) {
    u <- fit_atoms(lam, w, w0); lam <- u$lam; w <- u$w; w0 <- u$w0
    ## absorb atoms that collapsed to 0, merge near-duplicates, drop empty atoms
    small <- lam < 1e-4
    if (any(small)) { w0 <- w0 + sum(w[small]); lam <- lam[!small]; w <- w[!small] }
    if (length(lam) > 1) {
      o <- order(lam); lam <- lam[o]; w <- w[o]
      grp <- cumsum(c(1, diff(log(lam)) > 1e-3))
      lam <- as.vector(tapply(lam * w, grp, sum) / tapply(w, grp, sum)); w <- as.vector(tapply(w, grp, sum))
    }
    keep <- w > 1e-12; lam <- lam[keep]; w <- w[keep]
    if (!length(lam)) { lam <- 1e-3; w <- 1e-6 }
    s <- sum(w) + w0; w <- w / s; w0 <- w0 / s
    ## exact EM updates of the weights for fixed atoms (monotone in loglik)
    K0 <- cbind(xv == 0, outer(xv, lam, dpois)); ww <- c(w0, w)
    if (ww[1] == 0) ww[1] <- 1e-3
    ww <- ww / sum(ww)
    for (e in 1:3000) {
      pp <- as.vector(K0 %*% ww); wn <- ww * colSums(cnt * K0 / pp) / n
      if (max(abs(wn - ww)) < 1e-11) { ww <- wn; break }
      ww <- wn
    }
    w0 <- ww[1]; w <- ww[-1]
    p <- pmf_of(lam, w, w0)
    D <- colSums(cnt * Kg / p) / n - 1
    if (max(D) * n < tol) break
    j <- which.max(D)
    if (j == 1) { w0 <- w0 + 0.05; s <- sum(w) + w0; w <- w / s; w0 <- w0 / s }
    else { lam <- c(lam, lam_grid[j]); w <- c(w * 0.95, 0.05); w0 <- w0 * 0.95 }
  }
  p <- pmf_of(lam, w, w0)
  D <- colSums(cnt * Kg / p) / n - 1
  lam_all <- c(if (w0 > 0) 0 else NULL, lam); w_all <- c(if (w0 > 0) w0 else NULL, w)
  ord <- order(lam_all)
  list(loglik = sum(cnt * log(p)), lambda = lam_all[ord], weight = w_all[ord],
       bound = max(0, max(D)) * n, iter = it)
}
## atoms with weight >= min_w, merging atoms closer than 0.05
npmle_peaks <- function(fit, min_w = 0.02) {
  lam <- fit$lambda; w <- fit$weight
  grp <- cumsum(c(1, diff(lam) > 0.05))
  out <- data.frame(lambda = tapply(lam * w, grp, sum) / tapply(w, grp, sum),
                    weight = as.vector(tapply(w, grp, sum)))
  out[out$weight >= min_w, ]
}

## ---------- NPMLE gap test (parametric bootstrap) ---------------------------
## gap = loglik(NPMLE) - loglik(family G); large gap = G misfits
gap_test <- function(f, G, B = 200, seed = 1, tol = 1e-5) {
  n <- sum(f); d <- tab_of(f)
  fitG <- project(d$w, G, d$x)
  gap <- npmle(f, tol = tol)$loglik - fitG$loglik
  pG <- exp(fam_all[[G]]$lp(fitG$par, xs)); pG <- pG / sum(pG)
  set.seed(seed)
  gaps <- replicate(B, tryCatch({
    fb <- tabulate(sample(xs, n, TRUE, pG) + 1)
    db <- tab_of(fb)
    npmle(fb, tol = tol)$loglik - project(db$w, G, db$x)$loglik
  }, error = function(e) NA_real_))
  nb <- sum(!is.na(gaps)); gaps <- gaps[!is.na(gaps)]
  c(gap = gap, p = (1 + sum(gaps >= gap)) / (nb + 1), failed = B - nb)
}

## ---------- model confidence set (Hansen, Lunde & Nason, 2011), T_max -----
## Loss per observation = -log p(x_i) + (number of parameters)/n,
## i.e. AIC/(2n). Losses are computed at the fitted parameters and the
## block bootstrap of Hansen et al. reduces to the iid bootstrap here.
mcs <- function(f, fams = names(fam), B = 1000, alpha_mcs = 0.10, seed = 1, fits = NULL) {
  d <- tab_of(f); n <- sum(f); x <- rep(d$x, d$w)
  if (is.null(fits)) fits <- fit_all(f, fams)
  Lmat <- sapply(fams, function(G) -fits[[G]]$lp[x + 1] + fam_all[[G]]$np / n)
  ## bootstrap means of each column via multinomial resampling of the counts
  set.seed(seed)
  cnt <- rmultinom(B, n, d$w / n)                  # (#values) x B
  Lval <- sapply(fams, function(G) -fits[[G]]$lp + fam_all[[G]]$np / n)  # per value
  Lval[!is.finite(Lval)] <- 0                      # values never observed
  boot <- crossprod(cnt, Lval) / n                 # B x |fams|
  Lbar <- colMeans(Lmat)
  S <- fams; pvals <- setNames(rep(NA_real_, length(fams)), fams); running <- 0
  while (length(S) > 1) {
    d0 <- Lbar[S] - mean(Lbar[S])
    dB <- boot[, S, drop = FALSE] - rowMeans(boot[, S, drop = FALSE])
    dc <- sweep(dB, 2, d0)
    se <- sqrt(colMeans(dc^2))
    tt <- d0 / se
    Tstar <- apply(sweep(dc, 2, se, "/"), 1, max)
    p <- mean(Tstar >= max(tt))
    running <- max(running, p)
    worst <- S[which.max(tt)]
    pvals[worst] <- running
    S <- setdiff(S, worst)
  }
  pvals[S] <- 1
  data.frame(family = fams, MCS_p = pvals, in_set = pvals >= alpha_mcs, row.names = NULL)
}

## ---------- Vuong statistic (positive favours F) ----------------------------
vuong_z <- function(f, F, G, fits = NULL) {
  d <- tab_of(f); x <- rep(d$x, d$w); n <- length(x)
  if (is.null(fits)) fits <- fit_all(f, c(F, G))
  l <- fits[[F]]$lp[x + 1] - fits[[G]]$lp[x + 1]
  l <- l - (fam_all[[F]]$np - fam_all[[G]]$np) * log(n) / (2 * n)   # BIC-type correction
  sqrt(n) * mean(l) / sd(l)
}
