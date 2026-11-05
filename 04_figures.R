## =====================================================================
##  Figures for the manuscript (PDF, base graphics)
##  Colours: validated categorical slots 1-4 (blue, orange, aqua, yellow);
##  every series also has its own line type / symbol for greyscale print.
## =====================================================================
source("frontier_core.R"); source("datasets.R")
dir.create("figures", showWarnings = FALSE)
col4 <- c("#2a78d6", "#eb6834", "#1baf7a", "#eda100")
ink <- "#0b0b0b"; ink2 <- "#52514e"; gridc <- "#dcdbd5"
lty4 <- c(1, 2, 4, 5); pch4 <- c(16, 17, 15, 18)
setpar <- function(...) par(family = "sans", col.axis = ink2, col.lab = ink, fg = ink2,
                            las = 1, mgp = c(2.4, 0.6, 0), tcl = -0.3, cex.axis = 0.85, ...)
Dk <- function(m, k) pgd_disp(theta_for_mean(m, k), k)
dat <- t(sapply(benchmarks, function(b) { x <- rep(seq_along(b$f) - 1, b$f); n <- length(x)
  c(mean = mean(x), disp = var(x) * (n - 1) / n / mean(x), n = n) }))
cross <- readRDS("results/T2_crossings.rds")

## ---- Figure 1: mean-dispersion map ----
pdf("figures/fig1_mean_dispersion.pdf", width = 7.2, height = 3.6)
setpar(mfrow = c(1, 2), mar = c(3.6, 3.8, 2, 0.8))
for (panel in 1:2) {
  xl <- if (panel == 1) c(0, 5) else c(0.35, 1.25)
  yl <- if (panel == 1) c(1, 6.2) else c(1.1, 2.3)
  mg <- seq(max(0.01, xl[1]), xl[2], length.out = 300)
  plot(NA, xlim = xl, ylim = yl, xlab = "Mean", ylab = "Index of dispersion",
       main = if (panel == 1) "(a) Four one-parameter curves" else "(b) Region of the benchmark data",
       cex.main = 0.95, font.main = 1, col.main = ink)
  abline(h = pretty(yl), v = pretty(xl), col = gridc, lwd = 0.6)
  for (k in 1:4) lines(mg, sapply(mg, Dk, k = k), col = col4[k], lty = lty4[k], lwd = 2)
  if (panel == 1) {
    legend("topleft", bty = "n", lwd = 2, col = col4, lty = lty4, cex = 0.8, text.col = ink,
           legend = c("Geometric (k=1)", "Poisson-Lindley (k=2)", "Poisson-Sujatha (k=3)",
                      "Poisson-Amarendra (k=4)"))
  } else {
    cc <- cross[cross$k1 >= 2, ]
    points(cc$mu_cross, cc$D_cross, pch = 21, bg = "white", col = ink, cex = 1.1)
    points(dat[, "mean"], dat[, "disp"], pch = 19, col = ink, cex = 0.8)
    lp <- ifelse(rownames(dat) %in% c("accidents", "yeast"), 2, 4)
    text(dat[, "mean"], dat[, "disp"], seq_len(nrow(dat)), pos = lp, cex = 0.7, col = ink, offset = 0.3)
  }
}
dev.off()

## ---- Figure 2: frontier as a function of the mean ----
curves <- readRDS("results/F2_curves.rds")
pdf("figures/fig2_frontier.pdf", width = 6.4, height = 4)
setpar(mar = c(3.6, 4.6, 0.8, 0.8))
prs <- list(c("PGD2", "PGD3"), c("PGD2", "PGD4"), c("PGD3", "PGD4"), c("PGD2", "NB"))
nm <- c("P-Lindley vs P-Sujatha", "P-Lindley vs P-Amarendra", "P-Sujatha vs P-Amarendra",
        "P-Lindley vs neg. binomial")
plot(NA, xlim = c(0.1, 5), ylim = c(1e2, 1e9), log = "xy", xlab = "Mean of the true distribution",
     ylab = "", yaxt = "n")
at <- 10^(2:9); axis(2, at = at, labels = parse(text = paste0("10^", 2:9)))
mtext("Sample size n* for 80% power", side = 2, line = 3.3, col = ink, las = 0)
abline(h = at, v = c(0.1, 0.2, 0.5, 1, 2, 5), col = gridc, lwd = 0.6)
rect(0.1, 80, 5, 1012, col = adjustcolor("#9a9890", 0.18), border = NA)
text(0.105, 300, "n of the ten\nbenchmark datasets", adj = c(0, 0.5), cex = 0.72, col = ink2)
for (i in 1:4) {
  cdat <- curves[curves$true == prs[[i]][1] & curves$rival == prs[[i]][2], ]
  lines(cdat$mean, cdat$n_star, col = col4[i], lty = lty4[i], lwd = 2)
}
cc <- cross[cross$k1 >= 2, ]
abline(v = cc$mu_cross, lty = 3, col = ink2, lwd = 0.8)
text(cc$mu_cross, 4e8, c("2|3", "2|4", "3|4"), cex = 0.7, col = ink2, pos = 4, offset = 0.15)
legend("topright", bty = "n", lwd = 2, col = col4, lty = lty4, legend = nm, cex = 0.78,
       text.col = ink, bg = "white")
dev.off()

## ---- Figure 3: predicted vs empirical power ----
s1 <- readRDS("results/S1_power.rds")
pdf("figures/fig3_power.pdf", width = 7.2, height = 2.6)
setpar(mfrow = c(1, 4), mar = c(3.4, 3.2, 2.2, 0.4), oma = c(0, 1, 0, 0))
cs <- unique(s1[, c("F", "G", "mu")])
ttl <- c(P = "")
lbl <- c(PGD1 = "geometric", PGD2 = "P-Lindley", PGD3 = "P-Sujatha", PGD4 = "P-Amarendra", NB = "neg. binomial")
for (i in seq_len(nrow(cs))) {
  d <- s1[s1$F == cs$F[i] & s1$G == cs$G[i] & s1$mu == cs$mu[i], ]
  plot(NA, xlim = c(0.2, 2.4), ylim = c(0, 1), log = "x", xlab = "n / n*", ylab = if (i == 1) "Power" else "",
       main = sprintf("%s vs %s\nmean %g, n* = %s", lbl[cs$F[i]], lbl[cs$G[i]], cs$mu[i],
                      format(d$n_star[1], big.mark = ",")), cex.main = 0.78, font.main = 1, col.main = ink,
       xaxt = "n")
  axis(1, at = c(0.25, 0.5, 1, 2), labels = c("1/4", "1/2", "1", "2"))
  abline(h = seq(0, 1, 0.2), v = c(0.25, 0.5, 1, 2), col = gridc, lwd = 0.6)
  mm <- exp(seq(log(0.2), log(2.4), length.out = 100))
  fr <- frontier(cs$F[i], cs$G[i], cs$mu[i])
  lines(mm, pnorm(sqrt(mm * fr$n_star) * fr$KL / fr$omega - z_a), col = col4[1], lwd = 2)
  arrows(d$mult, d$power_emp - 2 * d$mc_se, d$mult, d$power_emp + 2 * d$mc_se, angle = 90,
         code = 3, length = 0.02, col = ink)
  points(d$mult, d$power_emp, pch = 21, bg = col4[2], col = ink, cex = 1.1)
  if (i == 1) legend("bottomright", bty = "n", cex = 0.72, text.col = ink,
                     legend = c("Predicted", "Simulated"), lwd = c(2, NA), pch = c(NA, 21),
                     col = c(col4[1], ink), pt.bg = c(NA, col4[2]))
}
dev.off()

## ---- Figure 4: AIC selects the true family ----
s2 <- readRDS("results/S2_aic.rds")
pdf("figures/fig4_aic.pdf", width = 7.2, height = 2.8)
setpar(mfrow = c(1, 3), mar = c(3.4, 3.4, 1.8, 0.4))
kn <- c("2" = "Truth: Poisson-Lindley", "3" = "Truth: Poisson-Sujatha", "4" = "Truth: Poisson-Amarendra")
mus <- c(0.5, 1, 2, 5)
for (k in 2:4) {
  plot(NA, xlim = c(100, 10000), ylim = c(0, 1), log = "x", xlab = "Sample size n",
       ylab = if (k == 2) "P(AIC picks the true family)" else "", main = kn[as.character(k)],
       cex.main = 0.85, font.main = 1, col.main = ink, xaxt = "n")
  axis(1, at = c(100, 300, 1000, 3000, 10000), labels = c("100", "300", "1k", "3k", "10k"))
  abline(h = seq(0, 1, 0.2), v = c(100, 300, 1000, 3000, 10000), col = gridc, lwd = 0.6)
  abline(h = 0.2, lty = 3, col = ink2)
  for (j in seq_along(mus)) {
    d <- s2[s2$k == k & s2$mu == mus[j], ]; d <- d[order(d$n), ]
    lines(d$n, d$p_correct, col = col4[j], lty = lty4[j], lwd = 2)
    points(d$n, d$p_correct, col = col4[j], pch = pch4[j], cex = 0.9)
  }
  if (k == 2) {
    legend("topleft", bty = "n", cex = 0.72, text.col = ink, title = "Mean", title.col = ink,
           col = col4, lty = lty4, pch = pch4, lwd = 2, legend = mus)
    text(2000, 0.2, "chance (1 in 5)", pos = 1, cex = 0.68, col = ink2, offset = 0.3)
  }
}
dev.off()
cat("figures written\n")
