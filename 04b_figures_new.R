## =====================================================================
##  Additional figures (revision): AIC selection with seven candidates,
##  published claims, and Lindley-type lifetime laws.
##  Colours: validated categorical slots in fixed order; every series also
##  carries its own line type or symbol for greyscale print.
## =====================================================================
dir.create("figures", showWarnings = FALSE)
col6 <- c("#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300")
ink <- "#0b0b0b"; ink2 <- "#52514e"; gridc <- "#dcdbd5"
lty6 <- c(1, 2, 4, 5, 6, 3); pch4 <- c(16, 17, 15, 18)
setpar <- function(...) par(family = "sans", col.axis = ink2, col.lab = ink, fg = ink2,
                            las = 1, mgp = c(2.4, 0.6, 0), tcl = -0.3, cex.axis = 0.85, ...)

## ---- Figure 4 (revised): AIC selection among seven candidates ----
s2 <- readRDS("results/S2b_aic.rds")
nm <- c(PL = "Poisson-Lindley", PSh = "Poisson-Shanker", PAk = "Poisson-Akash",
        PS = "Poisson-Sujatha", PA = "Poisson-Amarendra")
mus <- c(0.5, 1, 2, 5)
pdf("figures/fig4_aic.pdf", width = 7.2, height = 4.8)
setpar(mfrow = c(2, 3), mar = c(3.4, 3.4, 1.8, 0.4))
for (tr in names(nm)) {
  plot(NA, xlim = c(100, 10000), ylim = c(0, 1), log = "x", xlab = "Sample size n",
       ylab = "P(AIC picks the truth)", main = paste("Truth:", nm[tr]),
       cex.main = 0.85, font.main = 1, col.main = ink, xaxt = "n")
  axis(1, at = c(100, 300, 1000, 3000, 10000), labels = c("100", "300", "1k", "3k", "10k"))
  abline(h = seq(0, 1, 0.2), v = c(100, 300, 1000, 3000, 10000), col = gridc, lwd = 0.6)
  abline(h = 1 / 7, lty = 3, col = ink2)
  for (j in seq_along(mus)) {
    d <- s2[s2$truth == tr & s2$mu == mus[j], ]; d <- d[order(d$n), ]
    lines(d$n, d$p_correct, col = col6[j], lty = lty6[j], lwd = 2)
    points(d$n, d$p_correct, col = col6[j], pch = pch4[j], cex = 0.9)
  }
}
plot.new()
legend("center", bty = "n", cex = 0.9, text.col = ink, title = "Mean of the truth", title.col = ink,
       col = col6[1:4], lty = lty6[1:4], pch = pch4, lwd = 2, legend = mus)
legend("bottom", bty = "n", cex = 0.8, text.col = ink2, lty = 3, col = ink2,
       legend = "random choice among seven (1/7)")
dev.off()

## ---- Figure 5: published claims ----
rv <- readRDS("results/T_review.rds")$claims
rv <- rv[!is.na(rv$dAIC), ]
rv$lab <- paste0(rv$paper, ": ", rv$data)
rv <- rv[order(rv$dAIC), ]
k <- nrow(rv)
pdf("figures/fig5_claims.pdf", width = 6.8, height = 7.2)
setpar(mar = c(3.6, 13.5, 2.6, 1))
xl <- c(-15, 8)
plot(NA, xlim = xl, ylim = c(0.5, k + 0.5), yaxt = "n", ylab = "",
     xlab = "AIC of best rival minus AIC of claimed winner")
rect(0, 0, 2, k + 1, col = adjustcolor("#9a9890", 0.20), border = NA)
abline(v = seq(-15, 5, 5), col = gridc, lwd = 0.6); abline(v = c(0, 2), col = ink2, lty = c(1, 2))
recomputed <- !is.na(rv$n_star)
x <- pmax(pmin(rv$dAIC, xl[2]), xl[1])
segments(0, seq_len(k), x, seq_len(k), col = gridc)
points(x, seq_len(k), pch = ifelse(recomputed, 19, 21), bg = "white", col = ink, cex = 0.9)
clipped <- rv$dAIC < xl[1]
if (any(clipped)) text(xl[1], which(clipped), sprintf("%.1f", rv$dAIC[clipped]), pos = 4, cex = 0.65, col = ink2)
axis(2, at = seq_len(k), labels = rv$lab, cex.axis = 0.62, tick = FALSE)
mtext("0 to 2:", side = 3, line = 1.0, at = 1, cex = 0.62, col = ink2); mtext("not decisive", side = 3, line = 0.3, at = 1, cex = 0.62, col = ink2)
mtext("a rival fits better", side = 3, line = 0.3, at = -7.5, cex = 0.62, col = ink2)
mtext("winner ahead", side = 3, line = 0.3, at = 6, cex = 0.62, col = ink2)
legend("topleft", bty = "n", cex = 0.7, text.col = ink, pch = c(19, 21), pt.bg = "white",
       legend = c("recomputed here", "as printed in the source"))
dev.off()

## ---- Figure 6: lifetime laws, frontier against the mean ----
cl <- readRDS("results/F6_lifetime_curves.rds")
rn <- c(PS = "Sujatha", PSh = "Shanker", PAk = "Akash", EXP = "exponential", GAM = "gamma", WEI = "Weibull")
pdf("figures/fig6_lifetime.pdf", width = 6.6, height = 4.2)
setpar(mar = c(5.2, 4.6, 0.8, 0.8))
plot(NA, xlim = c(0.05, 1000), ylim = c(10, 1e7), log = "xy", yaxt = "n", xaxt = "n",
     xlab = "", ylab = "")
axis(2, at = 10^(1:7), labels = parse(text = paste0("10^", 1:7)))
axis(1, at = c(0.1, 1, 10, 100, 1000), labels = c("0.1", "1", "10", "100", "1000"))
mtext("Mean lifetime, in the unit of measurement", side = 1, line = 2.2, col = ink)
mtext("n* for 80% power (Lindley true)", side = 2, line = 3.3, col = ink, las = 0)
abline(h = 10^(1:7), v = c(0.1, 1, 10, 100, 1000), col = gridc, lwd = 0.6)
i <- 0
for (G in names(rn)) {
  i <- i + 1; d <- cl[cl$rival == G, ]
  lines(d$mean, d$n_star, col = col6[i], lty = lty6[i], lwd = 2)
}
legend("topright", bty = "n", cex = 0.75, text.col = ink, lwd = 2, col = col6, lty = lty6,
       legend = paste("Lindley vs", rn), bg = "white")
um <- c("bank, hours" = 9.877 / 60, "bank, minutes" = 9.877, "bank, seconds" = 592.62)
axis(1, at = um, labels = FALSE, tcl = 0.5, col.ticks = ink)
text(um, 13, names(um), srt = 90, adj = 0, cex = 0.62, col = ink2)
dev.off()
cat("new figures written\n")
