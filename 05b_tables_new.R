## =====================================================================
##  LaTeX tables for the revision, generated from saved results
## =====================================================================
dir.create("tables", showWarnings = FALSE)
fmtn <- function(x) ifelse(is.na(x), "--", formatC(x, format = "d", big.mark = "{,}"))
w <- function(lines, file) writeLines(lines, file.path("tables", file))
lbls <- c(PGD1 = "Geometric", PL = "P--Lindley", PSh = "P--Shanker", PAk = "P--Akash",
          PS = "P--Sujatha", PA = "P--Amarendra", PAr = "P--Aradhana", PIs = "P--Ishita",
          PSu = "P--Suja", NB = "Neg.\\ binomial", POIS = "Poisson", ZIP = "ZIP",
          PGD2 = "PGD $k=2$", PGD3 = "PGD $k=3$", PGD4 = "PGD $k=4$")

## ---- kappa table ----
kap <- data.frame(
  fam = c("Geometric (exponential mixing)", "Poisson--Lindley", "Poisson--Shanker", "Poisson--Akash",
          "Poisson--Sujatha", "Poisson--Amarendra", "Poisson--Aradhana", "Poisson--Ishita"),
  dens = c("$1$", "$1+\\lambda$", "$\\theta+\\lambda$", "$1+\\lambda^2$", "$1+\\lambda+\\lambda^2$",
           "$1+\\lambda+\\lambda^2+\\lambda^3$", "$(1+\\lambda)^2$", "$\\theta+\\lambda^2$"),
  r1 = c("0", "1", "0", "0", "1", "1", "2", "0"),
  r2 = c("0", "0", "0", "1", "1", "1", "1", "0"),
  kappa = c(0, -2, 0, 4, 2, 2, -4, 0))
num <- readRDS("results/T_kappa.rds")
rows <- sprintf("%s & %s & %s & %s & $%d$ \\\\", kap$fam, kap$dens, kap$r1, kap$r2, kap$kappa)
w(c("\\begin{tabular}{llccr}", "\\toprule",
    "Family & Mixing density $\\propto p_\\theta(\\lambda)\\,e^{-\\theta\\lambda}$, $p_\\theta(\\lambda)=$ & $\\rho_1$ & $\\rho_2$ & $\\kappa$\\\\",
    "\\midrule", rows, "\\bottomrule", "\\end{tabular}"), "tab_kappa.tex")

## ---- verification of the general law: summary ----
g <- readRDS("results/T_general_summary.rds")
cat(sprintf("ratio at 0.005: %.3f-%.3f; at 0.05: %.3f-%.3f; H2/KL %.3f-%.3f; om2/KL %.3f-%.3f; reg max %.1f\n",
            min(g$ratio_KL_005), max(g$ratio_KL_005), min(g$ratio_KL_05), max(g$ratio_KL_05),
            min(g$H2_over_KL), max(g$H2_over_KL), min(g$om2_over_KL), max(g$om2_over_KL), max(g$reg_max)))
## predicted vs numerical constants of the dispersion gap
kv <- c(PL = -2, PSh = 0, PAk = 4, PS = 2, PA = 2, PAr = -4, PIs = 0)
g$pred <- kv[g$F] - kv[g$G]
g3 <- g[g$order == 3, ]
cat(sprintf("order-3 pairs: %d; max |numerical - predicted| = %.3f\n", nrow(g3), max(abs(g3$const - g3$pred))))
saveRDS(g, "results/T_general_summary2.rds")

## ---- lookup table (appendix) ----
lk <- readRDS("results/T_lookup.rds")
fams <- c("PGD1", "PL", "PSh", "PAk", "PS", "PA")
prs <- list()
for (i in 1:5) for (j in (i + 1):6) prs[[length(prs) + 1]] <- c(fams[i], fams[j])
for (f in fams) prs[[length(prs) + 1]] <- c(f, "NB")
ms <- c(0.25, 0.5, 1, 2, 5)
rows <- sapply(prs, function(p) {
  if (p[1] == "PGD1" && p[2] == "NB")
    return(sprintf("%s & %s & \\multicolumn{5}{c}{nested: not separable at any $n$} \\\\", lbls[p[1]], lbls[p[2]]))
  v <- sapply(ms, function(m) { r <- lk[lk$true == p[1] & lk$rival == p[2] & lk$mean == m, ]
    if (r$KL < 1e-12) "$\\infty$" else fmtn(r$n_star) })
  sprintf("%s & %s & %s \\\\", lbls[p[1]], lbls[p[2]], paste(v, collapse = " & "))
})
w(c("\\begin{tabular}{llrrrrr}", "\\toprule",
    "True & Rival & $\\mu=0.25$ & $\\mu=0.5$ & $\\mu=1$ & $\\mu=2$ & $\\mu=5$\\\\", "\\midrule",
    rows[1:15], "\\addlinespace", rows[16:21], "\\bottomrule", "\\end{tabular}"), "tab_lookup.tex")
lk2 <- lk[!(lk$true == "PGD1" & lk$rival == "NB") & lk$KL > 1e-12, ]
r14 <- lk2$n_star / lk2$floor
cat(sprintf("n*/floor over the lookup table: median %.2f, range %.2f-%.2f\n", median(r14), min(r14), max(r14)))

## ---- review summary by paper ----
rv <- readRDS("results/T_review.rds"); P <- rv$papers; C <- rv$claims
rows <- sapply(seq_len(nrow(P)), function(i) {
  d <- C[C$paper == P$id[i], ]; e <- d[!is.na(d$dAIC), ]
  sprintf("%s & %s & %d & %s & %d & %d & %s \\\\", P$authors[i], P$proposed[i], nrow(d),
          if (nrow(e)) sprintf("%.2f", median(e$dAIC)) else "--",
          sum(e$dAIC >= 0 & e$dAIC < 2), sum(e$dAIC < 0),
          if (any(!is.na(d$floor))) sprintf("%d/%d", sum(d$n < d$floor, na.rm = TRUE), sum(!is.na(d$floor))) else "--")
})
w(c("\\begin{tabular}{p{4.6cm}p{3.4cm}rrrrr}", "\\toprule",
    "Source & Model proposed & Claims & Median $\\Delta$AIC & $0\\le\\Delta<2$ & $\\Delta<0$ & Below floor\\\\",
    "\\midrule", rows, "\\midrule",
    sprintf("Total & & %d & %.2f & %d & %d & %d/%d \\\\", nrow(C), median(C$dAIC, na.rm = TRUE),
            sum(C$dAIC >= 0 & C$dAIC < 2, na.rm = TRUE), sum(C$dAIC < 0, na.rm = TRUE),
            sum(C$n < C$floor, na.rm = TRUE), sum(!is.na(C$floor))),
    "\\bottomrule", "\\end{tabular}"), "tab_review.tex")

## ---- lifetime: units ----
u <- readRDS("results/T_lifetime_units.rds")
ln <- c(EXP = "Exponential", PL = "Lindley", PSh = "Shanker", PAk = "Akash", PS = "Sujatha",
        PA = "Amarendra", GAM = "Gamma", WEI = "Weibull")
rows <- sapply(seq_len(nrow(u)), function(i) with(u[i, ], sprintf(
  "%s & %s & %.3g & %s & %s & %.2f & %.2f & %s & %s \\\\",
  if (unit %in% c("minutes", "months")) (if (data == "bank") "Bank waiting times" else "Bladder cancer remission") else "",
  unit, mean, ln[best], paste0(ln[best_1p], "/", ln[second_1p]), dAIC_1p, range_AIC_1p,
  fmtn(n_star), fmtn(floor))))
w(c("\\begin{tabular}{llrllrrrr}", "\\toprule",
    "Data & Unit & Mean & Best overall & Best/second Lindley-type & $\\Delta$AIC & AIC range & $n^*$ & Floor\\\\",
    "\\midrule", rows[1:3], "\\addlinespace", rows[4:6], "\\bottomrule", "\\end{tabular}"), "tab_lifetime.tex")

## ---- S4b gap test (if available) ----
if (file.exists("results/S4b_gap.rds")) {
  s4 <- readRDS("results/S4b_gap.rds")
  rows <- sapply(seq_len(nrow(s4)), function(i) with(s4[i, ], sprintf(
    "%s & %s & %d & %.3f & %.3f \\\\", if (truth == "TWO") "Two-point mixture (outside all families)" else
      "Poisson--Lindley, mean 0.8 (inside)", fmtn(n), reps, reject_05, reject_10)))
  w(c("\\begin{tabular}{lrrcc}", "\\toprule", "Truth & $n$ & Replicates & Rejection at 5\\% & Rejection at 10\\%\\\\",
      "\\midrule", rows, "\\bottomrule", "\\end{tabular}"), "tab_sim_gap.tex")
}
cat("new tables written\n")
