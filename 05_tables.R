## =====================================================================
##  LaTeX tables generated from saved results (no hand-copied numbers)
## =====================================================================
dir.create("tables", showWarnings = FALSE)
source("datasets.R")
fmtn <- function(x) formatC(x, format = "d", big.mark = "{,}")
lbl <- c(PGD1 = "Geometric", PGD2 = "Poisson--Lindley", PGD3 = "Poisson--Sujatha",
         PGD4 = "Poisson--Amarendra", NB = "Negative binomial", ZIP = "Zero-inflated Poisson")
lbls <- c(PGD1 = "Geom.", PGD2 = "PL", PGD3 = "PS", PGD4 = "PA", NB = "NB", ZIP = "ZIP")
w <- function(lines, file) writeLines(lines, file.path("tables", file))

## ---- Table: crossings ----
cr <- readRDS("results/T2_crossings.rds")
rows <- apply(cr, 1, function(r) {
  pair <- sprintf("%s ($k=%d$) and %s ($k=%d$)", lbl[paste0("PGD", r[["k1"]])], as.integer(r[["k1"]]),
                  lbl[paste0("PGD", r[["k2"]])], as.integer(r[["k2"]]))
  if (is.na(r[["mu_cross"]])) sprintf("%s & \\multicolumn{4}{c}{no crossing} \\\\", pair) else
    sprintf("%s & %.4f & %.4f & %.3f & %.3f \\\\", pair, r[["mu_cross"]], r[["D_cross"]],
            r[["theta1"]], r[["theta2"]])
})
w(c("\\begin{tabular}{lcccc}", "\\toprule",
    "Pair & $\\mu_c$ & $D(\\mu_c)$ & $\\theta_k(\\mu_c)$ & $\\theta_{k'}(\\mu_c)$\\\\",
    "\\midrule", rows, "\\bottomrule", "\\end{tabular}"), "tab_crossings.tex")

## ---- Table: frontier ----
t1 <- readRDS("results/T1_frontier.rds")
pairs <- unique(t1[, c("true", "rival")])
keep <- with(pairs, paste(true, rival)) %in% c("PGD2 PGD1", "PGD2 PGD3", "PGD3 PGD2", "PGD2 PGD4",
                                              "PGD3 PGD4", "PGD4 PGD3", "PGD2 NB", "PGD4 NB")
pairs <- pairs[keep, ]
rows <- apply(pairs, 1, function(p) {
  d <- t1[t1$true == p[["true"]] & t1$rival == p[["rival"]], ]; d <- d[order(d$mean), ]
  paste0(lbl[p[["true"]]], " & ", lbl[p[["rival"]]], " & ",
         paste(sprintf("%s & %s", fmtn(d$n_star), fmtn(d$n_floor)), collapse = " & "), " \\\\")
})
w(c("\\begin{tabular}{llrrrrrrrr}", "\\toprule",
    "& & \\multicolumn{2}{c}{$\\mu=0.5$} & \\multicolumn{2}{c}{$\\mu=1$} & \\multicolumn{2}{c}{$\\mu=2$} & \\multicolumn{2}{c}{$\\mu=5$}\\\\",
    "\\cmidrule(lr){3-4}\\cmidrule(lr){5-6}\\cmidrule(lr){7-8}\\cmidrule(lr){9-10}",
    "True family & Rival family & $n^*$ & floor & $n^*$ & floor & $n^*$ & floor & $n^*$ & floor\\\\",
    "\\midrule", rows, "\\bottomrule", "\\end{tabular}"), "tab_frontier.tex")

## ---- Table: simulation S1 ----
s1 <- readRDS("results/S1_power.rds")
rows <- sapply(seq_len(nrow(s1)), function(i) with(s1[i, ], sprintf(
  "%s & %s & %g & %s & %s & %s & %.3f & %.3f (%.3f) \\\\",
  if (mult == 0.25) lbl[F] else "", if (mult == 0.25) lbl[G] else "", if (mult == 0.25) mu else NA,
  if (mult == 0.25) fmtn(n_star) else "", c("1/4", "1/2", "1", "2")[match(mult, c(0.25, 0.5, 1, 2))],
  fmtn(n), power_pred, power_emp, mc_se)))
rows <- gsub(" & NA & ", " &  & ", rows)
idx <- which(s1$mult == 0.25)[-1]
rows2 <- character(0)
for (i in seq_along(rows)) { if (i %in% idx) rows2 <- c(rows2, "\\addlinespace"); rows2 <- c(rows2, rows[i]) }
w(c("\\begin{tabular}{llrrcrcc}", "\\toprule",
    "True & Rival & $\\mu$ & $n^*$ & $n/n^*$ & $n$ & Predicted & Simulated (s.e.)\\\\", "\\midrule",
    rows2, "\\bottomrule", "\\end{tabular}"), "tab_sim_power.tex")

## ---- Table: simulation S3 (model confidence set) ----
s3 <- readRDS("results/S3_mcs.rds"); s3 <- s3[order(s3$truth, s3$mu, s3$n), ]
rows <- sapply(seq_len(nrow(s3)), function(i) with(s3[i, ], sprintf(
  "%s & %.2f & %s & %s & %.2f & %.2f \\\\",
  if (truth == "TWO") "Two-point mixture" else "Poisson--Lindley", mu, fmtn(n),
  if (is.na(coverage)) "--" else sprintf("%.3f", coverage), mean_size, p_all6)))
w(c("\\begin{tabular}{lrrccc}", "\\toprule",
    "Truth & Mean & $n$ & Coverage & Mean set size & P(all six kept)\\\\", "\\midrule",
    rows, "\\bottomrule", "\\end{tabular}"), "tab_sim_mcs.tex")

## ---- Table: simulation S4 (gap test) ----
if (file.exists("results/S4_gap.rds")) {
  s4 <- readRDS("results/S4_gap.rds")
  rows <- sapply(seq_len(nrow(s4)), function(i) with(s4[i, ], sprintf(
    "%s & %s & %.3f & %.3f \\\\", if (truth == "TWO") "Two-point mixture (outside all families)" else
      "Poisson--Lindley, mean 0.8 (inside)", fmtn(n), reject_05, reject_10)))
  w(c("\\begin{tabular}{lrcc}", "\\toprule", "Truth & $n$ & Rejection at 5\\% & Rejection at 10\\%\\\\",
      "\\midrule", rows, "\\bottomrule", "\\end{tabular}"), "tab_sim_gap.tex")
}

## ---- Table: audit of benchmark data ----
if (file.exists("results/T_audit.rds")) {
  a <- readRDS("results/T_audit.rds"); s <- a$summary; det <- a$detail
  s$claimed <- sapply(s$id, function(i) benchmarks[[i]]$claimed)
  short_names <- c(yeast = "Yeast cells", mites = "Red mites", cornborer = "Corn borer",
                   chromatid = "Chromatid aberrations", lesions60 = "Lesions, 60 $\\mu$g/kg",
                   lesions70 = "Lesions, 70 $\\mu$g/kg", lesions90 = "Lesions, 90 $\\mu$g/kg",
                   accidents = "Accidents", thunder_jun = "Thunderstorms, June",
                   thunder_sum = "Thunderstorms, summer")
  s$name <- short_names[s$id]
  pf <- readRDS("results/T_pgd_frontier.rds")
  rows <- sapply(seq_len(nrow(s)), function(i) with(s[i, ], {
    q <- pf[pf$id == id, ]
    sprintf("%d & %s & %s & %.2f & %.2f & %d & %s & %s/%s & %s & %s \\\\",
            i, name, fmtn(n), mean, disp, mcs_size, sprintf("%d/6", n_gap_reject),
            lbls[q$best_pgd], lbls[q$second_pgd], fmtn(q$n_star_pgd), fmtn(q$floor_pgd))
  }))
  w(c("\\begin{tabular}{rlrrrccl rr}", "\\toprule",
      "& & & & & \\multicolumn{2}{c}{Audit} & \\multicolumn{3}{c}{Best Poisson--Gold member vs runner-up}\\\\",
      "\\cmidrule(lr){6-7}\\cmidrule(lr){8-10}",
      "& Dataset & $n$ & Mean & $D$ & MCS & Gap rej. & Pair & $n^*$ & Floor\\\\",
      "\\midrule", rows, "\\bottomrule", "\\end{tabular}"), "tab_audit.tex")
  ## detail: Delta AIC per family per dataset
  fams <- c("PGD1", "PGD2", "PGD3", "PGD4", "NB", "ZIP")
  rows <- sapply(seq_len(nrow(s)), function(i) {
    dd <- det[det$id == s$id[i], ]; dd <- dd[match(fams, dd$family), ]
    cells <- sprintf("%.2f%s", dd$dAIC, ifelse(dd$gap_p <= 0.05, "$^{\\dagger}$", ""))
    paste0(i, " & ", paste(cells, collapse = " & "), " & ", s$claimed[i], " \\\\")
  })
  w(c("\\begin{tabular}{rrrrrrrl}", "\\toprule",
      paste0("& ", paste(lbls[fams], collapse = " & "), " & Favoured in the source\\\\"),
      "\\midrule", rows, "\\bottomrule", "\\end{tabular}"), "tab_audit_detail.tex")
}
cat("tables written\n")
