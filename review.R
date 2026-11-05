## =====================================================================
##  Review of published model comparisons (open-access sources)
##  Each record: a dataset and a published claim that one model fits best.
##  For each claim we report
##   - dAIC: AIC of the best rival minus AIC of the claimed winner
##           (recomputed by us when both models are implemented here;
##            otherwise taken from the AIC or -2 log L printed in the source)
##   - n* and the universal floor for the claimed winner against the best
##     rival, when both are implemented here (at the winner's fitted values)
## =====================================================================
source("frontier_core.R")

## ---------- sources ----------
papers <- data.frame(rbind(
  c("SF16", "Shanker and Fesshaye (2016)", "Biometrics & Biostatistics International Journal 3(4): 00069", "Poisson-Sujatha"),
  c("SF17", "Shanker, Fesshaye, Shanker, Leonida and Sium (2017)", "Biometrics & Biostatistics International Journal 5(1): 6-14", "Poisson-Shanker"),
  c("TS25", "Tabassum and Shanker (2025)", "Hungarian Statistical Review 8(2): 94-110", "Poisson-Sujatha"),
  c("SS19", "Shanker and Shukla (2018)", "Journal of Applied Quantitative Methods 13(2)", "Poisson-weighted Akash"),
  c("SS20", "Shanker and Shukla (2020)", "Hungarian Statistical Review 3(1): 3-17", "quasi Poisson-Aradhana"),
  c("RVP24", "Radhakrishnan, Venkatesan and Prasanth (2024)", "The Bioscan 19(2): 928-935", "Poisson-Suja"),
  c("MP25", "Morifat and Para (2025)", "Reliability: Theory & Applications 20(3): 390-405", "two-parameter Poisson-Garima"),
  c("HE26", "Haj Ahmad and El-Awady (2026)", "AIMS Mathematics 11(9): 28121-28160", "Poisson two-sum XLindley"),
  c("ATBA24", "Alomair, Tajuddin, Bakouch and Almohisen (2024)", "Axioms 13(2): 125", "Poisson-improved second-degree Lindley"),
  c("MBO20", "Mazucheli, Bertoli and Oliveira (2020)", "Revista Colombiana de Estadistica 43(1): 21-48", "two discrete Shanker models"),
  c("DPA26", "Divya, Prasanth and Muhammed Anvar (2026)", "Statistics and Applications 24(1): 127-153", "zero-one-inflated Poisson-Garima"),
  c("AS26", "Ahmad and Skinder (2026)", "Thailand Statistician 24(1): 137-153", "zero-inflated Poisson-Garima"),
  c("HA21", "Hanandeh and Al-Nasser (2021)", "Austrian Journal of Statistics 50(4): 19-35", "Poisson-Gold")
), stringsAsFactors = FALSE)
names(papers) <- c("id", "authors", "where", "proposed")

## ---------- datasets (frequencies as printed; NA if not fully printed) ----------
D <- list(
  yeast187  = c(128, 37, 18, 3, 1),
  yeast400  = c(213, 128, 37, 18, 3, 1),
  mites80   = c(38, 17, 10, 9, 3, 2, 1),
  mites150  = c(70, 38, 17, 10, 9, 3, 2, 1),
  borer136  = c(83, 36, 14, 2, 1),
  borer324  = c(188, 83, 36, 14, 2, 1),
  chrom400  = c(268, 87, 26, 9, 4, 2, 1, 3),
  les60     = c(413, 124, 42, 15, 5, 0, 2),
  les70     = c(200, 57, 30, 7, 4, 0, 2),
  les90     = c(155, 83, 33, 14, 11, 3, 1),
  thJun     = c(187, 77, 40, 17, 6, 2, 1),
  thJul     = c(177, 80, 47, 26, 9, 2),
  thAug     = c(185, 89, 30, 24, 10, 3),
  thSum     = c(549, 246, 117, 67, 25, 7, 1),
  elbw60    = c(35, 11, 8, 4, 2),
  greece111 = c(39, 26, 17, 9, 6, 7, 6, 0, 0, 1),
  dic190    = c(92, 55, 27, 11, 4, 1),
  dic520    = c(437, 66, 15, 1, 1),
  dic631    = c(473, 119, 34, 3, 2)
)

## ---------- claims ----------
## type "ours": both winner and rival implemented here; dAIC recomputed.
##   cand: the models compared in the source (implemented ones).
## type "printed": dAIC from values printed in the source.
cl <- function(paper, data, winner, cand = NA, rival = NA, aic_w = NA, aic_r = NA, note = "")
  data.frame(paper = paper, data = data, winner = winner, cand = paste(cand, collapse = ","),
             rival = rival, aic_w = aic_w, aic_r = aic_r, note = note, stringsAsFactors = FALSE)
C3 <- c("POIS", "PL", "PS"); C3s <- c("POIS", "PL", "PSh")
claims <- rbind(
  ## SF16: Poisson, Poisson-Lindley, Poisson-Sujatha (claims by chi-square)
  cl("SF16", "yeast187", "POIS", C3, note = "Poisson stated as best"),
  cl("SF16", "mites80", "PL", C3), cl("SF16", "borer136", "PS", C3),
  cl("SF16", "chrom400", "PL", C3), cl("SF16", "les60", "PS", C3),
  cl("SF16", "les70", "PS", C3), cl("SF16", "les90", "PL", C3),
  ## SF17: Poisson, Poisson-Lindley, Poisson-Shanker
  cl("SF17", "yeast400", "POIS", C3s), cl("SF17", "mites80", "PSh", C3s),
  cl("SF17", "borer324", "PL", C3s), cl("SF17", "chrom400", "PSh", C3s),
  cl("SF17", "les60", "PSh", C3s), cl("SF17", "les70", "PSh", C3s),
  cl("SF17", "les90", "PSh", C3s), cl("SF17", "thJun", "PL", C3s),
  cl("SF17", "thJul", "PL", C3s), cl("SF17", "thAug", "PSh", C3s),
  cl("SF17", "thSum", "PL", C3s),
  ## TS25: Poisson, PL, P-Garima, P-Komal(?), P-Akash, P-Sujatha; -2logL printed
  cl("TS25", "borer324", "PS", rival = "PL (printed)", aic_w = 713.837 + 2, aic_r = 714.0992 + 2),
  cl("TS25", "thJun", "PS", rival = "PL (printed)", aic_w = 788.6261 + 2, aic_r = 788.8836 + 2),
  cl("TS25", "thJul", "PS", rival = "PL (printed)", aic_w = 879.6258 + 2, aic_r = 880.5501 + 2),
  ## SS19: Poisson-weighted Akash (two parameters) -- AIC printed
  cl("SS19", "mites150", "P-wAkash", rival = "P-wLindley (printed)", aic_w = 443.41, aic_r = 429.35,
     note = "claim by chi-square; printed AIC favours the rival"),
  cl("SS19", "thJun", "P-wAkash", rival = "PL (printed)", aic_w = 788.84 + 4, aic_r = 790.88,
     note = "printed AIC of the winner (788.73) inconsistent with its -2logL"),
  ## SS20: quasi Poisson-Aradhana -- chi-square only
  cl("SS20", "yeast400", "QPAD", note = "no likelihood printed"),
  cl("SS20", "borer324", "QPAD", note = "no likelihood printed"),
  cl("SS20", "mites80", "QPAD", note = "no likelihood printed"),
  ## RVP24: Poisson, PL, P-Akash, P-Suja
  cl("RVP24", "elbw60", "PSu", c("POIS", "PL", "PAk", "PSu")),
  ## MP25: two-parameter Poisson-Garima vs NB, discrete Weibull, ... (-loglik printed)
  cl("MP25", "thJun", "TPPGD", rival = "discrete Weibull (printed)", aic_w = 792.8389, aic_r = 793.1275),
  cl("MP25", "thJul", "TPPGD", rival = "discrete Weibull (printed)", aic_w = 883.9069, aic_r = 884.0968),
  cl("MP25", "greece111", "TPPGD", rival = "discrete Weibull (printed)", aic_w = 400.3046, aic_r = 400.5058),
  ## HE26: Poisson two-sum XLindley -- AIC printed
  cl("HE26", "dic190", "PTS-XL", rival = "COM-Poisson (printed)", aic_w = 483.001, aic_r = 484.303),
  ## ATBA24: authors conclude PL and PISDL are equally good
  cl("ATBA24", "dic520", "PL", rival = "PISDL (printed)", aic_w = 560.80, aic_r = 560.90,
     note = "authors call the two models equally competent"),
  cl("ATBA24", "dic631", "PISDL", rival = "PL (printed)", aic_w = 929.91, aic_r = 930.17,
     note = "authors call the two models equally competent"),
  ## MBO20: AICc printed (censored tables)
  cl("MBO20", "borers/hill (n=120)", "T2DS", rival = "NB (printed AICc)", aic_w = 403.65, aic_r = 404.71),
  cl("MBO20", "strikes (n=108)", "T1DS", rival = "T2DS (printed AICc)", aic_w = 569.42, aic_r = 570.55),
  ## DPA26: zero-one-inflated models, AIC printed
  cl("DPA26", "polio (n=168)", "ZOIPG", rival = "ZOIPL (printed)", aic_w = 536.0119, aic_r = 536.387),
  cl("DPA26", "syphilis (n=209)", "ZIPG", rival = "ZIPL (printed)", aic_w = 452.5824, aic_r = 453.8479),
  cl("DPA26", "drug offences (n=144)", "ZOIPG", rival = "ZOIPL (printed)", aic_w = 559.5282, aic_r = 561.2334),
  ## AS26: zero-inflated Poisson-Garima, AIC printed
  cl("AS26", "rabbit stillbirths (n=402)", "ZIPGD", rival = "ZINB (printed)", aic_w = 669.50, aic_r = 675.62),
  cl("AS26", "migrants (n=402)", "ZIPGD", rival = "ZINB (printed)", aic_w = 988.10, aic_r = 993.62),
  ## HA21: Poisson-Gold k = 1..4, Poisson-Lindley, Poisson-Shanker
  cl("HA21", "thJun", "PGD4", c("POIS", "PGD1", "PGD2", "PGD3", "PGD4", "PSh")),
  cl("HA21", "thSum", "PGD4", c("POIS", "PGD1", "PGD2", "PGD3", "PGD4", "PSh"))
)

pmf_fit <- function(G, fit) { p <- exp(fam_all[[G]]$lp(fit$par, xs)); p / sum(p) }
res <- do.call(rbind, lapply(seq_len(nrow(claims)), function(i) {
  c0 <- claims[i, ]
  f <- D[[c0$data]]; n <- if (!is.null(f)) sum(f) else as.numeric(sub(".*n=([0-9]+).*", "\\1", c0$data))
  mu <- if (!is.null(f)) sum((seq_along(f) - 1) * f) / n else NA
  out <- data.frame(paper = c0$paper, data = c0$data, n = n, mean = mu, winner = c0$winner,
                    rival = NA, dAIC = NA, n_star = NA, floor = NA, note = c0$note)
  if (c0$cand != "NA") {
    cand <- strsplit(c0$cand, ",")[[1]]
    fits <- fit_all(f, cand)
    aic <- sapply(cand, function(G) -2 * fits[[G]]$loglik + 2 * fam_all[[G]]$np)
    rv <- setdiff(cand, c0$winner); best_rival <- rv[which.min(aic[rv])]
    out$rival <- best_rival; out$dAIC <- aic[best_rival] - aic[c0$winner]
    fr <- frontier_pmf(pmf_fit(c0$winner, fits[[c0$winner]]), best_rival)
    out$n_star <- fr$n_star; out$floor <- fr$n_floor
  } else if (!is.na(c0$aic_w)) {
    out$rival <- c0$rival; out$dAIC <- c0$aic_r - c0$aic_w
  }
  out
}))
res$supported <- ifelse(is.na(res$dAIC), "no likelihood",
                  ifelse(res$dAIC < 0, "rival better", ifelse(res$dAIC < 2, "dAIC < 2", "dAIC >= 2")))
saveRDS(list(papers = papers, claims = res, data = D), "results/T_review.rds")
print(res[, c("paper", "data", "n", "winner", "rival", "dAIC", "n_star", "floor", "supported")], row.names = FALSE, digits = 4)
cat("\nClaims:", nrow(res), " papers:", length(unique(res$paper)), "\n")
print(table(res$supported))
k <- !is.na(res$floor)
cat("Implemented pairs:", sum(k), " with n below the floor:", sum(res$n[k] < res$floor[k]),
    " below n*:", sum(res$n[k] < res$n_star[k]), "\n")
