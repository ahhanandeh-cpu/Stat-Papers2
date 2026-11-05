## =====================================================================
##  Benchmark count datasets (frequency f[i] = number of units with count i-1)
##
##  Sources (frequencies copied from the cited tables; each sample mean was
##  checked against the mean or Poisson estimate printed in the source):
##   [SF16] Shanker, R. and Fesshaye, H. (2016). On Poisson-Sujatha
##          distribution and its applications to model count data from
##          biological sciences. Biometrics & Biostatistics International
##          Journal 3(4): 00069. Tables 1-7.
##   [B79]  Bulmer, M.G. (1979). Principles of Statistics. Dover;
##          data of Greenwood, M. and Yule, G.U. (1920).
##   [HA21] Hanandeh, A.A. and Al-Nasser, A.D. (2021). Austrian Journal of
##          Statistics 50(4): 19-35, Tables 3-4; data of Falls, Williford
##          and Carter (1971).
##  Open-ended last cells ("4+" in corn borer, "7+" in chromatid aberrations)
##  are coded at their lower bound, which reproduces the sample means printed
##  in [SF16].
##
##  claimed: the model favoured in the source among those it compared
##  (by chi-square), as reported there.
## =====================================================================

benchmarks <- list(
  list(id = "yeast", name = "Yeast cells per haemocytometer square",
       source = "SF16 Table 1 (Student)", f = c(128, 37, 18, 3, 1),
       claimed = "Poisson (PS lowest $\\chi^2$)"),
  list(id = "mites", name = "European red mites on apple leaves",
       source = "SF16 Table 2", f = c(38, 17, 10, 9, 3, 2, 1),
       claimed = "PL"),
  list(id = "cornborer", name = "European corn borer larvae per plant",
       source = "SF16 Table 3 (McGuire et al.)", f = c(83, 36, 14, 2, 1),
       claimed = "PS"),
  list(id = "chromatid", name = "Chromatid aberrations (chinon 1, 24 h)",
       source = "SF16 Table 4 (Catcheside et al.)", f = c(268, 87, 26, 9, 4, 2, 1, 3),
       claimed = "PL (marginally)"),
  list(id = "lesions60", name = "Lymphoblast lesions, 60 ug/kg",
       source = "SF16 Table 5", f = c(413, 124, 42, 15, 5, 0, 2),
       claimed = "PS"),
  list(id = "lesions70", name = "Lymphoblast lesions, 70 ug/kg",
       source = "SF16 Table 6", f = c(200, 57, 30, 7, 4, 0, 2),
       claimed = "PS"),
  list(id = "lesions90", name = "Lymphoblast lesions, 90 ug/kg",
       source = "SF16 Table 7", f = c(155, 83, 33, 14, 11, 3, 1),
       claimed = "PL"),
  list(id = "accidents", name = "Accidents to women munitions workers",
       source = "B79 (Greenwood and Yule, 1920)", f = c(447, 132, 42, 21, 3, 2),
       claimed = "--"),
  list(id = "thunder_jun", name = "Thunderstorm days, June",
       source = "HA21 Table 3 (Falls et al., 1971)", f = c(187, 77, 40, 17, 6, 2, 1),
       claimed = "PA (PGD $k=4$)"),
  list(id = "thunder_sum", name = "Thunderstorm days, summer",
       source = "HA21 Table 4 (Falls et al., 1971)", f = c(549, 246, 117, 67, 25, 7, 1),
       claimed = "PGD $k=3,4$ (slightly)")
)
names(benchmarks) <- vapply(benchmarks, `[[`, "", "id")

## sanity checks against the sources
local({
  chk <- c(yeast = 0.459893, mites = 1.15, cornborer = 0.544118, chromatid = 0.5475,
           lesions60 = 0.47421, lesions70 = 0.55333, lesions90 = 0.853333)
  for (id in names(chk)) {
    f <- benchmarks[[id]]$f; m <- sum((seq_along(f) - 1) * f) / sum(f)
    stopifnot(abs(m - chk[[id]]) < 1e-4)
  }
  stopifnot(sum(benchmarks$accidents$f) == 647, sum(benchmarks$lesions60$f) == 601,
            sum(benchmarks$thunder_jun$f) == 330, sum(benchmarks$thunder_sum$f) == 1012)
})
