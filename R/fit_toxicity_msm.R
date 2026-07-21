suppressMessages(library(msm))

INTERIM <- file.path("data", "interim")
OUTPUT  <- file.path("data", "output")
dir.create(OUTPUT, showWarnings = FALSE, recursive = TRUE)

d <- read.csv(file.path(INTERIM, "toxicity_panel.csv"))
d$da_cat <- as.integer(d$da_cat)
d <- d[order(d$subj, d$t_years), ]
cat(sprintf("data: %d DA obs, %d sites; category counts: %s\n",
            nrow(d), length(unique(d$subj)),
            paste(table(d$da_cat), collapse = "/")))

fit_K <- function(K) {
  Q <- matrix(0, K, K)
  for (i in 1:K) { if (i > 1) Q[i, i-1] <- 0.3; if (i < K) Q[i, i+1] <- 0.3 }
  hmodel <- lapply(1:K, function(k) {
    frac <- if (K == 1) 0 else (k - 1) / (K - 1)
    pr <- c((1 - frac) * 0.88 + 0.04, 0.10 + 0.10 * frac, frac * 0.6 + 0.02)
    msm::hmmCat(prob = pr / sum(pr))
  })
  try(msm(da_cat ~ t_years, subject = subj, data = d, qmatrix = Q,
          hmodel = hmodel, control = list(fnscale = 4000, maxit = 20000)),
      silent = TRUE)
}

rows <- list()
for (K in 2:4) {
  fit <- fit_K(K)
  if (inherits(fit, "try-error")) { cat(sprintf("K=%d FAILED\n", K)); next }
  ll <- as.numeric(logLik.msm(fit)); np <- attr(logLik.msm(fit), "df"); n <- nrow(d)
  cat(sprintf("\n== K=%d ==  logLik=%.2f npar=%d AIC=%.1f BIC=%.1f conv=%d\n",
              K, ll, np, -2*ll + 2*np, -2*ll + log(n)*np, fit$opt$convergence))
  print(round(qmatrix.msm(fit, ci = "none"), 3))
  cat("mean sojourn (yr):", round(sojourn.msm(fit)$estimates, 3), "\n")
  rows[[as.character(K)]] <- data.frame(
    K = K, logLik = round(ll, 1), npar = np,
    AIC = round(-2*ll + 2*np, 1), BIC = round(-2*ll + log(n)*np, 1))
  saveRDS(fit, file.path(OUTPUT, sprintf("msm_tox_K%d.rds", K)))
}

tab <- do.call(rbind, rows)
cat("\n==== model selection ====\n"); print(tab, row.names = FALSE)
cat(sprintf("Best AIC: K=%d   Best BIC: K=%d\n",
            tab$K[which.min(tab$AIC)], tab$K[which.min(tab$BIC)]))
write.csv(tab, file.path(OUTPUT, "msm_tox_modelsel.csv"), row.names = FALSE)
