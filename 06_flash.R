# FLASH on the pure, semi-real and null-GWAS data sets, with the noise variance
# treated as known (var_type = NULL). With var_type = 0 flashier would fit one
# extra variance, shared by every entry of the matrix, on top of the standard
# errors.
#
# Two calls with known noise are run on every data set:
#   FLASH-Z    Z-scores with unit noise variance
#   FLASH-SE   effect estimates with their standard errors as the noise
# flashier 0.2.51, greedy fit capped at the true number of factors, backfit and
# null check. In the null design no cap is given and only the number of
# factors is recorded.
#
# Usage: Rscript 06_flash.R pure TASK      TASK = 1..60, five data sets each
#        Rscript 06_flash.R semi TASK      TASK = 1..200, one data set each
#        Rscript 06_flash.R null TASK      TASK = 1..20, 42 data sets each
source(file.path(Sys.getenv("FLARE_SIMULATION", "."), "R/common.R"))
.libPaths(unique(c(BENCHMARK_LIBRARY, .libPaths())))
suppressPackageStartupMessages(library(flashier))
args <- commandArgs(TRUE)
family <- match.arg(args[1], c("pure", "semi", "null"))
task <- as.integer(args[2])

CALLS <- c(`FLASH-Z` = "z", `FLASH-SE` = "se")

flash_known_noise <- function(B, SE, input, K = NULL) {
  a <- if (input == "z") list(data = B / SE, S = 1) else list(data = B, S = SE)
  a <- c(a, list(var_type = NULL, backfit = TRUE, nullcheck = TRUE, verbose = 0L))
  if (!is.null(K)) a$greedy_Kmax <- K
  do.call(flashier::flash, a)
}

# Same scale convention as the other comparison methods: each factor carries
# equal norm in U and V, which leaves U V' unchanged.
balance <- function(U, V) {
  for (k in seq_len(ncol(U))) {
    un <- sqrt(sum(U[, k]^2))
    vn <- sqrt(sum(V[, k]^2))
    if (un < 1e-14 || vn < 1e-14) {
      U[, k] <- 0
      V[, k] <- 0
    } else {
      U[, k] <- U[, k] * sqrt(vn / un)
      V[, k] <- V[, k] * sqrt(un / vn)
    }
  }
  list(U = U, V = V)
}

# Fit one call, bring the result to exactly K columns (strongest first,
# zero-padded) and score it against the truth.
run_one <- function(bundle, K, id, label) {
  input <- CALLS[[label]]
  started <- Sys.time()
  f <- tryCatch(flash_known_noise(bundle$B, bundle$SE, input, K), error = function(e) conditionMessage(e))
  seconds <- as.numeric(difftime(Sys.time(), started, units = "secs"))
  if (is.character(f)) {
    return(list(row = cbind(id, method_label = label, status = paste("error:", f), seconds = seconds), fit = NULL))
  }
  N <- nrow(bundle$B)
  M <- ncol(bundle$B)
  K_found <- as.integer(f$n_factors)
  U <- matrix(0, N, K)
  V <- matrix(0, M, K)
  if (K_found > 0) {
    b <- balance(as.matrix(f$L_pm), as.matrix(f$F_pm))
    strongest <- order(-colSums(b$U^2) * colSums(b$V^2))[seq_len(min(K, K_found))]
    U[, seq_along(strongest)] <- b$U[, strongest, drop = FALSE]
    V[, seq_along(strongest)] <- b$V[, strongest, drop = FALSE]
  }
  dimnames(U) <- list(rownames(bundle$B), NULL)
  dimnames(V) <- list(colnames(bundle$B), NULL)
  # Fitted effects: Z-scale factors are multiplied back by the standard errors.
  Xhat <- U %*% t(V)
  if (input == "z") Xhat <- Xhat * bundle$SE
  scores <- score_all(bundle, U, V, list(), Xhat = Xhat, v_is_z = input == "z")
  list(row = cbind(id, method_label = label, status = "ok", seconds = seconds, K_found = K_found, scores),
       fit = list(U = U, V = V))
}

run_both <- function(bundle, K, id, run_id, fit_dir) {
  rows <- list()
  for (label in names(CALLS)) {
    r <- run_one(bundle, K, id, label)
    if (!is.null(r$fit)) saveRDS(r$fit, out_path(fit_dir, sprintf("%s__%s.rds", run_id, gsub("[^A-Za-z0-9]+", "_", label))))
    rows[[label]] <- r$row
  }
  rbindlist(rows, fill = TRUE)
}

if (family == "pure") {
  stopifnot(task >= 1L, task <= 60L)
  load_scoring()
  tasks <- fread(PURE_TASKS)[family == "pure"]
  stopifnot(nrow(tasks) == 300L)
  limit <- as.integer(Sys.getenv("FLARE_TEST_LIMIT", "5"))
  for (i in ((task - 1L) * 5L + 1:5)[seq_len(limit)]) {
    job <- tasks[i]
    out <- out_path("results/pure_flash", paste0(job$run_id, ".tsv"))
    if (file.exists(out)) next
    bundle <- readRDS(job$bundle_path)
    id <- data.table(run_id = job$run_id, architecture = job$architecture, condition = job$condition)
    fwrite(run_both(bundle, ncol(bundle$U_true), id, job$run_id, "fits/pure_flash"), out, sep = "\t")
    cat("FLASH_PURE_DONE", job$run_id, "\n")
  }
} else if (family == "semi") {
  load_scoring()
  design <- CJ(truth = c("A", "B"), rho = c(0.25, 0.5, 1, 2, 4), seed = 101:120, sorted = FALSE)
  stopifnot(task >= 1L, task <= nrow(design))
  d <- design[task]
  run_id <- sprintf("semi_v4_%s_rho_%03d_seed_%03d", d$truth, round(100 * d$rho), d$seed)
  out <- out_path("results/semi_flash", paste0(run_id, ".tsv"))
  if (!file.exists(out)) {
    bundle <- readRDS(file.path(SEMI_DATA, d$truth, sprintf("rho_%03d", round(100 * d$rho)), sprintf("seed_%03d.rds", d$seed)))
    id <- data.table(run_id = run_id, truth = d$truth, rho = d$rho, seed = d$seed,
                     replicate_id = d$seed - 100L, signal_scale = bundle$meta$signal_scale)
    fwrite(run_both(bundle, 6L, id, run_id, "fits/semi_flash"), out, sep = "\t")
    cat("FLASH_SEMI_DONE", run_id, "\n")
  }
} else {
  stopifnot(task >= 1L, task <= 20L)
  sims <- file.path(NULL_GWAS, "sims")
  stems <- sort(sub("\\.BETA\\.csv$", "", list.files(sims, pattern = "\\.BETA\\.csv$")))
  stopifnot(length(stems) == 840L)
  limit <- as.integer(Sys.getenv("FLARE_TEST_LIMIT", "42"))
  count <- function(B, SE, input) {
    tryCatch(as.integer(flash_known_noise(B, SE, input)$n_factors), error = function(e) NA_integer_)
  }
  for (stem in stems[((task - 1L) * 42L + 1:42)[seq_len(limit)]]) {
    out <- out_path("results/null_gwas_flash", paste0(stem, ".tsv"))
    if (file.exists(out)) next
    B <- as.matrix(fread(file.path(sims, paste0(stem, ".BETA.csv"))))
    SE <- as.matrix(fread(file.path(sims, paste0(stem, ".SE.csv"))))
    fwrite(data.table(case = stem, `FLASH-Z` = count(B, SE, "z"), `FLASH-SE` = count(B, SE, "se")), out, sep = "\t")
    cat("FLASH_NULL_DONE", stem, "\n")
  }
}
