# GFA (Morrison et al.) on the pure and semi-real simulation data sets, scored
# with the same functions as every other method.
#
# GFA is fitted to Z-scores with the per-trait sample sizes and the same
# sampling correlation the other methods receive. It selects its own number of
# factors; the true number is given as the maximum, as for FLASH. Factors are
# returned on the Z scale and handled like the other Z-scale methods (SVD,
# FactorGo): the fitted effects are U V' multiplied by the standard errors.
#
# GFA needs flashier >= 1.0, which is newer than the version used for the FLASH
# comparison, so it has its own library and its own driver.
#
# Usage: Rscript 05_gfa.R pure TASK      TASK = 1..60, five data sets each
#        Rscript 05_gfa.R semi TASK      TASK = 1..200, one data set each
source(file.path(Sys.getenv("FLARE_SIMULATION", "."), "R/common.R"))
if (nzchar(GFA_LIBRARY)) .libPaths(c(GFA_LIBRARY, .libPaths()))
suppressPackageStartupMessages(library(GFA))
args <- commandArgs(TRUE)
family <- match.arg(args[1], c("pure", "semi"))
task <- as.integer(args[2])
load_scoring()

# One GFA fit. Returns NULL when GFA finds no factor.
#
# gfa_fit() of GFA 1.0.0 fails in its duplicate check when the fit holds no
# free factor, so the fit is run without the wrap-up, free factors are
# counted, and the wrap-up steps are applied only when there is at least one.
fit_gfa <- function(bundle, K) {
  Z <- bundle$B / bundle$SE
  M <- ncol(Z)
  R <- bundle$C
  if (max(abs(R - diag(M))) < 1e-12) {
    R <- NULL
  } else {
    # GFA refuses an ill-conditioned R; project it as its error message asks.
    ev <- eigen(R, symmetric = TRUE, only.values = TRUE)$values
    if (min(ev) <= 0 || max(ev) / min(ev) > 1000) {
      R <- as.matrix(Matrix::nearPD(R, corr = TRUE, posd.tol = 1 / 1000)$mat)
      dimnames(R) <- dimnames(bundle$C)
    }
  }
  params <- gfa_default_parameters()
  params$kmax <- K
  set.seed(23)
  g <- gfa_fit(Z_hat = Z, N = as.numeric(bundle$sample_sizes), R = R, params = params, no_wrapup = TRUE)
  fit <- g$fit
  method <- fit$method     # flash_nullcheck() returns an object without this field
  n_free <- function(f) {
    fixed <- f$flash_fit$fix.dim
    f$n_factors - if (length(fixed)) sum(!vapply(fixed, is.null, logical(1))) else 0L
  }
  if (n_free(fit) == 0) return(NULL)
  fit <- flashier::flash_nullcheck(fit, remove = TRUE)
  if (n_free(fit) == 0) return(NULL)
  fit <- GFA:::gfa_duplicate_check(fit, dim = 2, check_thresh = g$params$duplicate_check_thresh)
  out <- GFA:::gfa_wrapup(fit, method = method, scale = g$scale, nullcheck = TRUE)
  # The wrap-up divides the trait loadings by sqrt(N); put them back on the Z scale.
  list(U = as.matrix(out$L_hat), V = as.matrix(out$F_hat) * g$scale,
       iteration_limit = !is.null(fit$flash_fit$maxiter.reached))
}

# Fit, bring to exactly K columns (strongest first, zero-padded) and score.
run_one <- function(bundle, K, id) {
  started <- Sys.time()
  res <- tryCatch(suppressMessages(suppressWarnings(fit_gfa(bundle, K))), error = function(e) conditionMessage(e))
  seconds <- as.numeric(difftime(Sys.time(), started, units = "secs"))
  if (is.character(res)) {
    return(list(row = cbind(id, method_label = "GFA", status = paste("error:", res), seconds = seconds), fit = NULL))
  }
  N <- nrow(bundle$B)
  M <- ncol(bundle$B)
  K_found <- if (is.null(res)) 0L else ncol(res$U)
  U <- matrix(0, N, K)
  V <- matrix(0, M, K)
  if (K_found > 0) {
    strongest <- order(-colSums(res$U^2) * colSums(res$V^2))[seq_len(min(K, K_found))]
    U[, seq_along(strongest)] <- res$U[, strongest, drop = FALSE]
    V[, seq_along(strongest)] <- res$V[, strongest, drop = FALSE]
  }
  dimnames(U) <- list(rownames(bundle$B), NULL)
  dimnames(V) <- list(colnames(bundle$B), NULL)
  Xhat <- (U %*% t(V)) * bundle$SE
  scores <- score_all(bundle, U, V, list(), Xhat = Xhat, v_is_z = TRUE)
  list(row = cbind(id, method_label = "GFA", status = "ok", seconds = seconds, K_found = K_found,
                   iteration_limit = isTRUE(res$iteration_limit), scores),
       fit = list(U = U, V = V))
}

if (family == "pure") {
  stopifnot(task >= 1L, task <= 60L)
  tasks <- fread(PURE_TASKS)[family == "pure"]
  stopifnot(nrow(tasks) == 300L)
  limit <- as.integer(Sys.getenv("FLARE_TEST_LIMIT", "5"))
  for (i in ((task - 1L) * 5L + 1:5)[seq_len(limit)]) {
    job <- tasks[i]
    out <- out_path("results/pure_gfa", paste0(job$run_id, ".tsv"))
    if (file.exists(out)) next
    bundle <- readRDS(job$bundle_path)
    id <- data.table(run_id = job$run_id, architecture = job$architecture, condition = job$condition)
    r <- run_one(bundle, ncol(bundle$U_true), id)
    if (!is.null(r$fit)) saveRDS(r$fit, out_path("fits/pure_gfa", paste0(job$run_id, "__GFA.rds")))
    fwrite(r$row, out, sep = "\t")
    cat("GFA_PURE_DONE", job$run_id, r$row$status, "\n")
  }
} else {
  design <- CJ(truth = c("A", "B"), rho = c(0.25, 0.5, 1, 2, 4), seed = 101:120, sorted = FALSE)
  stopifnot(task >= 1L, task <= nrow(design))
  d <- design[task]
  run_id <- sprintf("semi_v4_%s_rho_%03d_seed_%03d", d$truth, round(100 * d$rho), d$seed)
  out <- out_path("results/semi_gfa", paste0(run_id, ".tsv"))
  if (!file.exists(out)) {
    bundle <- readRDS(file.path(SEMI_DATA, d$truth, sprintf("rho_%03d", round(100 * d$rho)), sprintf("seed_%03d.rds", d$seed)))
    id <- data.table(run_id = run_id, truth = d$truth, rho = d$rho, seed = d$seed,
                     replicate_id = d$seed - 100L, signal_scale = bundle$meta$signal_scale)
    r <- run_one(bundle, 6L, id)
    if (!is.null(r$fit)) saveRDS(r$fit, out_path("fits/semi_gfa", paste0(run_id, "__GFA.rds")))
    fwrite(r$row, out, sep = "\t")
    cat("GFA_SEMI_DONE", run_id, r$row$status, "\n")
  }
}
