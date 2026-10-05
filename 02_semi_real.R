# Semi-real simulation (Figure 2c): factors estimated from All of Us, noise at
# five precision levels, two truths, 20 seeds = 200 data sets. All methods are
# run on each data set: FLARE (three variants), GLEANR, FactorGo, SVD, SVD-adj
# and SSVD. GFA and FLASH have their own drivers (05_gfa.R, 06_flash.R).
#
# Usage: Rscript 02_semi_real.R TASK        TASK = 1..200
source(file.path(Sys.getenv("FLARE_SIMULATION", "."), "R/common.R"))
design <- CJ(truth = c("A", "B"), rho = c(0.25, 0.5, 1, 2, 4), seed = 101:120, sorted = FALSE)
task <- as.integer(commandArgs(TRUE)[1])
stopifnot(task >= 1L, task <= nrow(design))
d <- design[task]
run_id <- sprintf("semi_v4_%s_rho_%03d_seed_%03d", d$truth, round(100 * d$rho), d$seed)
out <- out_path("results/semi", paste0(run_id, ".tsv"))
if (file.exists(out)) quit(status = 0)

load_scoring()
reference <- readRDS(SEMI_REFERENCE)
.libPaths(unique(c(BENCHMARK_LIBRARY, .libPaths())))
source(file.path(REPO, "R/benchmark_methods.R"))
Sys.setenv(JAX_PLATFORMS = "cpu")
bundle <- readRDS(file.path(SEMI_DATA, d$truth, sprintf("rho_%03d", round(100 * d$rho)), sprintf("seed_%03d.rds", d$seed)))
stopifnot(identical(rownames(reference$A_true), rownames(bundle$B)))
K <- 6L

# GLEANR with its own grid search, as in the original benchmark.
gleanr_fit <- function() {
  X <- as.matrix(bundle$B)
  W <- 1 / as.matrix(bundle$SE)
  work <- out_path("work/semi", run_id, "gleanr")
  g <- gleanr::initializeGLEANR(X = X, W = W, C = bundle$C, snp.ids = rownames(X), trait.names = colnames(X), K = K,
    init.mat = "V", covar_shrinkage = -1, is.sim = TRUE, save_out = TRUE, min_bic_search_iter = 5L, conv_objective = 0.001, verbosity = 0)
  o <- g$options
  o$scale <- FALSE; o$svd_init <- TRUE; o$bic.var <- "sklearn_eBIC"; o$regression_method <- "glmnet"
  o$iter <- 200L; o$save_out <- TRUE; o$ncores <- 1L; o$nsplits <- 1L; o$parallel <- FALSE
  o$Kmin <- K; o$K <- "GRID"; o$conv0 <- 0.001; o$min.bicsearch.iter <- 5L
  reg <- gleanr::prepRegressionElements(X, W, g$W_c, o)
  bic <- gleanr::getBICMatricesGLMNET(work, o, X, W, g$W_c, g$all_ids, g$names, min.iter = 5L, max.iter = 200L, reg.elements = reg)
  U <- matrix(0, nrow(X), K)
  V <- matrix(0, ncol(X), K)
  if (!is.null(bic$K) && is.finite(bic$K) && !is.na(bic$alpha) && !all(bic$optimal.v == 0)) {
    o <- bic$options
    o$alpha1 <- bic$alpha; o$lambda1 <- bic$lambda; o$K <- bic$K; o$iter <- 200L
    o$ncores <- 1L; o$nsplits <- 1L; o$parallel <- FALSE
    f <- gleanr::gwasML_ALS_Routine(o, X, W, g$W_c, bic$optimal.v, maxK = K, reg.elements = reg)
    k <- min(ncol(f$U), K)
    if (k > 0) {
      U[, seq_len(k)] <- as.matrix(f$U)[, seq_len(k)]
      V[, seq_len(k)] <- as.matrix(f$V)[, seq_len(k)]
    }
  }
  dimnames(U) <- list(rownames(X), NULL)
  dimnames(V) <- list(colnames(X), NULL)
  list(U = U, V = V, Xhat = U %*% t(V))
}
flare_all <- function() {
  fits <- flare_variants(bundle$B, bundle$SE, reference$A_true, C = bundle$C, K = K, control = CONTROL$semi)
  lapply(fits, function(x) list(U = x$U, V = x$V, Xhat = x$U %*% t(x$V)))
}
bench <- function(m) {
  setNames(list(fit_formal20_benchmark(m, bundle, K, sample_sizes = bundle$sample_sizes,
    factorgo_executable = FACTORGO, work_dir = out_path("work/semi", run_id, "factorgo"))), m)
}
methods <- list(FLARE = flare_all, GLEANR = function() list(GLEANR = gleanr_fit()),
                FactorGo = function() bench("FactorGo"),
                SVD = function() bench("SVD"), `SVD-adj` = function() bench("SVD-adj"), SSVD = function() bench("SSVD"))
only <- Sys.getenv("FLARE_TEST_METHODS", "")
if (nzchar(only)) methods <- methods[strsplit(only, ",")[[1]]]
z_scale <- c("SVD", "SVD-adj", "SSVD", "FactorGo")   # methods that return V on the Z scale

rows <- list()
traits <- list()
for (m in names(methods)) {
  started <- Sys.time()
  res <- tryCatch(methods[[m]](), error = function(e) conditionMessage(e))
  seconds <- as.numeric(difftime(Sys.time(), started, units = "secs"))
  id <- data.table(run_id = run_id, truth = d$truth, rho = d$rho, seed = d$seed,
                   replicate_id = d$seed - 100L, signal_scale = bundle$meta$signal_scale)
  if (is.character(res)) {
    rows[[length(rows) + 1L]] <- cbind(id, method_label = m, status = paste("error:", res), seconds = seconds)
    next
  }
  for (label in names(res)) {
    x <- res[[label]]
    V <- as.matrix(x$V)
    if (is.null(rownames(V))) rownames(V) <- colnames(bundle$B)
    rows[[length(rows) + 1L]] <- cbind(id, method_label = label, status = "ok", seconds = seconds,
      score_all(bundle, x$U, V, list(), Xhat = x$Xhat, v_is_z = label %in% z_scale))
    traits[[length(traits) + 1L]] <- cbind(id, method_label = label, trait_table(bundle, V, x$Xhat, bundle$group))
    saveRDS(list(U = x$U, V = V), out_path("fits/semi", sprintf("%s__%s.rds", run_id, gsub("[^A-Za-z0-9]+", "_", label))))
  }
}
fwrite(rbindlist(traits, fill = TRUE), out_path("results/semi_traits", paste0(run_id, ".tsv")), sep = "\t")
fwrite(rbindlist(rows, fill = TRUE), out, sep = "\t")
cat("SEMI_DONE", run_id, "\n")
