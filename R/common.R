# Shared settings and helpers of the simulation drivers.
#
# Every driver fits FLARE through the package (FLARE >= 3.0.0). The simulated
# data sets are not part of this repository; their location and the location
# of the comparison methods are given through environment variables.
suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(FLARE)
})
setDTthreads(1L)

REPO <- Sys.getenv("FLARE_SIMULATION", ".")                       # this repository
OUT <- Sys.getenv("FLARE_SIM_OUT", file.path(REPO, "output"))     # results and fits
DATA <- Sys.getenv("FLARE_SIM_DATA", file.path(REPO, "data"))     # simulated data sets

# Inputs produced by the data-generation steps.
PURE_TASKS <- file.path(DATA, "pure", "tasks.tsv")
SEMI_DATA <- file.path(DATA, "semi")
SEMI_REFERENCE <- file.path(DATA, "semi", "reference.rds")
NULL_GWAS <- file.path(DATA, "null_gwas")
FDR_GRID <- file.path(DATA, "false_discovery", "null_grid.tsv")
NOISE_REFERENCES <- file.path(DATA, "noise_reference")
# R libraries with the comparison methods (colon-separated), the library with
# GFA and its flashier version, and the FactorGo executable.
BENCHMARK_LIBRARY <- strsplit(Sys.getenv("FLARE_BENCHMARK_LIBRARY", ""), ":", fixed = TRUE)[[1]]
GFA_LIBRARY <- Sys.getenv("FLARE_GFA_LIBRARY", "")
FACTORGO <- Sys.getenv("FACTORGO", "factorgo")

out_path <- function(...) {
  p <- file.path(OUT, ...)
  dir.create(dirname(p), recursive = TRUE, showWarnings = FALSE)
  p
}

# Alignment and R2 scores of the benchmark (score_factorization,
# gleanr_exact_alignment, load_graph_pair, metric_meta, ...).
load_scoring <- function() {
  for (f in c("scoring.R", "factorwise_scores.R", "pure_helpers.R")) source(file.path(REPO, "R", f), local = globalenv())
}

# Sweep limits under which each family of simulations was run for the paper.
# They only matter for the few fits that reach the limit.
CONTROL <- list(
  pure = flare_control(max_cycles = c(1000, 500)),
  semi = flare_control(max_cycles = c(150, 200)),
  null = flare_control(max_cycles = 500)
)

# Noise reference stored for an N x M problem (20,000 draws, seed 1); it is
# simulated and cached on first use.
noise_reference <- function(N, M) {
  f <- file.path(NOISE_REFERENCES, sprintf("white_top_N%d_M%d_draws20000_seed1.rds", N, M))
  if (file.exists(f)) return(readRDS(f))
  cache <- out_path("reference", basename(f))
  if (!file.exists(cache)) saveRDS(flare_noise_reference(N, M), cache)
  readRDS(cache)
}

# The three FLARE variants reported in the figures, from two calls:
#   FLARE          both stages with the annotation graph
#   FLARE-NG       both stages with an empty graph (the matched comparison)
#   FLARE stage 1  the graph-free solution both of them start from
# `...` is passed to flare_fit() (C, C_se, sample_sd, K, control, ...).
flare_variants <- function(B, SE, graph, ..., no_graph = TRUE) {
  with_graph <- flare_fit(B, SE, graph = graph, ..., screen = FALSE, verbose = FALSE)
  out <- list(
    FLARE = list(U = with_graph$U, V = with_graph$V, fit = with_graph),
    `FLARE stage 1` = list(U = with_graph$stage1$U, V = with_graph$stage1$V, fit = with_graph)
  )
  if (no_graph) {
    without <- flare_fit(B, SE, graph = NULL, ..., screen = FALSE, verbose = FALSE)
    out$`FLARE-NG` <- list(U = without$U, V = without$V, fit = without)
  }
  out
}

trait_scale <- function(SE) apply(as.matrix(SE), 2, median)

# Scores of one factorization against the simulated truth.
score_all <- function(bundle, U, V, meta, Xhat = NULL, v_is_z = FALSE) {
  r <- append_formal20_scores(score_factorization(bundle, U, V, meta, Xhat_override = Xhat), bundle, U, V)
  b <- trait_scale(bundle$SE)
  Vz <- if (v_is_z) V else V / b
  r$V_R2_noise_units <- gleanr_exact_alignment(bundle$U_true, bundle$V_true / b, U, Vz)$V_r2
  r$zero_loading_traits <- sum(rowSums(abs(as.matrix(V)) > 1e-10) == 0)
  as.data.table(r)
}

# Per-trait recovery table of the semi-real simulation.
trait_table <- function(bundle, V, Xhat, group = NULL) {
  V <- as.matrix(V)
  V[!is.finite(V)] <- 0
  data.table(
    trait = colnames(bundle$B), group = if (is.null(group)) NA_character_ else group,
    b = trait_scale(bundle$SE), true_signal = colMeans((bundle$X_true / bundle$SE)^2),
    truly_loaded = rowSums(bundle$V_true != 0) > 0,
    loaded = if (ncol(V)) rowSums(abs(V) > 1e-10) > 0 else FALSE,
    recovery = vapply(seq_len(ncol(Xhat)), function(j) safe_cor(Xhat[, j], bundle$X_true[, j])^2, numeric(1)),
    fitted_energy = colMeans((Xhat / bundle$SE)^2)
  )
}

# Shrinkage of a sampling correlation without the block step, as used for the
# null GWAS design (nine traits, one block by construction).
shrink_only <- function(C, C_se) {
  off <- C
  diag(off) <- 0
  se <- C_se
  diag(se) <- 0
  g <- if (sum(off^2) > 0) min(1, sum(se^2) / sum(off^2)) else 0
  out <- g * diag(nrow(C)) + (1 - g) * C
  dimnames(out) <- dimnames(C)
  stopifnot(min(eigen(out, symmetric = TRUE, only.values = TRUE)$values) > 0)
  out
}
