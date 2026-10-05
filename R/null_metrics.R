# Shared helpers for the false-discovery simulations. No truth-based alignment is
# needed under the null, so every quantity below is defined for any method output.

ZERO_TOL <- 1e-10

# Factors that are actually present: both the SNP column and the trait column are non-zero.
active_factors <- function(U, V, tol = ZERO_TOL) {
  if (!length(U) || !length(V) || !ncol(U)) return(integer())
  which(colSums(abs(U) > tol) > 0 & colSums(abs(V) > tol) > 0)
}

# Share of the standardized data sum of squares reproduced by the fitted effects.
# Under the global null all of it is noise, so this is spurious variance explained.
standardized_energy <- function(B, SE, Xhat) {
  sum((Xhat / SE)^2) / sum((B / SE)^2)
}

null_metrics <- function(B, SE, U, V, Xhat = NULL, null_traits = character(),
                         V_true = NULL, tol = ZERO_TOL) {
  U <- as.matrix(U); V <- as.matrix(V)
  U[!is.finite(U)] <- 0; V[!is.finite(V)] <- 0
  if (is.null(Xhat)) Xhat <- if (ncol(U)) U %*% t(V) else matrix(0, nrow(B), ncol(B))
  act <- active_factors(U, V, tol)
  nzV <- abs(V) > tol
  # A factor makes a cross-trait claim only when it is multi-trait: after scaling the
  # loading vector to unit length, no single trait exceeds 0.98 in absolute value.
  # This is the rule used for the real-data analyses.
  top <- if (length(act)) vapply(act, function(k) max(abs(V[, k])) / sqrt(sum(V[, k]^2)), numeric(1L)) else numeric()
  out <- list(
    K_returned = ncol(U),
    K_active = length(act),
    K_multi_trait = sum(top <= 0.98),
    K_single_trait = sum(top > 0.98),
    V_nonzero_entries = sum(nzV),
    V_nonzero_fraction = if (length(nzV)) mean(nzV) else 0,
    traits_with_loading = sum(rowSums(nzV) > 0),
    U_nonzero_fraction = if (length(U)) mean(abs(U) > tol) else 0,
    spurious_energy = standardized_energy(B, SE, Xhat)
  )
  if (length(null_traits)) {
    idx <- match(null_traits, rownames(V))
    stopifnot(!anyNA(idx))
    energy <- sum(V^2)
    out$null_traits <- length(idx)
    out$null_traits_loaded <- sum(rowSums(nzV[idx, , drop = FALSE]) > 0)
    out$null_trait_false_positive_rate <- out$null_traits_loaded / length(idx)
    out$loading_energy_on_null_traits <- if (energy > 0) sum(V[idx, ]^2) / energy else 0
    # Relative-magnitude version for dense methods, which never return exact zeros:
    # a null trait counts as loaded when its largest normalized loading exceeds 0.05.
    Vn <- sweep(V, 2L, pmax(sqrt(colSums(V^2)), 1e-300), "/")
    out$null_traits_loaded_rel005 <- sum(apply(abs(Vn[idx, , drop = FALSE]), 1L, max) > 0.05)
  }
  out
}

# Entry-level false discovery proportion and power against a known support.
# `est` must already be aligned to `truth` (same column order).
support_fdp <- function(truth, est, tol = ZERO_TOL, relative = NULL) {
  truth <- as.matrix(truth); est <- as.matrix(est)
  stopifnot(identical(dim(truth), dim(est)))
  S <- truth != 0
  if (is.null(relative)) {
    P <- abs(est) > tol
  } else {
    scale <- pmax(apply(abs(est), 2L, max), 1e-300)
    P <- sweep(abs(est), 2L, scale, "/") > relative
  }
  discoveries <- sum(P)
  c(discoveries = discoveries,
    false_discoveries = sum(P & !S),
    FDP = if (discoveries) sum(P & !S) / discoveries else 0,
    power = if (sum(S)) sum(P & S) / sum(S) else NA_real_,
    truth_nonzero_fraction = mean(S))
}
