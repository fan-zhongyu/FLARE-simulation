# Adapters for comparison methods run through one interface: truncated SVD, SVD after whitening,
# sparse SVD and FactorGo.
# All returned fits use a canonical, product-preserving equal-L2 scale for U/V.

balance_factor_pairs <- function(U, V, tol = 1e-14) {
  U <- as.matrix(U); V <- as.matrix(V)
  if (ncol(U) != ncol(V)) stop("U/V ranks differ.")
  for (k in seq_len(ncol(U))) {
    un <- sqrt(sum(U[, k]^2)); vn <- sqrt(sum(V[, k]^2))
    if (!is.finite(un) || !is.finite(vn) || un <= tol || vn <= tol) {
      U[, k] <- 0; V[, k] <- 0
    } else {
      scale <- sqrt(vn / un)
      U[, k] <- U[, k] * scale
      V[, k] <- V[, k] / scale
    }
  }
  list(U = U, V = V)
}

pad_factor_pairs <- function(U, V, K) {
  U <- as.matrix(U); V <- as.matrix(V)
  if (ncol(U) != ncol(V)) stop("U/V ranks differ.")
  if (ncol(U) > K) stop("Estimated rank exceeds K cap.")
  K_active <- ncol(U)
  if (K_active < K) {
    U <- cbind(U, matrix(0, nrow(U), K - K_active))
    V <- cbind(V, matrix(0, nrow(V), K - K_active))
  }
  list(U = U[, seq_len(K), drop = FALSE],
       V = V[, seq_len(K), drop = FALSE], K_active = K_active)
}

name_factor_pairs <- function(fit, bundle) {
  K <- ncol(fit$U)
  dimnames(fit$U) <- list(rownames(bundle$B), paste0("F", seq_len(K)))
  dimnames(fit$V) <- list(colnames(bundle$B), paste0("F", seq_len(K)))
  fit
}

# Raw-coordinate truncated SVD.  Unlike fit_svd_z(), the returned U and V
# factorize the effect matrix B itself, so U/V recovery is directly comparable
# with the raw-coordinate simulation truth.
fit_svd_raw <- function(bundle, K) {
  s <- svd(bundle$B, nu = K, nv = K)
  d <- s$d[seq_len(K)]
  U <- sweep(s$u[, seq_len(K), drop = FALSE], 2L, sqrt(d), "*")
  V <- sweep(s$v[, seq_len(K), drop = FALSE], 2L, sqrt(d), "*")
  uv <- balance_factor_pairs(U, V)
  name_factor_pairs(list(
    U = uv$U, V = uv$V, Xhat = uv$U %*% t(uv$V),
    K_active = K, input_coordinate = "raw_B",
    reconstruction_coordinate = "raw_B=U V^T",
    rank_contract = "fixed_K", singular_values = d
  ), bundle)
}

fit_svd_z <- function(bundle, K) {
  Z <- bundle$B / bundle$SE
  s <- svd(Z, nu = K, nv = K)
  d <- s$d[seq_len(K)]
  U <- sweep(s$u[, seq_len(K), drop = FALSE], 2L, sqrt(d), "*")
  V <- sweep(s$v[, seq_len(K), drop = FALSE], 2L, sqrt(d), "*")
  uv <- balance_factor_pairs(U, V)
  Xhat_z <- uv$U %*% t(uv$V)
  name_factor_pairs(list(
    U = uv$U, V = uv$V, Xhat = Xhat_z * bundle$SE,
    K_active = K, input_coordinate = "Z=B/SE",
    reconstruction_coordinate = "raw_B=(U_Z V_Z^T)*SE",
    rank_contract = "fixed_K", singular_values = d
  ), bundle)
}

blockify_covariance <- function(C, threshold = 0.2) {
  C <- as.matrix(C)
  if (!isSymmetric(C, tol = 1e-10)) C <- (C + t(C)) / 2
  ids <- colnames(C)
  if (is.null(ids)) ids <- as.character(seq_len(ncol(C)))
  rownames(C) <- colnames(C) <- ids
  hc <- stats::hclust(stats::as.dist(1 - abs(C)))
  group <- stats::cutree(hc, h = 1 - threshold)
  keep <- as.integer(names(table(group))[table(group) >= 2L])
  Cb <- diag(diag(C), nrow(C), ncol(C))
  dimnames(Cb) <- dimnames(C)
  for (g in keep) {
    idx <- which(group == g)
    Cb[idx, idx] <- C[idx, idx]
  }
  Cb
}

fit_svd_adj_z <- function(bundle, K, block_threshold = 0.2) {
  Z <- bundle$B / bundle$SE
  Cb <- blockify_covariance(bundle$C, block_threshold)
  Wc <- solve(chol(Cb))
  Zw <- t(Wc %*% t(Z))
  s <- svd(Zw, nu = K, nv = K)
  d <- s$d[seq_len(K)]
  U <- sweep(s$u[, seq_len(K), drop = FALSE], 2L, sqrt(d), "*")
  Vw <- sweep(s$v[, seq_len(K), drop = FALSE], 2L, sqrt(d), "*")
  # Zw = Z %*% t(Wc), so V must be mapped back to the original Z traits.
  V <- solve(Wc) %*% Vw
  uv <- balance_factor_pairs(U, V)
  Xhat_z <- uv$U %*% t(uv$V)
  name_factor_pairs(list(
    U = uv$U, V = uv$V, Xhat = Xhat_z * bundle$SE,
    K_active = K, input_coordinate = "whitened Z=(B/SE)%*%t(Wc)",
    reconstruction_coordinate = "unwhitened raw_B=(U_Z V_Z^T)*SE",
    rank_contract = "fixed_K", singular_values = d,
    covariance_block_threshold = block_threshold,
    covariance_block = Cb,
    maps_V_back_to_original_trait_coordinate = TRUE
  ), bundle)
}

# Exact reproduction of the SVD_whiten comparator in the public GLEANR
# manuscript workflow.  That workflow factorizes the covariance-whitened
# Z-score matrix and reports its singular vectors directly, without mapping
# the trait singular vectors back through the inverse whitener.  Keep this
# separate from fit_svd_adj_z(): the latter is an algebraically unwhitened
# extension, whereas this function is the publication-faithful comparator.
fit_svd_adj_gleanr_official <- function(bundle, K, block_threshold = 0.2) {
  Z <- bundle$B / bundle$SE
  Cb <- blockify_covariance(bundle$C, block_threshold)
  Wc <- solve(chol(Cb))
  Zw <- t(Wc %*% t(Z))
  s <- svd(Zw, nu = K, nv = K)
  d <- s$d[seq_len(K)]
  U <- sweep(s$u[, seq_len(K), drop = FALSE], 2L, sqrt(d), "*")
  V <- sweep(s$v[, seq_len(K), drop = FALSE], 2L, sqrt(d), "*")
  uv <- balance_factor_pairs(U, V)
  Xhat_zw <- uv$U %*% t(uv$V)
  name_factor_pairs(list(
    U = uv$U, V = uv$V, Xhat = Xhat_zw * bundle$SE,
    K_active = K,
    input_coordinate = "GLEANR-official whitened Z=(B/SE)%*%t(Wc)",
    reconstruction_coordinate =
      "GLEANR-official comparator: whitened-Z reconstruction multiplied by SE",
    rank_contract = "fixed_K",
    singular_values = d,
    covariance_block_threshold = block_threshold,
    covariance_block = Cb,
    maps_V_back_to_original_trait_coordinate = FALSE
  ), bundle)
}

fit_yang_ssvd_z <- function(bundle, K) {
  if (!requireNamespace("ssvd", quietly = TRUE)) stop("ssvd is not installed.")
  # ssvd::ssvd uses random starts. Tie them to the simulated dataset so the
  # formal benchmark is reproducible across processes and reruns.
  set.seed(as.integer(bundle$meta$seed) + 731L)
  Z <- bundle$B / bundle$SE
  warnings_seen <- character()
  s <- withCallingHandlers(
    ssvd::ssvd(Z, method = "method", r = K),
    warning = function(w) {
      warnings_seen <<- c(warnings_seen, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  if (ncol(s$u) != K || ncol(s$v) != K || length(s$d) != K) {
    stop("Yang SSVD did not return exactly K factors.")
  }
  U <- sweep(as.matrix(s$u), 2L, sqrt(as.numeric(s$d)), "*")
  V <- sweep(as.matrix(s$v), 2L, sqrt(as.numeric(s$d)), "*")
  uv <- balance_factor_pairs(U, V)
  Xhat_z <- uv$U %*% t(uv$V)
  name_factor_pairs(list(
    U = uv$U, V = uv$V, Xhat = Xhat_z * bundle$SE,
    K_active = K, input_coordinate = "Z=B/SE",
    reconstruction_coordinate = "raw_B=(U_Z V_Z^T)*SE",
    rank_contract = "fixed_r", singular_values = as.numeric(s$d),
    ssvd_niter = s$niter, ssvd_sigma_hat = s$sigma.hat,
    warnings = unique(warnings_seen)
  ), bundle)
}

# Raw-coordinate Yang et al. sparse SVD.  This is the strict raw-U/V
# counterpart of fit_yang_ssvd_z().
fit_yang_ssvd_raw <- function(bundle, K) {
  if (!requireNamespace("ssvd", quietly = TRUE)) stop("ssvd is not installed.")
  set.seed(as.integer(bundle$meta$seed) + 1731L)
  warnings_seen <- character()
  s <- withCallingHandlers(
    ssvd::ssvd(bundle$B, method = "method", r = K),
    warning = function(w) {
      warnings_seen <<- c(warnings_seen, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  if (ncol(s$u) != K || ncol(s$v) != K || length(s$d) != K) {
    stop("Yang SSVD did not return exactly K factors.")
  }
  U <- sweep(as.matrix(s$u), 2L, sqrt(as.numeric(s$d)), "*")
  V <- sweep(as.matrix(s$v), 2L, sqrt(as.numeric(s$d)), "*")
  uv <- balance_factor_pairs(U, V)
  name_factor_pairs(list(
    U = uv$U, V = uv$V, Xhat = uv$U %*% t(uv$V),
    K_active = K, input_coordinate = "raw_B",
    reconstruction_coordinate = "raw_B=U V^T",
    rank_contract = "fixed_r", singular_values = as.numeric(s$d),
    ssvd_niter = s$niter, ssvd_sigma_hat = s$sigma.hat,
    warnings = unique(warnings_seen)
  ), bundle)
}

fit_factorgo <- function(bundle, K, sample_sizes, executable, work_dir,
                         max_iter = 10000L) {
  sample_sizes <- as.numeric(sample_sizes)
  if (length(sample_sizes) != ncol(bundle$B) || any(!is.finite(sample_sizes)) ||
      any(sample_sizes <= 0)) stop("Invalid FactorGo sample-size vector.")
  if (!file.exists(executable)) stop("FactorGo executable not found: ", executable)
  dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)
  Zscore <- bundle$B / bundle$SE
  z_path <- file.path(work_dir, "input_Z.tsv")
  n_path <- file.path(work_dir, "input_N.tsv")
  prefix <- file.path(work_dir, "factorgo")
  z_tab <- data.frame(SNP = rownames(Zscore), Zscore, check.names = FALSE)
  data.table::fwrite(z_tab, z_path, sep = "\t", quote = FALSE)
  data.table::fwrite(data.frame(N = sample_sizes), n_path, sep = "\t", quote = FALSE)
  log_path <- file.path(work_dir, "factorgo_console.log")
  status <- system2(
    executable,
    args = c(z_path, n_path, "-k", as.character(K),
             "--init-factor", "svd", "--max-iter", as.character(max_iter),
             "--elbo-tol", "0.001", "-p", "cpu", "-s",
             as.character(bundle$meta$seed), "-o", prefix),
    stdout = log_path, stderr = log_path
  )
  if (!identical(status, 0L)) stop("FactorGo exited with status ", status, ".")
  read_gzip_table <- function(path) {
    con <- gzfile(path, open = "rt")
    on.exit(close(con), add = TRUE)
    as.matrix(utils::read.table(
      con, sep = "\t", header = FALSE, check.names = FALSE
    ))
  }
  U <- read_gzip_table(paste0(prefix, ".Wm.tsv.gz"))
  V <- read_gzip_table(paste0(prefix, ".Zm.tsv.gz"))
  if (!identical(dim(U), c(nrow(bundle$B), K)) ||
      !identical(dim(V), c(ncol(bundle$B), K))) {
    stop("Unexpected FactorGo output dimensions.")
  }
  uv <- balance_factor_pairs(U, V)
  factor_info <- read_gzip_table(paste0(prefix, ".factor.tsv.gz"))
  name_factor_pairs(list(
    U = uv$U, V = uv$V, Xhat = uv$U %*% t(uv$V), K_active = K,
    input_coordinate = "Z=B/SE plus supplied trait-specific sample sizes",
    reconstruction_coordinate = "raw_effect=U_variant V_trait^T",
    rank_contract = "K_with_ARD_effective_rank_at_most_K",
    factor_info = factor_info, scale_input = FALSE,
    sample_sizes = sample_sizes,
    init_factor = "svd", max_iter = max_iter, elbo_tol = 0.001
  ), bundle)
}

fit_formal20_benchmark <- function(method, bundle, K, sample_sizes = NULL,
                                   factorgo_executable = NULL,
                                   work_dir = tempdir()) {
  switch(
    method,
    `SVD-raw` = fit_svd_raw(bundle, K),
    SVD = fit_svd_z(bundle, K),
    `SVD-adj` = fit_svd_adj_z(bundle, K),
    `SVD-adj-GLEANR-official` = fit_svd_adj_gleanr_official(bundle, K),
    FactorGo = fit_factorgo(
      bundle, K, sample_sizes, factorgo_executable, work_dir
    ),
    `SSVD-raw` = fit_yang_ssvd_raw(bundle, K),
    SSVD = fit_yang_ssvd_z(bundle, K),
    stop("Unsupported benchmark method: ", method)
  )
}
