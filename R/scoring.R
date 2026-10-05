# Scores of a factorization against the simulated truth: matching of estimated to true
# factors, R2 of U and V, subspace and support measures, recovery of the effect matrix.

solver_transform <- function(X, SE, C) {
  Wc <- solve(chol(C))
  t(Wc) %*% t(X / SE)
}

safe_cor <- function(x, y) {
  x <- as.numeric(x)
  y <- as.numeric(y)
  if (length(x) != length(y) || length(x) < 2L) return(0)
  if (any(!is.finite(x)) || any(!is.finite(y))) return(0)
  if (stats::sd(x) <= 0 || stats::sd(y) <= 0) return(0)
  as.numeric(stats::cor(x, y))
}

permutations_int <- function(x) {
  if (length(x) == 1L) return(list(x))
  out <- vector("list", factorial(length(x)))
  at <- 1L
  for (i in seq_along(x)) {
    tails <- permutations_int(x[-i])
    for (z in tails) {
      out[[at]] <- c(x[[i]], z)
      at <- at + 1L
    }
  }
  out
}

# Maximum-weight one-to-one matching of K truth factors to P >= K estimated
# factors. Dynamic programming is exact for the formal K=6, P=12 experiment
# and avoids adding a package dependency for rectangular assignment.
max_weight_factor_assignment <- function(score) {
  score <- as.matrix(score)
  K <- nrow(score); P <- ncol(score)
  if (K < 1L || P < K) stop("Assignment requires P >= K >= 1.")
  if (P > 20L) {
    # Exact rectangular Hungarian assignment. This branch supports the
    # trait-dimension-capped K_init=24 model-selection experiments without the
    # 2^P memory growth of the historical bitmask implementation.
    cost <- max(score) - score
    nr <- nrow(cost)
    nc <- ncol(cost)
    u <- numeric(nr)
    v <- numeric(nc + 1L)
    p <- integer(nc + 1L)
    way <- integer(nc + 1L)
    for (ii in seq_len(nr)) {
      p[[1L]] <- ii
      j0 <- 1L
      minv <- rep(Inf, nc + 1L)
      used <- rep(FALSE, nc + 1L)
      repeat {
        used[[j0]] <- TRUE
        i0 <- p[[j0]]
        delta <- Inf
        j1 <- 0L
        for (jj in 2:(nc + 1L)) {
          if (!used[[jj]]) {
            cur <- cost[i0, jj - 1L] - u[[i0]] - v[[jj]]
            if (cur < minv[[jj]]) {
              minv[[jj]] <- cur
              way[[jj]] <- j0
            }
            if (minv[[jj]] < delta) {
              delta <- minv[[jj]]
              j1 <- jj
            }
          }
        }
        for (jj in seq_len(nc + 1L)) {
          if (used[[jj]]) {
            u[[p[[jj]]]] <- u[[p[[jj]]]] + delta
            v[[jj]] <- v[[jj]] - delta
          } else {
            minv[[jj]] <- minv[[jj]] - delta
          }
        }
        j0 <- j1
        if (p[[j0]] == 0L) break
      }
      repeat {
        j1 <- way[[j0]]
        p[[j0]] <- p[[j1]]
        j0 <- j1
        if (j0 == 1L) break
      }
    }
    assignment <- integer(nr)
    for (jj in 2:(nc + 1L)) {
      if (p[[jj]] > 0L) assignment[[p[[jj]]]] <- jj - 1L
    }
    return(list(
      order = assignment,
      score = sum(score[cbind(seq_len(K), assignment)])
    ))
  }
  nstate <- bitwShiftL(1L, P)
  dp <- rep(-Inf, nstate)
  dp[[1L]] <- 0
  parent_mask <- matrix(NA_integer_, K, nstate)
  parent_col <- matrix(NA_integer_, K, nstate)
  active <- 0L
  for (i in seq_len(K)) {
    next_dp <- rep(-Inf, nstate)
    next_active <- integer()
    for (mask in active) {
      base <- dp[[mask + 1L]]
      for (j in seq_len(P)) {
        bit <- bitwShiftL(1L, j - 1L)
        if (bitwAnd(mask, bit) != 0L) next
        new_mask <- bitwOr(mask, bit)
        value <- base + score[i, j]
        idx <- new_mask + 1L
        if (value > next_dp[[idx]]) {
          if (!is.finite(next_dp[[idx]])) next_active <- c(next_active, new_mask)
          next_dp[[idx]] <- value
          parent_mask[i, idx] <- mask
          parent_col[i, idx] <- j
        }
      }
    }
    dp <- next_dp
    active <- unique(next_active)
  }
  final_mask <- active[[which.max(dp[active + 1L])]]
  assignment <- integer(K)
  mask <- final_mask
  for (i in K:1L) {
    idx <- mask + 1L
    assignment[[i]] <- parent_col[i, idx]
    mask <- parent_mask[i, idx]
  }
  list(order = assignment, score = dp[[final_mask + 1L]])
}

# Exact GLEANR/Yuan criterion for K <= 6, implemented with sufficient
# statistics to avoid materializing 46,080 full candidate matrices.
gleanr_exact_alignment <- function(true_U, true_V, pred_U, pred_V) {
  true_U <- as.matrix(true_U); true_V <- as.matrix(true_V)
  pred_U <- as.matrix(pred_U); pred_V <- as.matrix(pred_V)
  K <- ncol(true_U)
  if (ncol(true_V) != K) stop("Truth U/V ranks differ.")
  if (ncol(pred_U) < K) pred_U <- cbind(pred_U, matrix(0, nrow(pred_U), K - ncol(pred_U)))
  if (ncol(pred_V) < K) pred_V <- cbind(pred_V, matrix(0, nrow(pred_V), K - ncol(pred_V)))
  if (ncol(pred_U) != ncol(pred_V)) stop("Predicted U/V ranks differ.")
  P <- ncol(pred_U)
  pred_U[!is.finite(pred_U)] <- 0
  pred_V[!is.finite(pred_V)] <- 0
  if (all(pred_U == 0) || all(pred_V == 0)) {
    ord <- seq_len(K)
    return(list(
      U_r2 = 0, V_r2 = 0,
      U = pred_U[, ord, drop = FALSE],
      V = pred_V[, ord, drop = FALSE],
      order = ord, signs = rep(1, K), criterion = 0,
      selected_cor_U = 0, selected_cor_V = 0,
      unmatched = setdiff(seq_len(P), ord),
      alignment_method = "zero_prediction"
    ))
  }
  if (P > K) {
    pair_cor <- function(T, Pm) {
      out <- matrix(0, ncol(T), ncol(Pm))
      for (i in seq_len(ncol(T))) for (j in seq_len(ncol(Pm))) {
        out[i, j] <- safe_cor(T[, i], Pm[, j])
      }
      out
    }
    cu <- pair_cor(true_U, pred_U)
    cv <- pair_cor(true_V, pred_V)
    matched <- max_weight_factor_assignment((cu^2 + cv^2) / 2)
    ord <- matched$order
    sg <- ifelse(cu[cbind(seq_len(K), ord)] +
                   cv[cbind(seq_len(K), ord)] < 0, -1, 1)
    aligned_U <- sweep(pred_U[, ord, drop = FALSE], 2L, sg, "*")
    aligned_V <- sweep(pred_V[, ord, drop = FALSE], 2L, sg, "*")
    return(list(
      U_r2 = safe_cor(true_U, aligned_U)^2,
      V_r2 = safe_cor(true_V, aligned_V)^2,
      U = aligned_U, V = aligned_V, order = ord, signs = sg,
      criterion = matched$score,
      selected_cor_U = mean(cu[cbind(seq_len(K), ord)] * sg),
      selected_cor_V = mean(cv[cbind(seq_len(K), ord)] * sg),
      unmatched = setdiff(seq_len(P), ord),
      alignment_method = "exact_rectangular_factorwise_assignment"
    ))
  }

  signs <- as.matrix(expand.grid(replicate(K, c(1, -1), simplify = FALSE)))
  perms <- permutations_int(seq_len(K))

  prep <- function(T, P) {
    nall <- length(T)
    st <- sum(T); sp2 <- sum(P^2)
    vt <- sum(T^2) - st^2 / nall
    list(n = nall, st = st, vt = max(vt, 0), sp2 = sp2,
         colsum = colSums(P), cross = crossprod(T, P))
  }
  u <- prep(true_U, pred_U)
  v <- prep(true_V, pred_V)
  correlations <- function(z, perm) {
    sp <- as.numeric(signs %*% z$colsum[perm])
    matched_cross <- z$cross[cbind(seq_len(K), perm)]
    cp <- as.numeric(signs %*% matched_cross)
    vp <- pmax(z$sp2 - sp^2 / z$n, 0)
    den <- sqrt(z$vt * vp)
    r <- (cp - z$st * sp / z$n) / den
    r[!is.finite(r)] <- 0
    pmax(-1, pmin(1, r))
  }

  best <- -Inf; best_perm <- seq_len(K); best_sign <- rep(1, K)
  best_ru <- 0; best_rv <- 0
  for (perm in perms) {
    ru <- correlations(u, perm)
    rv <- correlations(v, perm)
    score <- (ru + rv)^2
    ii <- which.max(score)
    if (score[[ii]] > best) {
      best <- score[[ii]]
      best_perm <- perm
      best_sign <- as.numeric(signs[ii, ])
      best_ru <- ru[[ii]]
      best_rv <- rv[[ii]]
    }
  }
  aligned_U <- sweep(pred_U[, best_perm, drop = FALSE], 2L, best_sign, "*")
  aligned_V <- sweep(pred_V[, best_perm, drop = FALSE], 2L, best_sign, "*")
  list(U_r2 = safe_cor(true_U, aligned_U)^2,
       V_r2 = safe_cor(true_V, aligned_V)^2,
       U = aligned_U, V = aligned_V, order = best_perm,
       signs = best_sign, criterion = best,
       selected_cor_U = best_ru, selected_cor_V = best_rv,
       unmatched = integer(), alignment_method = "exact_K_factor_alignment")
}

orthobasis <- function(X, tol = 1e-10) {
  s <- svd(as.matrix(X), nu = min(dim(X)), nv = 0)
  if (!length(s$d) || max(s$d) <= 0) return(matrix(0, nrow(X), 0L))
  r <- sum(s$d > tol * max(s$d))
  if (!r) return(matrix(0, nrow(X), 0L))
  s$u[, seq_len(r), drop = FALSE]
}

subspace_metrics <- function(T, P, prefix) {
  qt <- orthobasis(T); qp <- orthobasis(P)
  rt <- ncol(qt); rp <- ncol(qp)
  if (!rt || !rp) {
    overlap <- 0; min_cos2 <- 0; mean_cos2 <- 0
  } else {
    cs <- svd(crossprod(qt, qp), nu = 0, nv = 0)$d
    cs <- pmin(1, pmax(0, cs))
    overlap <- sum(cs^2) / max(rt, rp)
    mean_cos2 <- mean(cs^2)
    min_cos2 <- min(cs^2)
  }
  proj_dist <- sqrt(max(rt + rp - 2 * overlap * max(rt, rp), 0))
  stats::setNames(c(rt, rp, overlap, mean_cos2, min_cos2, proj_dist),
                  paste0(prefix, c("_rank_true", "_rank_hat", "_subspace_overlap",
                                   "_principal_cos2_mean", "_principal_cos2_min",
                                   "_projector_frobenius")))
}

support_metrics <- function(T, P) {
  K <- ncol(T)
  one <- lapply(seq_len(K), function(k) {
    truth <- T[, k] != 0
    pred <- abs(P[, k]) > 1e-10
    tp <- sum(truth & pred); fp <- sum(!truth & pred); fn <- sum(truth & !pred)
    precision <- if (tp + fp) tp / (tp + fp) else 0
    recall <- if (tp + fn) tp / (tp + fn) else 0
    f1 <- if (precision + recall) 2 * precision * recall / (precision + recall) else 0
    c(precision = precision, recall = recall, f1 = f1)
  })
  m <- do.call(rbind, one)
  c(U_support_precision = mean(m[, "precision"]),
    U_support_recall = mean(m[, "recall"]),
    U_support_f1 = mean(m[, "f1"]))
}

score_factorization <- function(bundle, U, V, fit_meta = list(),
                                Xhat_override = NULL) {
  U <- as.matrix(U); V <- as.matrix(V)
  TU <- bundle$U_true; TV <- bundle$V_true
  al <- gleanr_exact_alignment(TU, TV, U, V)
  Xhat <- if (is.null(Xhat_override)) U %*% t(V) else as.matrix(Xhat_override)
  Xtrue <- bundle$X_true
  if (!identical(dim(Xhat), dim(Xtrue))) {
    stop("Xhat dimensions differ from X_true.")
  }
  if (any(!is.finite(Xhat))) stop("Xhat contains non-finite values.")
  TXhat <- solver_transform(Xhat, bundle$SE, bundle$C)
  TXtrue <- solver_transform(Xtrue, bundle$SE, bundle$C)
  TY <- solver_transform(bundle$B, bundle$SE, bundle$C)
  raw_rmse <- sqrt(mean((Xhat - Xtrue)^2))
  solver_rmse <- sqrt(mean((TXhat - TXtrue)^2))
  truth_raw_rms <- sqrt(mean(Xtrue^2))
  truth_solver_rms <- sqrt(mean(TXtrue^2))

  un_t <- sqrt(colSums(TU^2)); vn_t <- sqrt(colSums(TV^2))
  un_h <- sqrt(colSums(al$U^2)); vn_h <- sqrt(colSums(al$V^2))
  safe_ratio <- function(a, b) ifelse(b > 0, a / b, NA_real_)
  rank_ratio <- safe_ratio(un_h * vn_h, un_t * vn_t)
  imbalance_t <- log10((un_t + 1e-300) / (vn_t + 1e-300))
  imbalance_h <- log10((un_h + 1e-300) / (vn_h + 1e-300))
  pred_un <- sqrt(colSums(U^2)); pred_vn <- sqrt(colSums(V^2))
  pred_strength <- pred_un * pred_vn
  pred_energy <- pred_strength^2
  strength_max <- if (length(pred_strength)) max(pred_strength) else 0
  strength_active <- if (is.finite(strength_max) && strength_max > 0) {
    pred_strength / strength_max > 1e-3
  } else {
    rep(FALSE, length(pred_strength))
  }
  unmatched_energy <- if (length(al$unmatched)) {
    sum(pred_energy[al$unmatched])
  } else 0
  total_energy <- sum(pred_energy)
  spurious_fraction <- if (total_energy > 0) unmatched_energy / total_energy else 0

  vals <- c(
    gleanr_R2_U = al$U_r2,
    gleanr_R2_V = al$V_r2,
    primary_mean_R2 = mean(c(al$U_r2, al$V_r2)),
    support_metrics(TU, al$U),
    subspace_metrics(TU, U, "U"),
    subspace_metrics(TV, V, "V"),
    X_pearson = safe_cor(Xtrue, Xhat),
    X_R2_correlation = safe_cor(Xtrue, Xhat)^2,
    X_raw_rmse = raw_rmse,
    X_raw_relative_rmse = raw_rmse / max(truth_raw_rms, 1e-300),
    X_solver_pearson = safe_cor(TXtrue, TXhat),
    X_solver_rmse = solver_rmse,
    X_solver_relative_rmse = solver_rmse / max(truth_solver_rms, 1e-300),
    fit_model_pve = 1 - sum((TY - TXhat)^2) / sum(TY^2),
    oracle_realized_pve = 1 - sum((TY - TXtrue)^2) / sum(TY^2),
    U_frobenius = sqrt(sum(U^2)),
    V_frobenius = sqrt(sum(V^2)),
    U_sparsity = mean(abs(U) <= 1e-10),
    V_sparsity = mean(abs(V) <= 1e-10),
    U_norm_ratio_geomean = exp(mean(log(pmax(safe_ratio(un_h, un_t), 1e-300)))),
    V_norm_ratio_geomean = exp(mean(log(pmax(safe_ratio(vn_h, vn_t), 1e-300)))),
    rank1_strength_ratio_geomean = exp(mean(log(pmax(rank_ratio, 1e-300)))),
    rank1_strength_log10_rmse = sqrt(mean(log10(pmax(rank_ratio, 1e-300))^2)),
    UV_log10_imbalance_rmse = sqrt(mean((imbalance_h - imbalance_t)^2)),
    UV_log10_imbalance_hat_mean = mean(imbalance_h),
    K_hat = ncol(U),
    K_effective_strength_1e3 = sum(strength_active),
    spurious_factor_energy_fraction = spurious_fraction,
    matched_factor_energy_fraction = if (total_energy > 0) 1 - spurious_fraction else 0
  )
  row <- as.data.frame(as.list(vals), check.names = FALSE)
  for (nm in names(fit_meta)) row[[nm]] <- fit_meta[[nm]]
  row$alignment_order <- paste(al$order, collapse = ",")
  row$alignment_signs <- paste(al$signs, collapse = ",")
  row$alignment_method <- al$alignment_method
  row$unmatched_factor_columns <- paste(al$unmatched, collapse = ",")
  row$factor_strengths <- paste(signif(pred_strength, 8), collapse = ",")
  row$U_column_norm_ratios <- paste(signif(safe_ratio(un_h, un_t), 8), collapse = ",")
  row$V_column_norm_ratios <- paste(signif(safe_ratio(vn_h, vn_t), 8), collapse = ",")
  row$rank1_strength_ratios <- paste(signif(rank_ratio, 8), collapse = ",")
  row
}
