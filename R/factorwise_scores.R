# Factor-by-factor R2 after matching.

formal20_factorwise_scores <- function(true_U, true_V, pred_U, pred_V) {
  true_U <- as.matrix(true_U); true_V <- as.matrix(true_V)
  pred_U <- as.matrix(pred_U); pred_V <- as.matrix(pred_V)
  K <- ncol(true_V)
  if (ncol(pred_U) < K) {
    pred_U <- cbind(pred_U, matrix(0, nrow(pred_U), K - ncol(pred_U)))
  }
  if (ncol(pred_V) < K) {
    pred_V <- cbind(pred_V, matrix(0, nrow(pred_V), K - ncol(pred_V)))
  }
  if (ncol(pred_U) != ncol(pred_V)) stop("Predicted U/V ranks differ.")
  if (ncol(pred_U) > K) {
    aligned <- gleanr_exact_alignment(true_U, true_V, pred_U, pred_V)
    pred_U <- aligned$U
    pred_V <- aligned$V
  } else {
    pred_U <- pred_U[, seq_len(K), drop = FALSE]
    pred_V <- pred_V[, seq_len(K), drop = FALSE]
  }
  cor2 <- function(A, B) {
    out <- matrix(0, ncol(A), ncol(B))
    for (i in seq_len(ncol(A))) for (j in seq_len(ncol(B))) {
      out[i, j] <- safe_cor(A[, i], B[, j])^2
    }
    out
  }
  cu <- cor2(true_U, pred_U)
  cv <- cor2(true_V, pred_V)
  perms <- permutations_int(seq_len(K))
  score_perm <- function(M, p) mean(M[cbind(seq_len(K), p)])
  paired_values <- vapply(
    perms, function(p) score_perm(cu + cv, p), numeric(1L)
  )
  v_values <- vapply(perms, function(p) score_perm(cv, p), numeric(1L))
  paired <- perms[[which.max(paired_values)]]
  vonly <- perms[[which.max(v_values)]]
  c(
    factorwise_paired_R2_U = score_perm(cu, paired),
    factorwise_paired_R2_V = score_perm(cv, paired),
    factorwise_Vmatched_R2_U = score_perm(cu, vonly),
    factorwise_Vmatched_R2_V = score_perm(cv, vonly)
  )
}

append_formal20_scores <- function(row, bundle, U, V) {
  z <- as.list(formal20_factorwise_scores(
    bundle$U_true, bundle$V_true, U, V
  ))
  for (nm in names(z)) row[[nm]] <- z[[nm]]
  row
}
