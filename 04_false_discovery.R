# False-discovery designs (Section 3.4): 660 data sets in three designs,
#   A  global null (no factor at all),
#   B  null traits (some traits carry no signal),
#   C  rank set above the truth.
# For each data set: the null metrics of the benchmark and the number of
# factors that pass the screen, for FLARE and for its stage 1 solution.
#
# Usage: Rscript 04_false_discovery.R TASK        TASK = 1..66, ten data sets each
source(file.path(Sys.getenv("FLARE_SIMULATION", "."), "R/common.R"))
source(file.path(REPO, "R/null_metrics.R"))
task <- as.integer(commandArgs(TRUE)[1])
stopifnot(task >= 1L, task <= 66L)
load_scoring()

grid <- fread(FDR_GRID)[method == "FLARE"]
stopifnot(nrow(grid) == 660L)
limit <- as.integer(Sys.getenv("FLARE_TEST_LIMIT", "10"))

for (i in ((task - 1L) * 10L + 1:10)[seq_len(limit)]) {
  job <- grid[i]
  out <- out_path("results/false_discovery", job$run_id, "metrics.tsv")
  if (file.exists(out)) next
  bundle <- readRDS(job$bundle)
  graph <- load_graph_pair(bundle)
  K <- as.integer(job$K)
  reference <- noise_reference(nrow(bundle$B), ncol(bundle$B))

  started <- Sys.time()
  fit <- flare_fit(bundle$B, bundle$SE, bundle$C, C_se = bundle$C_se, sample_sd = bundle$sample_sd,
                   graph = graph$A_true, K = K, noise_reference = reference,
                   control = CONTROL$null, verbose = FALSE)
  seconds <- as.numeric(difftime(Sys.time(), started, units = "secs"))
  variants <- list(
    FLARE = list(U = fit$U, V = fit$V, screen = fit$screen),
    `FLARE, no graph` = list(U = fit$stage1$U, V = fit$stage1$V,
      screen = flare_screen(bundle$B, bundle$SE, fit$stage1$V, fit$covariance$C, noise_reference = reference))
  )

  meta <- list(run_id = job$run_id, scenario_id = job$scenario_id, design = job$design,
               v_architecture = job$v_architecture, overlap_rho = job$overlap_rho,
               sample_size_label = job$sample_size_label, replicate_id = job$replicate_id,
               K = K, elapsed_seconds = seconds)
  null_traits <- if (identical(job$design, "null_traits")) bundle$null_traits else character()

  metrics <- rbindlist(lapply(names(variants), function(label) {
    x <- variants[[label]]
    V <- x$V
    rownames(V) <- colnames(bundle$B)
    m <- null_metrics(bundle$B, bundle$SE, x$U, V, null_traits = null_traits)
    m$factors_kept <- sum(x$screen$keep)
    m$screen_threshold <- x$screen$threshold[1]
    m$screen_largest <- max(x$screen$energy)
    if (identical(job$design, "overspecified_rank")) {
      # Which fitted factors correspond to true ones, and what the extra
      # columns look like (empty, single-trait or multi-trait).
      al <- gleanr_exact_alignment(bundle$U_true, bundle$V_true, x$U, V)
      energy <- colSums(x$U^2) * colSums(V^2)
      extra <- intersect(al$unmatched, active_factors(x$U, V))
      top <- vapply(extra, function(k) max(abs(V[, k])) / sqrt(sum(V[, k]^2)), numeric(1L))
      kept <- which(x$screen$keep)
      m <- c(m, list(
        K_truth = ncol(bundle$U_true), K_extra_active = length(extra),
        K_extra_multi_trait = sum(top <= 0.98), K_extra_single_trait = sum(top > 0.98),
        extra_factor_energy = if (sum(energy) > 0) sum(energy[extra]) / sum(energy) else 0,
        U_R2_matched = al$U_r2, V_R2_matched = al$V_r2,
        kept_true_factors = sum(kept %in% al$order[seq_len(ncol(bundle$U_true))]),
        kept_extra_factors = sum(kept %in% extra)))
    }
    if (identical(job$design, "null_traits")) {
      al <- gleanr_exact_alignment(bundle$U_true, bundle$V_true, x$U, V)
      m$U_R2_matched <- al$U_r2
      m$V_R2_matched <- al$V_r2
    }
    as.data.table(c(meta, list(method_label = label, status = "ok"), m))
  }), fill = TRUE)

  saveRDS(lapply(variants, function(x) x[c("U", "V")]), out_path("results/false_discovery", job$run_id, "fit.rds"))
  fwrite(metrics, out, sep = "\t")
  cat("FDR_DONE", job$run_id, "\n")
}
