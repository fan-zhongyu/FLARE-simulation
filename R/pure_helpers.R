# Helpers of the pure simulation: the graph of a data set and the labels written with every score.

load_graph_pair <- function(bundle) {
  graph <- readRDS(normalizePath(bundle$graph_path, mustWork = TRUE))
  ids <- rownames(bundle$B)
  for (name in c("A_true", "A_permuted")) {
    A <- Matrix::Matrix(graph[[name]], sparse = TRUE)
    if (!identical(rownames(A), ids) || !identical(colnames(A), ids)) {
      stop(name, " SNP identifiers do not match the bundle.")
    }
    if (max(abs(A - Matrix::t(A))) > 1e-10) {
      stop(name, " is asymmetric.")
    }
    graph[[name]] <- A
  }
  graph
}

metric_meta <- function(bundle, run_id, task_id, method, graph_arm) {
  meta <- bundle$meta
  list(
    run_id = run_id,
    task_id = as.integer(task_id),
    stage = "pure_v_architecture_formal20_mixedN",
    experiment_id = "flare_pure_v_architecture_formal20_mixedN_v1",
    seed = meta$seed,
    replicate_id = meta$replicate_id,
    scenario_id = meta$scenario_id,
    v_architecture = meta$v_architecture,
    v_architecture_label = meta$v_architecture_label,
    c_condition = meta$c_condition,
    sample_size = meta$sample_size,
    sample_size_label = meta$sample_size_label,
    n_condition = meta$n_condition,
    sample_size_design = meta$sample_size_design,
    sample_size_min = meta$sample_size_min,
    sample_size_median = meta$sample_size_median,
    sample_size_max = meta$sample_size_max,
    sample_size_geomean = meta$sample_size_geomean,
    w2 = meta$w2,
    graph_k = meta$graph_k,
    signal_level = meta$sample_size_label,
    method = method,
    graph_arm = graph_arm,
    truth_rank = meta$truth_rank,
    trait_h2_target = meta$trait_h2_target,
    truth_V_zero_fraction = meta$truth_V_zero_fraction,
    truth_V_dense_factors = meta$truth_V_dense_factors,
    covariance_input = "noisy_C_hat",
    covariance_truth_available_to_method = FALSE,
    C_se_offdiag = meta$C_se_offdiag,
    graph_truth_contract = "same_K5_graph_generation_and_fit",
    annotation_contract = "synthetic_scale_IDF_L2_weighted_cosine_OR",
    elapsed_seconds = NA_real_
  )
}
