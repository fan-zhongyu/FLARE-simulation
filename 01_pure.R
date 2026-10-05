# Pure simulation (Figure 2b): 3 loading architectures x 5 sample sizes x 20
# replicates = 300 data sets. For each data set FLARE, FLARE without the graph
# and the stage 1 solution are fitted and scored against the truth.
#
# Usage: Rscript 01_pure.R TASK        TASK = 1..60, five data sets each
source(file.path(Sys.getenv("FLARE_SIMULATION", "."), "R/common.R"))
task <- as.integer(commandArgs(TRUE)[1])
stopifnot(task >= 1L, task <= 60L)
load_scoring()

tasks <- fread(PURE_TASKS)[family == "pure"]
stopifnot(nrow(tasks) == 300L)
limit <- as.integer(Sys.getenv("FLARE_TEST_LIMIT", "5"))   # fewer data sets for a test run

for (i in ((task - 1L) * 5L + 1:5)[seq_len(limit)]) {
  job <- tasks[i]
  out <- out_path("results/pure", paste0(job$run_id, ".tsv"))
  if (file.exists(out)) next
  bundle <- readRDS(job$bundle_path)
  graph <- load_graph_pair(bundle)
  K <- ncol(bundle$U_true)

  started <- Sys.time()
  fits <- flare_variants(bundle$B, bundle$SE, graph$A_true, C = bundle$C,
                         C_se = bundle$C_se, sample_sd = bundle$sample_sd,
                         K = K, control = CONTROL$pure)
  seconds <- as.numeric(difftime(Sys.time(), started, units = "secs"))

  rows <- lapply(names(fits), function(m) {
    x <- fits[[m]]
    meta <- c(metric_meta(bundle, job$run_id, i, "FLARE", m),
              list(method_label = m, architecture = job$architecture,
                   condition = job$condition, converged = x$fit$converged,
                   seconds = seconds, FLARE_version = as.character(packageVersion("FLARE"))))
    saveRDS(x[c("U", "V")], out_path("fits/pure", sprintf("%s__%s.rds", job$run_id, gsub("[^A-Za-z0-9]+", "_", m))))
    score_all(bundle, x$U, x$V, meta)
  })
  fwrite(rbindlist(rows, fill = TRUE), out, sep = "\t")
  cat("PURE_DONE", job$run_id, "\n")
}
