# Null GWAS simulation (Figure 2a): 840 data sets with no genetic factor, built
# from permuted phenotypes with varying sample overlap (the design of the
# GLEANR paper). A method should return no factor. FLARE is fitted with K = 8
# and the number of factors that pass the screen is recorded, for the fit with
# the graph and for the graph-free stage 1 solution.
#
# Usage: Rscript 03_null_gwas.R TASK        TASK = 1..20, 42 data sets each
source(file.path(Sys.getenv("FLARE_SIMULATION", "."), "R/common.R"))
task <- as.integer(commandArgs(TRUE)[1])
stopifnot(task >= 1L, task <= 20L)

sims <- file.path(NULL_GWAS, "sims")
stems <- sort(sub("\\.BETA\\.csv$", "", list.files(sims, pattern = "\\.BETA\\.csv$")))
stopifnot(length(stems) == 840L)
read_matrix <- function(stem, suffix) {
  m <- as.matrix(fread(file.path(sims, paste0(stem, suffix))))
  if (nrow(m) == ncol(m)) rownames(m) <- colnames(m)
  m
}
with_ids <- function(m) {
  rownames(m) <- paste0("rs", seq_len(nrow(m)))
  m
}
graph <- readRDS(file.path(NULL_GWAS, "input", "annotation_graph_first1000.rds"))
K <- 8L
limit <- as.integer(Sys.getenv("FLARE_TEST_LIMIT", "42"))

for (stem in stems[((task - 1L) * 42L + 1:42)[seq_len(limit)]]) {
  out <- out_path("results/null_gwas", paste0(stem, ".tsv"))
  if (file.exists(out)) next
  started <- Sys.time()
  B <- with_ids(read_matrix(stem, ".BETA.csv"))
  SE <- with_ids(read_matrix(stem, ".SE.csv"))
  C_raw <- read_matrix(stem, ".COV.csv")
  C_se <- read_matrix(stem, ".COV_SE.csv")
  # Independent cohorts have an identity matrix; otherwise shrink the
  # estimated correlation (no block step: the nine traits form one block).
  independent <- max(abs(C_raw - diag(ncol(B)))) < 1e-12
  C <- if (independent) NULL else shrink_only(C_raw, C_se)
  reference <- noise_reference(nrow(B), ncol(B))

  fit <- flare_fit(B, SE, C, graph = graph, K = K, noise_reference = reference,
                   control = CONTROL$null, verbose = FALSE)
  stage1 <- flare_screen(B, SE, fit$stage1$V, fit$covariance$C, noise_reference = reference)

  parts <- strsplit(stem, ".", fixed = TRUE)[[1]]
  count_active <- function(U, V) sum(colSums(V != 0) > 0 & colSums(U != 0) > 0)
  row <- data.table(
    case = stem, seed = parts[1], overlap = parts[2], replicate = parts[3],
    FLARE = length(fit$keep), FLARE_largest = max(fit$screen$energy),
    FLARE_K_fitted = count_active(fit$U, fit$V),
    `FLARE-NG` = sum(stage1$keep), `FLARE-NG_largest` = max(stage1$energy),
    `FLARE-NG_K_fitted` = count_active(fit$stage1$U, fit$stage1$V),
    threshold = fit$screen$threshold[1],
    seconds = as.numeric(difftime(Sys.time(), started, units = "secs"))
  )
  saveRDS(list(FLARE = fit$V, `FLARE-NG` = fit$stage1$V, C = fit$covariance$C, unit = fit$unit),
          out_path("fits/null_gwas", paste0(stem, ".rds")))
  fwrite(row, out, sep = "\t")
  cat("NULL_GWAS_DONE", stem, "\n")
}
