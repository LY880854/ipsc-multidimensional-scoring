args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop(
    "Usage: Rscript train_esc_similarity.R ",
    "GTEX_COUNTS_CSV ESC_MANIFEST_TSV OUTPUT_DIR"
  )
}

gtex_file <- args[[1]]
manifest_file <- args[[2]]
output_dir <- args[[3]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

required_packages <- c("edgeR", "glmnet", "limma")
available <- vapply(required_packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))
if (!all(available)) {
  stop("Missing packages: ", paste(names(available)[!available], collapse = ", "))
}

read_counts <- function(path) {
  input <- utils::read.csv(
    path,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  if (ncol(input) < 2L) {
    stop("Count files must contain a gene column and at least one sample: ", path)
  }
  gene_id <- as.character(input[[1]])
  if (anyNA(gene_id) || any(gene_id == "") || anyDuplicated(gene_id)) {
    stop("Gene identifiers must be present and unique: ", path)
  }
  values <- as.matrix(input[, -1, drop = FALSE])
  storage.mode(values) <- "double"
  rownames(values) <- gene_id
  if (any(!is.finite(values)) || any(values < 0) || any(abs(values - round(values)) > 1e-8)) {
    stop("Count matrices must contain finite, non-negative integer counts: ", path)
  }
  if (!ncol(values) || anyDuplicated(colnames(values))) {
    stop("Sample identifiers must be present and unique: ", path)
  }
  values
}

manifest <- utils::read.delim(
  manifest_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
required_manifest <- c("dataset", "path")
missing_manifest <- setdiff(required_manifest, colnames(manifest))
if (length(missing_manifest)) {
  stop("Missing manifest columns: ", paste(missing_manifest, collapse = ", "))
}
if (!nrow(manifest) || anyDuplicated(manifest$dataset) || any(manifest$dataset == "")) {
  stop("The ESC manifest must contain unique, non-empty dataset names")
}
if (any(!file.exists(manifest$path))) {
  stop("Missing ESC count files: ", paste(manifest$path[!file.exists(manifest$path)], collapse = ", "))
}

esc_matrices <- setNames(lapply(manifest$path, read_counts), manifest$dataset)
gtex_counts <- read_counts(gtex_file)
all_sample_ids <- c(unlist(lapply(esc_matrices, colnames)), colnames(gtex_counts))
if (anyDuplicated(all_sample_ids)) {
  stop("Sample identifiers must be unique across ESC and GTEx inputs")
}

esc_common_genes <- Reduce(intersect, lapply(esc_matrices, rownames))
if (!length(esc_common_genes)) {
  stop("No gene identifiers are shared by all ESC datasets")
}
esc_counts <- do.call(
  cbind,
  lapply(esc_matrices, function(matrix) matrix[esc_common_genes, , drop = FALSE])
)
esc_batch <- rep(
  names(esc_matrices),
  times = vapply(esc_matrices, ncol, FUN.VALUE = integer(1))
)

esc_dge <- edgeR::DGEList(counts = esc_counts)
esc_dge <- edgeR::calcNormFactors(esc_dge, method = "TMM")
esc_logcpm <- edgeR::cpm(esc_dge, log = TRUE, prior.count = 1)
esc_corrected <- limma::removeBatchEffect(esc_logcpm, batch = esc_batch)

gtex_dge <- edgeR::DGEList(counts = gtex_counts)
gtex_dge <- edgeR::calcNormFactors(gtex_dge, method = "TMM")
gtex_logcpm <- edgeR::cpm(gtex_dge, log = TRUE, prior.count = 1)

features <- intersect(rownames(esc_corrected), rownames(gtex_logcpm))
if (!length(features)) {
  stop("No gene identifiers are shared by the ESC and GTEx matrices")
}
x <- t(cbind(
  esc_corrected[features, , drop = FALSE],
  gtex_logcpm[features, , drop = FALSE]
))
x_scaled <- scale(x)
center <- attr(x_scaled, "scaled:center")
scale_vector <- attr(x_scaled, "scaled:scale")
if (any(!is.finite(x_scaled)) || any(!is.finite(center)) || any(!is.finite(scale_vector))) {
  stop("Centering and scaling produced non-finite values")
}

y <- factor(
  c(rep("ESC", ncol(esc_counts)), rep("GTEx", ncol(gtex_counts))),
  levels = c("GTEx", "ESC")
)
set.seed(42)
fit <- glmnet::cv.glmnet(
  x = x_scaled,
  y = y,
  alpha = 0,
  family = "binomial"
)

model_bundle <- list(
  model = fit,
  center = center,
  scale = scale_vector,
  features = colnames(x_scaled),
  preprocessing = list(
    esc_normalization = "TMM; edgeR logCPM with prior.count = 1",
    esc_batch_adjustment = "limma removeBatchEffect across source datasets",
    gtex_normalization = "TMM; edgeR logCPM with prior.count = 1",
    model = "ridge logistic regression",
    alpha = 0,
    seed = 42
  )
)
saveRDS(model_bundle, file.path(output_dir, "esc_similarity_model.rds"))

training_scores <- as.numeric(
  stats::predict(fit, newx = x_scaled, s = "lambda.min", type = "link")
)
sample_metadata <- data.frame(
  sample_id = rownames(x_scaled),
  group = as.character(y),
  dataset = c(esc_batch, rep("GTEx", ncol(gtex_counts))),
  esc_similarity_score = training_scores,
  check.names = FALSE
)
utils::write.table(
  sample_metadata,
  file.path(output_dir, "esc_similarity_training_scores.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
utils::write.table(
  data.frame(feature = colnames(x_scaled)),
  file.path(output_dir, "esc_similarity_features.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
utils::write.table(
  data.frame(
    esc_samples = ncol(esc_counts),
    gtex_samples = ncol(gtex_counts),
    retained_features = ncol(x_scaled),
    lambda_min = fit$lambda.min,
    lambda_1se = fit$lambda.1se,
    check.names = FALSE
  ),
  file.path(output_dir, "esc_similarity_model_summary.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
