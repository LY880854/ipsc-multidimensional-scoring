args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop(
    "Usage: Rscript train_and_evaluate_naive_likeness.R ",
    "RAW_COUNT_TSV SAMPLE_METADATA_TSV OUTPUT_DIR"
  )
}

count_file <- args[[1]]
metadata_file <- args[[2]]
output_dir <- args[[3]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

required_packages <- c("edgeR", "glmnet", "pROC")
available <- vapply(required_packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))
if (!all(available)) {
  stop("Missing packages: ", paste(names(available)[!available], collapse = ", "))
}

counts_input <- utils::read.delim(
  count_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
if (!ncol(counts_input) || colnames(counts_input)[[1]] != "gene_id") {
  stop("The first count-matrix column must be named gene_id")
}
if (anyDuplicated(counts_input$gene_id)) {
  stop("gene_id values must be unique")
}

sample_ids <- colnames(counts_input)[-1]
if (!length(sample_ids) || anyDuplicated(sample_ids)) {
  stop("Count-matrix sample names must be present and unique")
}
counts <- as.matrix(counts_input[, -1, drop = FALSE])
storage.mode(counts) <- "double"
rownames(counts) <- counts_input$gene_id
if (any(!is.finite(counts)) || any(counts < 0) || any(abs(counts - round(counts)) > 1e-8)) {
  stop("The count matrix must contain finite, non-negative integer counts")
}

metadata <- utils::read.delim(
  metadata_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
required_metadata <- c("sample_id", "state", "role")
missing_metadata <- setdiff(required_metadata, colnames(metadata))
if (length(missing_metadata)) {
  stop("Missing metadata columns: ", paste(missing_metadata, collapse = ", "))
}
if (anyDuplicated(metadata$sample_id)) {
  stop("metadata sample_id values must be unique")
}
if (!all(sample_ids %in% metadata$sample_id)) {
  stop("Metadata does not cover every count-matrix sample")
}
metadata <- metadata[match(sample_ids, metadata$sample_id), , drop = FALSE]
if (!identical(metadata$sample_id, sample_ids)) {
  stop("Count matrix and metadata could not be aligned")
}

dge_all <- edgeR::DGEList(counts = counts)
dge_all <- edgeR::calcNormFactors(dge_all)
logcpm_all <- edgeR::cpm(dge_all, log = TRUE, prior.count = 1)

train_mask <- metadata$role == "train" & metadata$state %in% c("naive", "primed")
if (sum(train_mask) < 10L || length(unique(metadata$state[train_mask])) != 2L) {
  stop("Training metadata must contain both naive and primed samples")
}

logcpm_train <- logcpm_all[, train_mask, drop = FALSE]
keep_genes <- rowSums(logcpm_train > 1) >= ncol(logcpm_train) * 0.20
if (!any(keep_genes)) {
  stop("No genes passed the training-sample expression filter")
}

logcpm_all_filtered <- logcpm_all[keep_genes, , drop = FALSE]
x_train <- t(logcpm_all_filtered[, train_mask, drop = FALSE])
y_train <- factor(metadata$state[train_mask], levels = c("primed", "naive"))

set.seed(20260128)
fit <- glmnet::cv.glmnet(
  x = x_train,
  y = y_train,
  family = "binomial",
  alpha = 0.5,
  nfolds = 10
)

x_all <- t(logcpm_all_filtered)
fitted_scores <- as.numeric(
  stats::predict(fit, newx = x_all, s = "lambda.min", type = "response")
)

set.seed(20260128)
n_train <- nrow(x_train)
fold_id <- sample(rep(seq_len(5L), length.out = n_train))
oof_scores <- rep(NA_real_, n_train)

for (fold in seq_len(5L)) {
  test_index <- which(fold_id == fold)
  train_index <- setdiff(seq_len(n_train), test_index)
  fold_fit <- glmnet::cv.glmnet(
    x = x_train[train_index, , drop = FALSE],
    y = y_train[train_index],
    family = "binomial",
    alpha = 0.5,
    nfolds = 5
  )
  oof_scores[test_index] <- as.numeric(
    stats::predict(
      fold_fit,
      newx = x_train[test_index, , drop = FALSE],
      s = "lambda.min",
      type = "response"
    )
  )
}

roc_object <- pROC::roc(
  response = y_train,
  predictor = oof_scores,
  levels = c("primed", "naive"),
  direction = "<",
  quiet = TRUE
)
oof_auc <- as.numeric(pROC::auc(roc_object))

model_bundle <- list(
  model = fit,
  features = rownames(logcpm_all_filtered),
  preprocessing = list(
    normalization = "TMM",
    expression_scale = "edgeR logCPM with prior.count = 1",
    filter = "logCPM > 1 in at least 20% of training samples",
    positive_class = "naive",
    alpha = 0.5,
    seed = 20260128
  )
)
saveRDS(model_bundle, file.path(output_dir, "naive_likeness_model.rds"))

utils::write.table(
  data.frame(
    sample_id = metadata$sample_id,
    state = metadata$state,
    role = metadata$role,
    naive_likeness_score = fitted_scores,
    check.names = FALSE
  ),
  file.path(output_dir, "naive_likeness_scores.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = ""
)
utils::write.table(
  data.frame(
    sample_id = metadata$sample_id[train_mask],
    state = metadata$state[train_mask],
    fold = fold_id,
    oof_probability = oof_scores,
    check.names = FALSE
  ),
  file.path(output_dir, "naive_likeness_oof_predictions.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = ""
)
utils::write.table(
  data.frame(feature = rownames(logcpm_all_filtered)),
  file.path(output_dir, "naive_likeness_features.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
utils::write.table(
  data.frame(
    training_samples = n_train,
    retained_features = ncol(x_train),
    lambda_min = fit$lambda.min,
    lambda_1se = fit$lambda.1se,
    sample_level_five_fold_oof_auc = oof_auc,
    check.names = FALSE
  ),
  file.path(output_dir, "naive_likeness_model_summary.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
