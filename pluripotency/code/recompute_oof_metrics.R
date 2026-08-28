args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: Rscript recompute_oof_metrics.R OOF_PREDICTIONS_TSV OUTPUT_DIR")
}

input_file <- args[[1]]
output_dir <- args[[2]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

options(stringsAsFactors = FALSE, digits = 17)

pred <- utils::read.delim(
  input_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

required_columns <- c(
  "order", "sample_id", "group", "fold", "model",
  "observed_label", "oof_prediction", "residual"
)
missing_columns <- setdiff(required_columns, colnames(pred))
if (length(missing_columns)) {
  stop("Missing columns: ", paste(missing_columns, collapse = ", "))
}

if (any(!is.finite(pred$observed_label)) || any(!is.finite(pred$oof_prediction))) {
  stop("Observed labels and predictions must be finite")
}

model_order <- c("Elastic Net", "LightGBM", "Random Forest", "SVM", "XGBoost")
if (!setequal(unique(pred$model), model_order)) {
  stop("Unexpected model names in the prediction table")
}

counts <- table(factor(pred$model, levels = model_order))
if (any(counts != 734L)) {
  stop("Each model must contain 734 out-of-fold predictions")
}

for (model_name in model_order) {
  x <- pred[pred$model == model_name, , drop = FALSE]
  if (anyDuplicated(x$sample_id) || !identical(sort(x$order), seq_len(734L))) {
    stop("Sample identity or order failure for ", model_name)
  }
}

metric_values <- function(observed, predicted) {
  denominator <- sum((observed - mean(observed))^2)
  if (!is.finite(denominator) || denominator <= 0) {
    stop("R-squared denominator is not positive")
  }
  c(
    MAE = mean(abs(observed - predicted)),
    MSE = mean((predicted - observed)^2),
    RMSE = sqrt(mean((predicted - observed)^2)),
    R_squared = 1 - sum((predicted - observed)^2) / denominator,
    Pearson = as.numeric(stats::cor(predicted, observed, method = "pearson")),
    Spearman = as.numeric(stats::cor(predicted, observed, method = "spearman"))
  )
}

aggregate_rows <- lapply(model_order, function(model_name) {
  x <- pred[pred$model == model_name, , drop = FALSE]
  values <- metric_values(x$observed_label, x$oof_prediction)
  data.frame(
    model = model_name,
    n = nrow(x),
    as.list(values),
    check.names = FALSE
  )
})
aggregate <- do.call(rbind, aggregate_rows)

metric_order <- c("MAE", "MSE", "RMSE", "R_squared", "Pearson", "Spearman")
direction <- c(
  MAE = "lower_is_better",
  MSE = "lower_is_better",
  RMSE = "lower_is_better",
  R_squared = "higher_is_better",
  Pearson = "higher_is_better",
  Spearman = "higher_is_better"
)

rank_rows <- lapply(metric_order, function(metric_name) {
  values <- aggregate[[metric_name]]
  rank_input <- if (direction[[metric_name]] == "lower_is_better") values else -values
  data.frame(
    metric = metric_name,
    direction = direction[[metric_name]],
    model = aggregate$model,
    value = values,
    rank = rank(rank_input, ties.method = "min"),
    check.names = FALSE
  )
})
ranks <- do.call(rbind, rank_rows)

utils::write.table(
  aggregate,
  file.path(output_dir, "model_metrics_aggregate.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = ""
)
utils::write.table(
  ranks,
  file.path(output_dir, "model_metric_ranks.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = ""
)

message("Metrics written to ", normalizePath(output_dir))
