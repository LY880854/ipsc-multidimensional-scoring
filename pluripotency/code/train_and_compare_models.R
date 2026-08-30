args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) {
  stop(
    "Usage: Rscript train_and_compare_models.R ",
    "PREPARED_EXPRESSION_TSV SAMPLE_METADATA_TSV FEATURES_TXT OUTPUT_DIR"
  )
}

expression_file <- args[[1]]
metadata_file <- args[[2]]
feature_file <- args[[3]]
output_dir <- args[[4]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

required_packages <- c(
  "caret", "e1071", "glmnet", "lightgbm", "randomForest", "xgboost"
)
available <- vapply(required_packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))
if (!all(available)) {
  stop("Missing packages: ", paste(names(available)[!available], collapse = ", "))
}

metric_values <- function(observed, predicted) {
  denominator <- sum((observed - mean(observed))^2)
  if (!is.finite(denominator) || denominator <= 0) {
    stop("R-squared is undefined because the observed labels have zero variance")
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

features <- readLines(feature_file, warn = FALSE)
features <- features[nzchar(features)]
if (!length(features) || anyDuplicated(features)) {
  stop("The feature list must contain unique, non-empty gene identifiers")
}

expression_input <- utils::read.delim(
  expression_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
if (!ncol(expression_input) || colnames(expression_input)[[1]] != "sample_id") {
  stop("The first expression-matrix column must be named sample_id")
}
if (anyDuplicated(expression_input$sample_id)) {
  stop("Expression-matrix sample_id values must be unique")
}
missing_features <- setdiff(features, colnames(expression_input))
if (length(missing_features)) {
  stop("Missing required features: ", paste(missing_features, collapse = ", "))
}
x_raw <- as.matrix(expression_input[, features, drop = FALSE])
storage.mode(x_raw) <- "double"
rownames(x_raw) <- expression_input$sample_id
if (any(!is.finite(x_raw))) {
  stop("The prepared expression matrix contains non-finite values")
}

metadata <- utils::read.delim(
  metadata_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
required_metadata <- c("sample_id", "observed_label")
missing_metadata <- setdiff(required_metadata, colnames(metadata))
if (length(missing_metadata)) {
  stop("Missing metadata columns: ", paste(missing_metadata, collapse = ", "))
}
if (anyDuplicated(metadata$sample_id) || !all(rownames(x_raw) %in% metadata$sample_id)) {
  stop("Metadata sample_id values must uniquely cover the expression matrix")
}
metadata <- metadata[match(rownames(x_raw), metadata$sample_id), , drop = FALSE]
y <- as.numeric(metadata$observed_label)
if (any(!is.finite(y))) {
  stop("observed_label must contain finite numeric values")
}

if ("fold" %in% colnames(metadata)) {
  fold_id <- as.integer(metadata$fold)
  fold_levels <- sort(unique(fold_id))
  if (!identical(fold_levels, seq_len(5L))) {
    stop("A supplied fold column must use the integers 1 through 5")
  }
  folds <- lapply(fold_levels, function(fold) which(fold_id == fold))
} else {
  set.seed(123)
  folds <- caret::createFolds(y, k = 5L, list = TRUE)
  fold_id <- integer(length(y))
  for (fold in seq_along(folds)) {
    fold_id[folds[[fold]]] <- fold
  }
}

model_names <- c("Elastic Net", "LightGBM", "Random Forest", "SVM", "XGBoost")
oof <- matrix(
  NA_real_,
  nrow = nrow(x_raw),
  ncol = length(model_names),
  dimnames = list(rownames(x_raw), model_names)
)

for (fold in seq_along(folds)) {
  validation_index <- as.integer(folds[[fold]])
  training_index <- setdiff(seq_len(nrow(x_raw)), validation_index)
  x_training_raw <- x_raw[training_index, , drop = FALSE]
  x_validation_raw <- x_raw[validation_index, , drop = FALSE]
  y_training <- y[training_index]

  training_center <- colMeans(x_training_raw)
  training_scale <- apply(x_training_raw, 2L, stats::sd)
  active <- is.finite(training_scale) & training_scale > 0
  x_training <- matrix(0, nrow(x_training_raw), ncol(x_training_raw), dimnames = dimnames(x_training_raw))
  x_validation <- matrix(0, nrow(x_validation_raw), ncol(x_validation_raw), dimnames = dimnames(x_validation_raw))
  x_training[, active] <- sweep(
    sweep(x_training_raw[, active, drop = FALSE], 2L, training_center[active], "-"),
    2L,
    training_scale[active],
    "/"
  )
  x_validation[, active] <- sweep(
    sweep(x_validation_raw[, active, drop = FALSE], 2L, training_center[active], "-"),
    2L,
    training_scale[active],
    "/"
  )

  set.seed(123)
  lightgbm_data <- lightgbm::lgb.Dataset(
    data = x_training,
    label = y_training,
    free_raw_data = FALSE
  )
  lightgbm_model <- lightgbm::lgb.train(
    params = list(
      objective = "regression",
      metric = "rmse",
      learning_rate = 0.05,
      num_leaves = 31L,
      seed = 123L,
      bagging_seed = 123L,
      feature_fraction_seed = 123L,
      data_random_seed = 123L,
      deterministic = TRUE,
      force_col_wise = TRUE,
      verbosity = -1L
    ),
    data = lightgbm_data,
    nrounds = 200L,
    verbose = -1L
  )
  oof[validation_index, "LightGBM"] <- as.numeric(
    stats::predict(lightgbm_model, x_validation)
  )

  set.seed(123)
  xgboost_model <- xgboost::xgb.train(
    params = list(
      objective = "reg:squarederror",
      eta = 0.05,
      max_depth = 6L,
      seed = 123L
    ),
    data = xgboost::xgb.DMatrix(x_training, label = y_training),
    nrounds = 200L,
    verbose = 0L
  )
  oof[validation_index, "XGBoost"] <- as.numeric(
    stats::predict(xgboost_model, xgboost::xgb.DMatrix(x_validation))
  )

  set.seed(123)
  random_forest_model <- randomForest::randomForest(
    x = x_training,
    y = y_training,
    ntree = 200L
  )
  oof[validation_index, "Random Forest"] <- as.numeric(
    stats::predict(random_forest_model, x_validation)
  )

  set.seed(123)
  elastic_net_model <- glmnet::cv.glmnet(
    x = x_training,
    y = y_training,
    alpha = 0.5,
    standardize = FALSE
  )
  oof[validation_index, "Elastic Net"] <- as.numeric(
    stats::predict(
      elastic_net_model,
      newx = x_validation,
      s = "lambda.min"
    )
  )

  set.seed(123)
  svm_model <- e1071::svm(
    x = x_training,
    y = y_training,
    kernel = "linear",
    scale = FALSE
  )
  oof[validation_index, "SVM"] <- as.numeric(
    stats::predict(svm_model, x_validation)
  )
}

if (any(!is.finite(oof))) {
  stop("One or more out-of-fold predictions are missing or non-finite")
}

oof_rows <- lapply(model_names, function(model_name) {
  data.frame(
    sample_id = rownames(x_raw),
    group = if ("group" %in% colnames(metadata)) metadata$group else "",
    fold = fold_id,
    model = model_name,
    observed_label = y,
    oof_prediction = oof[, model_name],
    residual = y - oof[, model_name],
    check.names = FALSE
  )
})
oof_table <- do.call(rbind, oof_rows)

aggregate_rows <- lapply(model_names, function(model_name) {
  values <- metric_values(y, oof[, model_name])
  data.frame(
    model = model_name,
    n = length(y),
    as.list(values),
    check.names = FALSE
  )
})
aggregate_metrics <- do.call(rbind, aggregate_rows)

fold_rows <- list()
row_index <- 1L
for (model_name in model_names) {
  for (fold in seq_along(folds)) {
    index <- which(fold_id == fold)
    values <- metric_values(y[index], oof[index, model_name])
    fold_rows[[row_index]] <- data.frame(
      model = model_name,
      fold = fold,
      n = length(index),
      as.list(values),
      check.names = FALSE
    )
    row_index <- row_index + 1L
  }
}
fold_metrics <- do.call(rbind, fold_rows)

utils::write.table(
  data.frame(
    sample_id = rownames(x_raw),
    group = if ("group" %in% colnames(metadata)) metadata$group else "",
    observed_label = y,
    fold = fold_id,
    check.names = FALSE
  ),
  file.path(output_dir, "fold_assignment.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
utils::write.table(
  oof_table,
  file.path(output_dir, "oof_predictions.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
utils::write.table(
  aggregate_metrics,
  file.path(output_dir, "model_metrics_aggregate.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
utils::write.table(
  fold_metrics,
  file.path(output_dir, "model_metrics_by_fold.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
