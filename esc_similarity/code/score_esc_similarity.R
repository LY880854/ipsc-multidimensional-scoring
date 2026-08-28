args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop(
    "Usage: Rscript score_esc_similarity.R MODEL_BUNDLE_RDS ",
    "GENE_BY_SAMPLE_LOGCPM_TSV OUTPUT_TSV"
  )
}

model_file <- args[[1]]
input_file <- args[[2]]
output_file <- args[[3]]

if (!requireNamespace("glmnet", quietly = TRUE)) {
  stop("The glmnet package is required")
}

bundle <- readRDS(model_file)
required_objects <- c("model", "center", "scale", "features")
if (!all(required_objects %in% names(bundle))) {
  stop("The model bundle must contain model, center, scale, and features")
}

features <- as.character(bundle$features)
if (length(features) != 21807L || anyDuplicated(features)) {
  stop("The ESC model must contain 21,807 unique ordered features")
}
if (length(bundle$center) != length(features) || length(bundle$scale) != length(features)) {
  stop("Center and scale vectors do not match the model features")
}
if (any(!is.finite(bundle$center)) || any(!is.finite(bundle$scale)) || any(bundle$scale == 0)) {
  stop("Center and scale vectors are not valid")
}

input <- utils::read.delim(
  input_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
if (!ncol(input) || colnames(input)[[1]] != "gene_id") {
  stop("The first input column must be named gene_id")
}
if (anyDuplicated(input$gene_id)) {
  stop("gene_id values must be unique")
}

missing_features <- setdiff(features, input$gene_id)
if (length(missing_features)) {
  stop(
    "Input matrix is missing ", length(missing_features),
    " required model features; no imputation was performed"
  )
}

expression <- as.matrix(input[match(features, input$gene_id), -1, drop = FALSE])
storage.mode(expression) <- "double"
if (any(!is.finite(expression))) {
  stop("The expression matrix contains non-finite values")
}

x <- t(expression)
center <- as.numeric(bundle$center)
scale <- as.numeric(bundle$scale)
names(center) <- features
names(scale) <- features
x_scaled <- sweep(x, 2, center, FUN = "-")
x_scaled <- sweep(x_scaled, 2, scale, FUN = "/")

coefficient_matrix <- as.matrix(stats::coef(bundle$model, s = "lambda.min"))
if (!"(Intercept)" %in% rownames(coefficient_matrix)) {
  stop("The fitted intercept is missing")
}
if (!all(features %in% rownames(coefficient_matrix))) {
  stop("Model coefficients do not cover all ordered features")
}
intercept <- coefficient_matrix["(Intercept)", 1]
weights <- coefficient_matrix[features, 1]
scores <- as.numeric(intercept + x_scaled %*% weights)

output <- data.frame(
  sample_id = rownames(x_scaled),
  esc_similarity_score = scores,
  check.names = FALSE
)
utils::write.table(
  output,
  output_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = ""
)
