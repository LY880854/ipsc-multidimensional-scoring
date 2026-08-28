args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) {
  stop(
    "Usage: Rscript score_pluripotency.R MODEL_RDS FEATURES_TXT ",
    "PREPARED_SAMPLE_BY_FEATURE_TSV OUTPUT_TSV"
  )
}

model_file <- args[[1]]
feature_file <- args[[2]]
input_file <- args[[3]]
output_file <- args[[4]]

if (!requireNamespace("lightgbm", quietly = TRUE)) {
  stop("The lightgbm package is required")
}

features <- readLines(feature_file, warn = FALSE)
features <- features[nzchar(features)]
if (length(features) != 276L || anyDuplicated(features)) {
  stop("The feature file must contain 276 unique ordered genes")
}

input <- utils::read.delim(
  input_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
if (!ncol(input) || colnames(input)[[1]] != "sample_id") {
  stop("The first input column must be named sample_id")
}
if (anyDuplicated(input$sample_id)) {
  stop("sample_id values must be unique")
}
if (!identical(colnames(input)[-1], features)) {
  stop("Input feature names and order must exactly match features_276.txt")
}

x <- as.matrix(input[, -1, drop = FALSE])
storage.mode(x) <- "double"
if (any(!is.finite(x))) {
  stop("The prepared expression matrix contains non-finite values")
}

model <- readRDS(model_file)
scores <- as.numeric(stats::predict(model, x))
if (length(scores) != nrow(input) || any(!is.finite(scores))) {
  stop("Model prediction failed")
}

output <- data.frame(
  sample_id = input$sample_id,
  pluripotency_score = scores,
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
