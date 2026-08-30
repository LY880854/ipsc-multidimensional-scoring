args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop(
    "Usage: Rscript score_naive_likeness.R ",
    "MODEL_RDS RAW_COUNT_TSV OUTPUT_SCORES_TSV"
  )
}

model_file <- args[[1]]
count_file <- args[[2]]
output_file <- args[[3]]

required_packages <- c("edgeR", "glmnet")
available <- vapply(required_packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))
if (!all(available)) {
  stop("Missing packages: ", paste(names(available)[!available], collapse = ", "))
}

bundle <- readRDS(model_file)
required_bundle <- c("model", "features")
missing_bundle <- setdiff(required_bundle, names(bundle))
if (length(missing_bundle)) {
  stop("Model bundle is missing: ", paste(missing_bundle, collapse = ", "))
}
if (!length(bundle$features) || anyDuplicated(bundle$features)) {
  stop("Model features must be present and unique")
}

input <- utils::read.delim(
  count_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
if (!ncol(input) || colnames(input)[[1]] != "gene_id") {
  stop("The first count-matrix column must be named gene_id")
}
if (anyDuplicated(input$gene_id)) {
  stop("gene_id values must be unique")
}
counts <- as.matrix(input[, -1, drop = FALSE])
storage.mode(counts) <- "double"
rownames(counts) <- as.character(input$gene_id)
if (any(!is.finite(counts)) || any(counts < 0) || any(abs(counts - round(counts)) > 1e-8)) {
  stop("The count matrix must contain finite, non-negative integer counts")
}
if (!ncol(counts) || anyDuplicated(colnames(counts))) {
  stop("Sample identifiers must be present and unique")
}
missing_features <- setdiff(bundle$features, rownames(counts))
if (length(missing_features)) {
  stop("Missing model features: ", paste(missing_features, collapse = ", "))
}

dge <- edgeR::DGEList(counts = counts)
dge <- edgeR::calcNormFactors(dge, method = "TMM")
logcpm <- edgeR::cpm(dge, log = TRUE, prior.count = 1)
x <- t(logcpm[bundle$features, , drop = FALSE])
scores <- as.numeric(
  stats::predict(bundle$model, newx = x, s = "lambda.min", type = "response")
)

dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
utils::write.table(
  data.frame(
    sample_id = rownames(x),
    naive_likeness_score = scores,
    check.names = FALSE
  ),
  output_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
