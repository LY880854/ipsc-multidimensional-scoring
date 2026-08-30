args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 7L) {
  stop(
    "Usage: Rscript prepare_pluripotency_inputs.R ",
    "GTEX_TPM_CSV IPSC_TPM_CSV ASC_TPM_CSV ASC_METADATA_CSV ",
    "LABEL_TSV FEATURES_TXT OUTPUT_DIR"
  )
}

gtex_file <- args[[1]]
ipsc_file <- args[[2]]
asc_file <- args[[3]]
asc_metadata_file <- args[[4]]
label_file <- args[[5]]
feature_file <- args[[6]]
output_dir <- args[[7]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

if (!requireNamespace("sva", quietly = TRUE)) {
  stop("Missing package: sva")
}

read_expression <- function(path) {
  input <- utils::read.csv(
    path,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  if (ncol(input) < 2L) {
    stop("Expression files must contain a gene column and at least one sample: ", path)
  }
  gene_id <- as.character(input[[1]])
  if (anyNA(gene_id) || any(gene_id == "") || anyDuplicated(gene_id)) {
    stop("Gene identifiers must be present and unique: ", path)
  }
  values <- as.matrix(input[, -1, drop = FALSE])
  storage.mode(values) <- "double"
  rownames(values) <- gene_id
  if (any(!is.finite(values)) || any(values < 0)) {
    stop("Expression values must be finite and non-negative: ", path)
  }
  if (!ncol(values) || anyDuplicated(colnames(values))) {
    stop("Sample identifiers must be present and unique: ", path)
  }
  values
}

complete_genes <- function(expression, target_genes) {
  missing_genes <- setdiff(target_genes, rownames(expression))
  if (length(missing_genes)) {
    zeros <- matrix(
      0,
      nrow = length(missing_genes),
      ncol = ncol(expression),
      dimnames = list(missing_genes, colnames(expression))
    )
    expression <- rbind(expression, zeros)
  }
  expression[target_genes, , drop = FALSE]
}

features <- readLines(feature_file, warn = FALSE)
features <- features[nzchar(features)]
if (!length(features) || anyDuplicated(features)) {
  stop("The feature list must contain unique, non-empty gene identifiers")
}

gtex <- read_expression(gtex_file)
ipsc <- read_expression(ipsc_file)
asc <- read_expression(asc_file)

asc_metadata <- utils::read.csv(
  asc_metadata_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
required_asc_columns <- c("sample", "cell_type", "source")
missing_asc_columns <- setdiff(required_asc_columns, colnames(asc_metadata))
if (length(missing_asc_columns)) {
  stop("Missing ASC metadata columns: ", paste(missing_asc_columns, collapse = ", "))
}
if (!identical(colnames(asc), asc_metadata$sample)) {
  stop("ASC expression columns and metadata rows must be in the same order")
}

all_samples <- c(colnames(gtex), colnames(ipsc), colnames(asc))
if (anyDuplicated(all_samples)) {
  stop("Sample identifiers must be unique across the three expression files")
}

all_genes <- unique(c(rownames(gtex), rownames(ipsc), rownames(asc), features))
gtex <- complete_genes(gtex, all_genes)
ipsc <- complete_genes(ipsc, all_genes)
asc <- complete_genes(asc, all_genes)
combined <- cbind(gtex, ipsc, asc)

metadata <- rbind(
  data.frame(sample_id = colnames(gtex), group = "GTEx", source = "GTEx"),
  data.frame(sample_id = colnames(ipsc), group = "iPSC", source = "iPSC"),
  data.frame(
    sample_id = asc_metadata$sample,
    group = asc_metadata$cell_type,
    source = asc_metadata$source
  )
)
metadata <- metadata[match(colnames(combined), metadata$sample_id), , drop = FALSE]
if (!identical(metadata$sample_id, colnames(combined))) {
  stop("Combined expression matrix and metadata could not be aligned")
}

labels <- utils::read.delim(
  label_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
required_label_columns <- c("sample_id", "observed_label")
missing_label_columns <- setdiff(required_label_columns, colnames(labels))
if (length(missing_label_columns)) {
  stop("Missing label columns: ", paste(missing_label_columns, collapse = ", "))
}
if (anyDuplicated(labels$sample_id) || !all(labels$sample_id %in% metadata$sample_id)) {
  stop("LABEL_TSV sample identifiers must be unique and present in the expression data")
}

log_expression <- log2(combined + 1)
adjusted <- sva::ComBat(
  dat = log_expression,
  batch = metadata$source,
  par.prior = TRUE
)

sample_order <- labels$sample_id
prepared <- t(adjusted[features, sample_order, drop = FALSE])
if (any(!is.finite(prepared))) {
  stop("Prepared expression matrix contains non-finite values")
}

prepared_table <- data.frame(
  sample_id = rownames(prepared),
  prepared,
  check.names = FALSE
)
output_metadata <- data.frame(
  sample_id = sample_order,
  group = if ("group" %in% colnames(labels)) labels$group else metadata$group[match(sample_order, metadata$sample_id)],
  observed_label = as.numeric(labels$observed_label),
  check.names = FALSE
)
if ("fold" %in% colnames(labels)) {
  output_metadata$fold <- labels$fold
}
if (any(!is.finite(output_metadata$observed_label))) {
  stop("Observed labels must be finite numeric values")
}

utils::write.table(
  prepared_table,
  file.path(output_dir, "pluripotency_prepared_expression.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
utils::write.table(
  output_metadata,
  file.path(output_dir, "pluripotency_sample_metadata.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
