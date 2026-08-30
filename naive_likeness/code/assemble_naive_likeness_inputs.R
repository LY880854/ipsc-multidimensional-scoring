args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) {
  stop(
    "Usage: Rscript assemble_naive_likeness_inputs.R ",
    "COUNT_MANIFEST_TSV SAMPLE_METADATA_TSV OUTPUT_COUNTS_TSV OUTPUT_METADATA_TSV"
  )
}

manifest_file <- args[[1]]
metadata_file <- args[[2]]
output_counts_file <- args[[3]]
output_metadata_file <- args[[4]]

read_counts <- function(path) {
  extension <- tolower(tools::file_ext(path))
  input <- if (extension == "tsv" || extension == "txt") {
    utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
  } else {
    utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  }
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
  stop("The manifest must contain unique, non-empty dataset names")
}
if (any(!file.exists(manifest$path))) {
  stop("Missing count files: ", paste(manifest$path[!file.exists(manifest$path)], collapse = ", "))
}

metadata <- utils::read.delim(
  metadata_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
required_metadata <- c("sample_id", "dataset", "state", "role")
missing_metadata <- setdiff(required_metadata, colnames(metadata))
if (length(missing_metadata)) {
  stop("Missing metadata columns: ", paste(missing_metadata, collapse = ", "))
}
if (anyDuplicated(metadata$sample_id) || any(metadata$sample_id == "")) {
  stop("Metadata sample_id values must be present and unique")
}
if (!all(metadata$dataset %in% manifest$dataset)) {
  stop("Every metadata dataset must have a row in the count manifest")
}

matrices <- setNames(lapply(manifest$path, read_counts), manifest$dataset)
common_genes <- Reduce(intersect, lapply(matrices, rownames))
if (!length(common_genes)) {
  stop("No gene identifiers are shared by all count matrices")
}

selected_matrices <- lapply(names(matrices), function(dataset) {
  expected_samples <- metadata$sample_id[metadata$dataset == dataset]
  if (!length(expected_samples)) {
    return(NULL)
  }
  missing_samples <- setdiff(expected_samples, colnames(matrices[[dataset]]))
  if (length(missing_samples)) {
    stop(
      "Samples listed in metadata are missing from ", dataset, ": ",
      paste(missing_samples, collapse = ", ")
    )
  }
  matrices[[dataset]][common_genes, expected_samples, drop = FALSE]
})
selected_matrices <- selected_matrices[!vapply(selected_matrices, is.null, logical(1))]
combined <- do.call(cbind, selected_matrices)
if (anyDuplicated(colnames(combined)) || !setequal(colnames(combined), metadata$sample_id)) {
  stop("Assembled sample identifiers do not uniquely match the metadata")
}
combined <- combined[, metadata$sample_id, drop = FALSE]

dir.create(dirname(output_counts_file), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(output_metadata_file), recursive = TRUE, showWarnings = FALSE)
utils::write.table(
  data.frame(gene_id = rownames(combined), combined, check.names = FALSE),
  output_counts_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
utils::write.table(
  metadata,
  output_metadata_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
