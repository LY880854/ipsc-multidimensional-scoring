# ESC similarity input schema

`train_esc_similarity.R` accepts one GTEx raw-count CSV and an ESC manifest TSV. Each count file has unique gene identifiers in the first column and unique sample identifiers in the remaining columns. The manifest contains `dataset` and `path`, with one row per ESC source dataset; see `esc_training_manifest_template.tsv`.

The training script applies TMM normalization and logCPM transformation, adjusts ESC source-dataset effects, fits the ridge logistic-regression model, and saves the model, ordered features, center vector, and scale vector as one RDS bundle.

`score_esc_similarity.R` accepts a gene-by-sample TSV on the model-compatible TMM-logCPM scale. Its first column is `gene_id`, and all model features must be present.
