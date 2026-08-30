# Naive-likeness input schema

`assemble_naive_likeness_inputs.R` accepts a count manifest TSV with `dataset` and `path`, plus a sample metadata TSV with `sample_id`, `dataset`, `state`, and `role`. Count files may be CSV or TSV; the first column contains unique gene identifiers and the remaining columns are samples. A manifest template is included.

The assembly script retains genes shared by all listed count matrices, selects samples named in the metadata, and writes the combined raw-count matrix and aligned metadata used by `train_and_evaluate_naive_likeness.R`.

The training script requires `state` values `naive` and `primed` for rows whose `role` is `train`. It saves a fitted model bundle and the internal out-of-fold outputs. `score_naive_likeness.R` applies a saved bundle to a new raw-count TSV with `gene_id` as the first column.
