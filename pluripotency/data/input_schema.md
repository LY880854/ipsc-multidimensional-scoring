# Pluripotency input schema

`prepare_pluripotency_inputs.R` accepts three non-negative TPM matrices in CSV format. The first column contains unique gene identifiers and the remaining columns are samples. The ASC metadata CSV contains `sample`, `cell_type`, and `source`.

The label TSV contains one row per selected sample, in the intended output order. It requires `sample_id` and `observed_label`; `group` and `fold` are retained when present. The supplied `results/fold_assignment.tsv` can be used as this table.

`train_and_compare_models.R` consumes the two files produced by the preparation script. The expression TSV has `sample_id` followed by the 276 columns listed in `features_276.txt`. The metadata TSV has `sample_id`, `group`, and `observed_label`, with an optional five-fold assignment in `fold`.
