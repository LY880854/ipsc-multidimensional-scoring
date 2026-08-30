# Transcriptome-Based Multidimensional Scoring of Human iPSC States

This repository accompanies the manuscript **A Transcriptome-Based Multidimensional Scoring Framework for Quantitative Characterization of Human iPSC States**.

It contains the code, model objects, feature definitions, and saved evaluation outputs for the three transcriptome-based scores.

## Data availability

Public datasets are identified in Supplementary Table S1 of the manuscript. In-house sequencing data are not distributed through this code repository.

## Contents

### `pluripotency/`

- `code/score_pluripotency.R`: applies the final 276-feature LightGBM model to a prepared sample-by-feature matrix.
- `code/prepare_pluripotency_inputs.R`: combines the three expression sources, applies the reported log transformation and batch adjustment, and writes the ordered 276-feature model input.
- `code/train_and_compare_models.R`: trains the five compared algorithms under common sample-level five-fold assignments and writes their out-of-fold predictions and metrics.
- `code/recompute_oof_metrics.R`: recomputes the six Figure 1b metrics from the saved five-model out-of-fold predictions. It does not fit any model.
- `data/features_276.txt`: ordered 276-gene input contract used by the current manuscript and Figure 1b.
- `data/input_schema.md`: input columns and matrix orientation for preparation and model comparison.
- `model/pluripotency_lightgbm_276_features.rds`: final LightGBM model matching the 276-feature contract.
- `results/`: saved fold assignments, out-of-fold predictions, aggregate metrics, and metric ranks used for the corrected 734-sample Figure 1b comparison.
- `environment/software_versions.txt`: software versions recorded for the current Figure 1b evaluation.

### `esc_similarity/`

- `code/score_esc_similarity.R`: calculates the fitted intercept-plus-coefficient logit score from an input matrix already expressed on the model-compatible TMM-logCPM scale.
- `code/train_esc_similarity.R`: performs the reported ESC/GTEx preprocessing and fits the ridge logistic-regression model.
- `data/`: model input schema and an ESC count-file manifest template.
- `model/esc_similarity_glmnet_with_scale.rds`: ridge logistic-regression bundle containing the fitted model, 21,807 ordered features, and the corresponding center and scale vectors.

### `naive_likeness/`

- `code/assemble_naive_likeness_inputs.R`: assembles the public raw-count inputs using a count-file manifest and sample metadata.
- `code/train_and_evaluate_naive_likeness.R`: clean implementation of the reported TMM/logCPM preprocessing, training-sample expression filter, elastic-net model, fitted scores, and sample-level five-fold out-of-fold evaluation.
- `code/score_naive_likeness.R`: applies a fitted naïve-likeness model bundle to new raw-count data.
- `data/`: input schema and a count-file manifest template.

The training script writes a fitted naïve-likeness model bundle. A pre-fitted naïve-likeness model object is not included.

## Input conventions

- Tab-separated text files are used for expression matrices, metadata, scores, and metrics.
- Pluripotency scoring expects rows to be samples and columns to be the 276 ordered gene features, with `sample_id` as the first column.
- ESC similarity scoring expects rows to be genes and columns to be samples, with `gene_id` as the first column. Values must already be compatible with the TMM-logCPM scale used to fit the model.
- Naive-likeness training expects a raw-count matrix with genes in rows and samples in columns, plus metadata containing `sample_id`, `state`, and `role`.
- File-level schemas and manifest templates are provided in each model directory under `data/`.

## Example commands

```bash
Rscript pluripotency/code/recompute_oof_metrics.R \
  pluripotency/results/oof_predictions.tsv output_metrics

Rscript pluripotency/code/prepare_pluripotency_inputs.R \
  gtex_tpm.csv ipsc_tpm.csv asc_tpm.csv asc_metadata.csv \
  pluripotency/results/fold_assignment.tsv \
  pluripotency/data/features_276.txt pluripotency_prepared

Rscript pluripotency/code/train_and_compare_models.R \
  pluripotency_prepared/pluripotency_prepared_expression.tsv \
  pluripotency_prepared/pluripotency_sample_metadata.tsv \
  pluripotency/data/features_276.txt pluripotency_oof

Rscript pluripotency/code/score_pluripotency.R \
  pluripotency/model/pluripotency_lightgbm_276_features.rds \
  pluripotency/data/features_276.txt prepared_pluripotency_matrix.tsv \
  pluripotency_scores.tsv

Rscript esc_similarity/code/score_esc_similarity.R \
  esc_similarity/model/esc_similarity_glmnet_with_scale.rds \
  prepared_esc_logcpm_matrix.tsv esc_similarity_scores.tsv

Rscript esc_similarity/code/train_esc_similarity.R \
  gtex_raw_counts.csv esc_count_manifest.tsv esc_training_output

Rscript naive_likeness/code/assemble_naive_likeness_inputs.R \
  naive_count_manifest.tsv naive_sample_metadata.tsv \
  naive_raw_counts.tsv naive_aligned_metadata.tsv

Rscript naive_likeness/code/train_and_evaluate_naive_likeness.R \
  naive_raw_counts.tsv naive_aligned_metadata.tsv naive_output

Rscript naive_likeness/code/score_naive_likeness.R \
  naive_output/naive_likeness_model.rds new_raw_counts.tsv \
  naive_likeness_scores.tsv
```

## Software requirements

- Pluripotency preparation, training, and scoring: R, `sva`, `caret`, `e1071`, `glmnet`, `lightgbm`, `randomForest`, and `xgboost`
- ESC similarity training and scoring: R, `edgeR`, `glmnet`, and `limma`
- Naive-likeness training, evaluation, and scoring: R, `edgeR`, `glmnet`, and `pROC`

The recorded environment for the Figure 1b evaluation is provided in `pluripotency/environment/software_versions.txt`.

## Citation

Please cite the associated article when using this code or the included model objects. The final publication citation will be added after publication.

## License

This repository is released under the MIT License. See `LICENSE`.
