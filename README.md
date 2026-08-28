# Transcriptome-Based Multidimensional Scoring of Human iPSC States

This repository accompanies the manuscript **A Transcriptome-Based Multidimensional Scoring Framework for Quantitative Characterization of Human iPSC States**.

It contains the code, model objects, feature definitions, and saved evaluation outputs used to document the three transcriptome-based scores. The release package was assembled without model fitting, prediction, or reanalysis.

## Data availability

Public datasets are identified in Supplementary Table S1 of the manuscript submission package. The in-house data supporting the study are not included here and are available from the corresponding author upon reasonable request.

## Contents

### `pluripotency/`

- `code/score_pluripotency.R`: applies the final 276-feature LightGBM model to a prepared sample-by-feature matrix.
- `code/recompute_oof_metrics.R`: recomputes the six Figure 1b metrics from the saved five-model out-of-fold predictions. It does not fit any model.
- `data/features_276.txt`: ordered 276-gene input contract used by the current manuscript and Figure 1b.
- `model/pluripotency_lightgbm_276_features.rds`: final LightGBM model matching the 276-feature contract.
- `results/`: saved fold assignments, out-of-fold predictions, aggregate metrics, and metric ranks used for the corrected 734-sample Figure 1b comparison.
- `environment/software_versions.txt`: software versions recorded for the current Figure 1b evaluation.

### `esc_similarity/`

- `code/score_esc_similarity.R`: calculates the fitted intercept-plus-coefficient logit score from an input matrix already expressed on the model-compatible TMM-logCPM scale.
- `model/esc_similarity_glmnet_with_scale.rds`: ridge logistic-regression bundle containing the fitted model, 21,807 ordered features, and the corresponding center and scale vectors.

### `naive_likeness/`

- `code/train_and_evaluate_naive_likeness.R`: clean implementation of the reported TMM/logCPM preprocessing, training-sample expression filter, elastic-net model, fitted scores, and sample-level five-fold out-of-fold evaluation.

The clean training and evaluation code is included, but no fitted naive-likeness model is included in this local package. No model was refitted while assembling the package.

## Input conventions

- Tab-separated text files are used for expression matrices, metadata, scores, and metrics.
- Pluripotency scoring expects rows to be samples and columns to be the 276 ordered gene features, with `sample_id` as the first column.
- ESC similarity scoring expects rows to be genes and columns to be samples, with `gene_id` as the first column. Values must already be compatible with the TMM-logCPM scale used to fit the model.
- Naive-likeness training expects a raw-count matrix with genes in rows and samples in columns, plus metadata containing `sample_id`, `state`, and `role`.

## Example commands

```bash
Rscript pluripotency/code/recompute_oof_metrics.R \
  pluripotency/results/oof_predictions.tsv output_metrics

Rscript pluripotency/code/score_pluripotency.R \
  pluripotency/model/pluripotency_lightgbm_276_features.rds \
  pluripotency/data/features_276.txt prepared_pluripotency_matrix.tsv \
  pluripotency_scores.tsv

Rscript esc_similarity/code/score_esc_similarity.R \
  esc_similarity/model/esc_similarity_glmnet_with_scale.rds \
  prepared_esc_logcpm_matrix.tsv esc_similarity_scores.tsv

Rscript naive_likeness/code/train_and_evaluate_naive_likeness.R \
  naive_raw_counts.tsv naive_sample_metadata.tsv naive_output
```

## Software requirements

- Pluripotency scoring: R and `lightgbm`
- ESC similarity scoring: R and `glmnet`
- Naive-likeness training and evaluation: R, `edgeR`, `glmnet`, and `pROC`

The recorded environment for the Figure 1b evaluation is provided in `pluripotency/environment/software_versions.txt`.

## Citation

Please cite the associated article when using this code or the included model objects. The final publication citation will be added after publication.

## License

This repository is released under the MIT License. See `LICENSE`.
