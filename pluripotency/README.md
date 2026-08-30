# Pluripotency score

The current manuscript uses an ordered 276-gene input contract. The included RDS is the final LightGBM model associated with that contract.

`prepare_pluripotency_inputs.R` combines the GTEx, healthy-iPSC, and adult-stem-cell TPM matrices, applies `log2(TPM + 1)` and ComBat adjustment by source, and writes the sample-by-feature matrix in the order defined by the label table and `features_276.txt`.

`train_and_compare_models.R` fits elastic net, LightGBM, random forest, linear SVM, and XGBoost using common sample-level five-fold assignments. Centering and scaling are learned within each training fold. The script writes fold assignments, sample-level out-of-fold predictions, aggregate metrics, and metrics by fold.

`score_pluripotency.R` applies the model to an already prepared sample-by-feature matrix. It does not perform identifier conversion, normalization, batch adjustment, missing-feature imputation, or model fitting. The input must therefore match the preprocessing and feature order described in Methods 5.2 and 5.3.

`recompute_oof_metrics.R` independently recomputes MAE, MSE, RMSE, R-squared, Pearson correlation, and Spearman correlation from the saved 734-sample out-of-fold predictions used for Figure 1b. It does not retrain the five compared models.

The required columns and matrix orientations are listed in `data/input_schema.md`.
