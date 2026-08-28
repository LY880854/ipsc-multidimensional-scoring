# Pluripotency score

The current manuscript uses an ordered 276-gene input contract. The included RDS is the final LightGBM model associated with that contract.

`score_pluripotency.R` applies the model to an already prepared sample-by-feature matrix. It does not perform identifier conversion, normalization, batch adjustment, missing-feature imputation, or model fitting. The input must therefore match the preprocessing and feature order described in Methods 5.2 and 5.3.

`recompute_oof_metrics.R` independently recomputes MAE, MSE, RMSE, R-squared, Pearson correlation, and Spearman correlation from the saved 734-sample out-of-fold predictions used for Figure 1b. It does not retrain the five compared models.
