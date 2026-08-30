# ESC similarity score

The model bundle contains the fitted ridge logistic-regression object, 21,807 ordered features, and the center and scale vectors estimated from the ESC/GTEx reference matrix.

`train_esc_similarity.R` reads raw ESC and GTEx count matrices, performs TMM normalization and logCPM transformation, adjusts source-dataset effects within the ESC reference samples, centers and scales the combined ESC/GTEx matrix, and fits the ridge logistic-regression model. The saved bundle contains the fitted model, ordered features, center vector, scale vector, and preprocessing description.

`score_esc_similarity.R` requires a gene-by-sample matrix on the same TMM-logCPM expression scale. It applies the saved center and scale vectors, then calculates the fitted intercept plus standardized expression multiplied by the coefficients selected at `lambda.min`. The output is the logit-scale ESC similarity score described in Methods 5.4, not a probability or an absolute measure of pluripotency.

All 21,807 model features are required. The script stops rather than silently dropping or imputing missing genes.

The input contract and a seven-dataset ESC manifest template are provided under `data/`.
