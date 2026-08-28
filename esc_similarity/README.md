# ESC similarity score

The model bundle contains the fitted ridge logistic-regression object, 21,807 ordered features, and the center and scale vectors estimated from the ESC/GTEx reference matrix.

`score_esc_similarity.R` requires a gene-by-sample matrix on the same TMM-logCPM expression scale. It applies the saved center and scale vectors, then calculates the fitted intercept plus standardized expression multiplied by the coefficients selected at `lambda.min`. The output is the logit-scale ESC similarity score described in Methods 5.4, not a probability or an absolute measure of pluripotency.

All 21,807 model features are required. The script stops rather than silently dropping or imputing missing genes.
