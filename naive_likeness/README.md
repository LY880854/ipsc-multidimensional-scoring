# Naive likeness score

The included script is a clean implementation of the model construction and sample-level out-of-fold evaluation reported in Methods 5.5. It performs TMM normalization and logCPM transformation, retains genes with logCPM greater than 1 in at least 20% of training samples, fits an elastic-net binomial model with `alpha = 0.5`, and calculates naive-class probabilities.

The final model uses 10-fold internal cross-validation for `lambda.min`. The reported internal ROC procedure is reproduced by an outer sample-level five-fold split, with a five-fold `cv.glmnet` fit inside each training partition.

No fitted naive-likeness model is included. The script is provided for transparent model construction and later controlled export from the existing training inputs; it was not executed while assembling this package.
