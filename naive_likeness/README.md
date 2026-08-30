# Naive likeness score

`assemble_naive_likeness_inputs.R` combines the public raw-count files according to a count manifest and aligns them to sample metadata. The input contract and a nine-dataset manifest template are provided under `data/`.

`train_and_evaluate_naive_likeness.R` implements the model construction and sample-level out-of-fold evaluation reported in Methods 5.5. It performs TMM normalization and logCPM transformation, retains genes with logCPM greater than 1 in at least 20% of training samples, fits an elastic-net binomial model with `alpha = 0.5`, calculates naive-class probabilities, and saves the fitted model bundle.

The final model uses 10-fold internal cross-validation for `lambda.min`. The reported internal ROC procedure is reproduced by an outer sample-level five-fold split, with a five-fold `cv.glmnet` fit inside each training partition.

`score_naive_likeness.R` applies the saved model bundle to a new raw-count matrix after the same TMM/logCPM transformation. A pre-fitted naïve-likeness model object is not included in the repository.
