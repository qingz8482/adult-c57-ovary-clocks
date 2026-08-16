# Integrated PCA/glmnet clock

The integrated PCA/glmnet model scales the training-selected genes, fits PCA within the training data, and uses Gaussian `glmnet` on retained principal components. The frozen model uses the Version B age/source-batch-penalized selector and retains PC1–PC7.

Scripts here calculate component contributions, component-level and integrated gene ranks, and pathway enrichment. PCA component signs are model directions, not direct biological up- or down-regulation.
