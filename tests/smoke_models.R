#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(glmnet))
set.seed(20260630)
n <- 90L
p <- 40L
x <- matrix(rnorm(n * p), nrow = n)
latent_age_axis <- rnorm(n)
x[, 1:10] <- latent_age_axis + matrix(rnorm(n * 10, sd = 0.15), nrow = n)
age <- 35 + 8 * latent_age_axis + rnorm(n, sd = 0.5)
train <- 1:60
test <- 61:90

fit_elasticnet <- cv.glmnet(x[train, ], age[train], alpha = 0.5, nfolds = 5)
pred_elasticnet <- as.numeric(predict(fit_elasticnet, x[test, ], s = "lambda.min"))
stopifnot(length(pred_elasticnet) == length(test), cor(pred_elasticnet, age[test]) > 0.8)

center <- colMeans(x[train, ])
scale <- apply(x[train, ], 2, sd)
x_train <- scale(x[train, ], center = center, scale = scale)
x_test <- scale(x[test, ], center = center, scale = scale)
pca <- prcomp(x_train, center = FALSE, scale. = FALSE)
fit_pca <- cv.glmnet(pca$x[, 1:8], age[train], alpha = 0.5, nfolds = 5)
test_scores <- predict(pca, x_test)[, 1:8]
pred_pca <- as.numeric(predict(fit_pca, test_scores, s = "lambda.min"))
stopifnot(length(pred_pca) == length(test), cor(pred_pca, age[test]) > 0.8)

cell_class <- rep(c("Granulosa", "Theca_Stroma"), each = n / 2)
predicted <- numeric(n)
for (cls in unique(cell_class)) {
  idx <- which(cell_class == cls)
  fold_test <- tail(idx, 10)
  fold_train <- setdiff(idx, fold_test)
  fit <- cv.glmnet(x[fold_train, ], age[fold_train], alpha = 0.5, nfolds = 3)
  predicted[fold_test] <- as.numeric(predict(fit, x[fold_test, ], s = "lambda.min"))
  stopifnot(all(is.finite(predicted[fold_test])))
}
cat("synthetic model smoke tests passed\n")
