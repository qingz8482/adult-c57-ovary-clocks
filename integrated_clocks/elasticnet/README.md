# Integrated ElasticNet clock

The integrated ElasticNet model is fitted directly to training-scaled selected genes with Gaussian `glmnet`. The frozen model uses the HVG selector. Scripts here interpret ElasticNet coefficients and test whether selected-gene summaries are sensitive to unequal dataset sampling.

Run shared preprocessing and training first. These scripts resolve the shared integrated-clock directory relative to the repository and consume its frozen outputs.
