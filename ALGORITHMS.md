# Algorithms and data contract

## Numeric interfaces

`fitAndEvaluate(X,y,batchId,heldOut,...)` accepts real numeric predictors in
observation rows and predictor columns, a single measured response column,
string-convertible batch IDs, and a logical or numeric zero/one test mask.
MATLAB's argument declarations convert numeric predictors/responses to double
and compatible labels to strings. Supply columns explicitly. Inputs must share
row counts and include at least one predictor. Batch and optional group labels
must be nonempty and nonmissing, including rows later discarded for missing data.

Nonfinite predictor/response rows are removed together. `RowIndex` refers to
the original input, with retained held-out observations in their original order.
At least three complete calibration rows and one complete test row are required.
Overlap between calibration and test batch IDs is rejected before filtering.
Each calibration group must contain two or more batches, with at least two
training observations remaining in every CV split. One batch may include several
groups; it stays intact within each fitted estimator's folds.

The default `Models` order is `PLS`, `Groupwise PLS`, `ElasticNet`, `Tree Ensemble`.
There must be at least one model and no duplicates. Prediction columns use
`matlab.lang.makeValidName`: `PLS`, `GroupwisePLS`, `ElasticNet`, `TreeEnsemble`.
Metric rows and saved models follow the requested order. Groupwise PLS fits one
separate PLS estimator per arbitrary string group; it is not a multilevel PLS
algorithm. A test/inference group absent from calibration is rejected.

Inference requires the fitted predictor count and order. No column-name matching
or feature construction occurs. Nonfinite rows yield NaN; an empty M-by-P matrix
with M=0 yields a 0-by-1 prediction. Grouped inference still validates its labels.

## Fold isolation and tuning

Calibration batches are enumerated in first-appearance order and allocated
round-robin to `min(NumFolds, numberOfBatches)` folds. `NumFolds` defaults to four.
Allocation is deterministic; input batch order affects folds. CV errors pool
individual observations, so larger batches receive more weight. These errors
select hyperparameters and are not an independent assessment of the chosen
model. Held-out batches provide the final metrics.

PLS candidates run from one to `min(MaxComponents,P,N-1)` requested components;
the default maximum is five. Actual extraction may stop earlier. The minimum
pooled CV RMSE selects a candidate; exact ties use the first candidate. Means,
scales, and any low-level optional selection/block balancing are learned solely
from each fold's training rows. Final fitting uses complete calibration rows.

Without variable reduction, `customPLSCV` fits once per fold at the maximum
requested count and reuses its coefficient path for candidate prefixes. With
`Reduce=true`, selection depends on the candidate count and uses separate fits.
Selection discards low coefficient magnitudes (32nd percentile), or jointly low
VIP (<0.999) and low selectivity (32nd percentile). At least one predictor remains.
The returned model expands selected predictor fields back into original column
positions. Omitted coefficient contributions are zero.

ElasticNet standardizes predictors and the response from fold-training data,
then calls `lasso` with `Standardize=false`. It tests lambdas
`[0.005,0.02,0.08,0.3]` in standardized response units and alphas `[0.2,0.5,0.8]`.
Tree Ensemble uses bagged regression trees (`fitrensemble`, `Method="Bag"`),
40 learning cycles, and minimum leaf sizes `[2,5,10]`. It is a bagged tree
ensemble; no claim is made that it reproduces the historical random forest.
Both families drop calibration-constant columns. A constant response or entirely
constant predictors use the calibration response mean. PLS handles these cases
with zero extracted components. Grouped means are fitted independently.

The fitter and generator temporarily seed the twister RNG and restore the caller
state on success or failure. The seed controls stochastic fits but not batch-fold
allocation. Hyperparameter grids are intentionally small and fixed.

## PLS1 equations

`customPLS` centers X and y, optionally applies sample-standard-deviation scaling,
and optionally divides each labeled block by its centered matrix spectral norm.
Block labels, including zero, identify actual blocks; one label is required per
predictor. Constant columns have scale one before any block scaling.

From residual matrices E and F, each NIPALS component uses
`w = E'*F / norm(E'*F)`, `t = E*w`, `p = E'*t/(t'*t)`, and
`q = F'*t/(t'*t)`, followed by deflation `E = E-t*p'`, `F = F-t*q`.
Near-zero covariance/score energy stops extraction. Rotations are
`W/(P'*W)`. Regression coefficients use these rotations, response loadings,
and fitted scales, with an intercept in original predictor/response units.
The coefficient path starts with the mean-only model. Prediction requests above
the fitted component count use the fitted count; component zero uses the mean.

VIP uses component-wise explained response energy. Selectivity ratios use
target-projection explained/residual predictor energies with a denominator guard.
Saved low-level PLS structs include calibration diagnostics; high-level fitted
models retain only fields needed for inference. Tree estimators are compacted.

`spectralSNV` centers and scales each complete spectrum independently using its
sample standard deviation. Constant spectra become zero; incomplete spectra
become NaN. No cross-row statistics are fitted.

## Metrics and artificial example

Metrics exclude nonfinite measured/predicted pairs together:

- RMSE: `sqrt(mean((y-yhat).^2))`.
- MAE: `mean(abs(y-yhat))`.
- WAPE: `sum(abs(y-yhat))/sum(abs(y))`, a fraction; NaN for an all-zero denominator.
- R2Y: `1-sum((y-yhat).^2)/sum((y-mean(y)).^2)`; NaN for a constant measured response.

Empty finite pairs yield NaN metrics. Negative held-out R2 is possible, including
for grouped models trained on few artificial batches, and is not a software error.

The generator uses an explicit artificial signal model, starting at 2020-01-01.
Defaults are seed 42, 18 batches, and 240 minute observations per batch. Batch
count must be a multiple of six and at least 18; batch length must be at least
72 minutes. Laboratory sampling starts at offset 60 every 12 minutes with random
missingness. Missing sensors/outages are also inserted. Short batches can have
too few complete observations for grouped CV; generation does not guarantee a
usable fitting split. Final six batches are marked test batches.

The demo uses the default batch size, aligns process/spectral features to a fixed
45-minute within-batch residence time, and uses fixed category indicators without
learning a vocabulary from test data. It includes no learned imputation or time
optimization. Finite floating-point inputs should have numerically reasonable
scales; these routines do not guarantee stability under arithmetic overflow or
all ill-conditioned configurations.
