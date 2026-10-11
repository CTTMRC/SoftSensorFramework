# Soft Sensor MATLAB Release

A standalone regression library with a numerical data interface, batch-held-out
evaluation, reusable fitted models, and an artificial process/spectral example.
Licensed under MIT. This source package contains no research datasets, fitted
research models, workbooks, or study figures.

## Requirements and first run

Verified on MATLAB R2025a Update 1 on Windows. Other versions and operating
systems have not been verified. Base MATLAB supports the custom PLS core,
Groupwise PLS, generator, metrics, preprocessing, and export helpers.
ElasticNet, Tree Ensemble, and reference comparisons require Statistics and
Machine Learning Toolbox. No PLSToolbox or downloaded dependencies are needed.

Extract the archive, open MATLAB in its `PublicationReady` folder, and run:

```matlab
addpath(pwd);                       % Add this release root once.
result = softsensor.runSyntheticDemo();
verification = softsensor.verifyRelease();
```

For base MATLAB:

```matlab
result = softsensor.runSyntheticDemo(Models=["PLS","Groupwise PLS"]);
verification = softsensor.verifyRelease(IncludeToolboxTests=false);
```

The demo uses artificial batches, seven fixed spectral bands, seventeen process
variables, and fixed artificial category indicators. It supports `Target="Mw"`
(default) and `Target="SmallMolecules"`. It prints metrics and writes nothing
by default. Its results demonstrate software behavior and do not reproduce or
validate the original scientific study.

## Fit your numeric data

```matlab
% Supply aligned data that you are authorized to use:
% X: N-by-P predictors; y: N-by-1 measured response.
% batchId: N-by-1 strings; heldOut: N-by-1 logical test mask.
result = softsensor.fitAndEvaluate(X,y,batchId,heldOut,Models="PLS");

% Fit all four families by also supplying an N-by-1 Group string column:
result = softsensor.fitAndEvaluate(X,y,batchId,heldOut,Group=groupLabels);

% Models are in result.modelNames order; Xnew uses the fitted predictor order.
prediction = softsensor.predictFittedModel(result.models{1},Xnew);
% A grouped fitted model requires Group=newGroupLabels at inference.
```

Predictors must be available at inference time. Fit any caller-side learned
imputation, spectral transformation, encoding, or feature selection within
calibration folds. This library cannot undo leakage already present in its inputs.

Incomplete predictor/response rows are dropped together and reported in
`result.droppedRows`. Batch IDs must be defined before row filtering; a batch
cannot occur in both calibration and test. Each fitted group needs at least two
complete calibration batches, and every validation fold must leave at least two
training rows. New groups at inference raise an error. The exact contract,
tuning grids, metrics, and model limitations are in [ALGORITHMS.md](ALGORITHMS.md).

## Saving and synthetic export

```matlab
% Choose a generated-output location outside the source release:
runFolder = softsensor.saveResult(result,"../generated/results");
loaded = load(fullfile(runFolder,"result.mat"),"result");
prediction = softsensor.predictFittedModel(loaded.result.models{1},Xnew);

% Or generate, fit, and save one artificial result:
result = softsensor.runSyntheticDemo(OutputFolder="../generated/results");

% Export an artificial dataset separately; MAT is the default:
[Xp,Xs,Y,metadata] = softsensor.generateSyntheticData( ...
    OutputFolder="../generated/synthetic");
% ExportFormat="csv" or "both" is opt-in.
```

Every requested save creates a new directory under the supplied parent. Writes
are staged in a sibling `.partial` directory and renamed only after success.
Consumers must ignore every `.partial` directory, including one left by an
interrupted rename. Completed result folders contain `result.mat`, three small
CSV tables, and `COMPLETE.json`; raw synthetic inputs are not duplicated for
each result. Dataset exports contain only the explicitly selected formats.

`saveResult` serializes the supplied struct in full. Real-data predictions,
identifiers, fitted coefficients, and additional caller fields can be private.
Keep such outputs separate from this public source package. `isSynthetic=true`
does not anonymize caller-supplied data. Save/reload verification covers this
MATLAB version; serialized toolbox models can require compatible MATLAB versions.

## Included interfaces

| Function in `+softsensor` | Purpose |
|---|---|
| `customPLS`, `predictCustomPLS` | PLS1 fit and saved predictions/scores |
| `customPLSCV` | Fold-isolated component tuning and optional predictor selection |
| `spectralSNV` | Independent row spectral normalization |
| `fitAndEvaluate`, `predictFittedModel` | Batch tuning, held-out metrics, compact inference |
| `evaluateMetrics` | RMSE, MAE, fractional WAPE, centered R2 |
| `generateSyntheticData`, `runSyntheticDemo` | Artificial data and two-target demo |
| `saveResult`, `exportSyntheticData` | Explicit result saves and dataset exports |
| `verifyCustomPLS`, `verifyRelease` | Numerical verification and release test suite |

See [RELEASE_NOTES.md](RELEASE_NOTES.md) for scope and exclusions, and
[VERIFICATION.md](VERIFICATION.md) for the recorded checks. The release contains
source, tests, documentation, and a SHA-256 manifest. Generated data and models
are deliberately absent.
