function [prediction,scores] = predictCustomPLS(model,X,nComponents)
% <header>  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%MACA%%
% PREDICTCUSTOMPLS Predict responses and scores using a saved custom PLS model.
%--------------------------------------------------------------------------
% Syntax
%   prediction = softsensor.predictCustomPLS(model, X)
%   [prediction, scores] = softsensor.predictCustomPLS(model, X, nComponents)
%
% Description
%   Reuse fitted regression coefficients, calibration means and scales.
%   Project predictors into latent scores using the stored score rotations.
%
% Inputs
%   model       - Scalar struct returned by customPLS or customPLSCV.
%   X           - Real double M-by-P matrix in the fitted predictor order.
%   nComponents - Optional nonnegative integer. Default: model.nComponents.
%                 Requests above the fitted count use all fitted components.
%
% Outputs
%   prediction - M-by-1 response predictions in the fitted response units.
%   scores     - M-by-C scores, where C is the component count actually used.
%
% Notes
%   - X must have model.nPredictors columns; column names are not checked.
%   - Any NaN or Inf in a row produces NaN predictions and scores for that row.
%   - Component zero gives the calibration mean and an M-by-0 score matrix.
%   - Selected models retain full predictor indexing; do not reduce X again.
%   - No file writes or changes to RNG, path or working folder.
%
% Example
%   X = (1:10)';
%   model = softsensor.customPLS(X, 3*X + 4, 1);
%   [prediction, scores] = softsensor.predictCustomPLS(model, [11; 12]);
%
% Dependencies
%   Base MATLAB; model fields populated by customPLS/customPLSCV.
%
% Metadata
%   Project:               SoftSensorFrameworkPP3
%   Documentation updated: 2026-10-07
%   Status:                Custom PLS core; artificial-data checks passed
%
% See also customPLS, customPLSCV.
%--------------------------------------------------------------------------
% </header> %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%MACA%%
    arguments
        model (1,1) struct
        X (:,:) double {mustBeReal}
        nComponents (1,1) double {mustBeInteger,mustBeNonnegative,mustBeFinite} = model.nComponents
    end
    if size(X,2) ~= model.nPredictors
        error("CustomPLS:Predictors","X must have the fitted predictor count and order.");
    end
    count = min(nComponents,model.nComponents);
    finite = all(isfinite(X),2);
    prediction = NaN(size(X,1),1);
    prediction(finite) = [ones(sum(finite),1),X(finite,:)]*model.betaPath(:,count+1);
    scores = NaN(size(X,1),count);
    scores(finite,:) = ((X(finite,:)-model.xMean)./model.xScale)*model.rotations(:,1:count);
end
