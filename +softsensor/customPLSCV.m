function [model,prediction] = customPLSCV(X,y,foldId,options)
% <header>  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%MACA%%
% CUSTOMPLSCV Tune custom PLS components with calibration-only validation folds.
%--------------------------------------------------------------------------
% Syntax
%   [model, prediction] = softsensor.customPLSCV(X, y, foldId)
%   [model, prediction] = softsensor.customPLSCV(X, y, foldId, Name=Value, ...)
%
% Description
%   Fit each component candidate within each training fold, predict its
%   validation rows, select minimum pooled RMSE, then refit on all calibration
%   rows. Optional preprocessing and variable selection are fitted within folds.
%
% Inputs
%   X      - Finite real double N-by-P calibration predictor matrix.
%   y      - Finite real double N-by-1 response column aligned with X.
%   foldId - Positive integer double N-by-1 fold assignment. At least two
%            distinct folds, each leaving at least two training observations.
%
% Name-Value Arguments
%   MaxComponents - Positive integer candidate limit. Default: 20; capped at
%                   min(MaxComponents, P, N-1). Fits can stop earlier by rank.
%   Reduce        - Logical coefficient/VIP/selectivity reduction flag. Default: false.
%   Scale         - Logical predictor standard-deviation scaling flag. Default: false.
%   Blocks        - Row vector of nonnegative integer predictor block labels.
%                   Default: []. One label per original predictor when supplied.
%
% Outputs
%   model      - Final customPLS struct with cv.foldId, cv.predictions,
%                cv.predictionPath, cv.rmse, cv.r2, cv.selectedRequestedComponents
%                and cv.type. rmsecv/r2cv contain the candidate tuning curves.
%   prediction - N-by-1 out-of-fold predictions at the selected candidate count.
%
% Notes
%   - Ties in pooled RMSE select the first candidate count.
%   - These errors tune the model; use separate data for final performance claims.
%   - Keep all repeated observations of a sample/batch in the same input fold.
%   - Reduction uses coefficient/selectivity 32nd percentiles and VIP < 0.999;
%     at least one predictor is retained. Selection is repeated in each fold.
%   - Reduced models expand predictor-indexed fields to the original P columns;
%     predictorMask marks retained columns, and omitted contributions are zero.
%   - No file writes or changes to RNG, path or working folder.
%
% Example
%   X = [(1:20)', sin((1:20)')];
%   y = 3*X(:,1) + X(:,2);
%   foldId = repelem((1:4)', 5);
%   [model, oof] = softsensor.customPLSCV(X, y, foldId, MaxComponents=2);
%
% Dependencies
%   Base MATLAB; customPLS, predictCustomPLS.
%
% Metadata
%   Project:               SoftSensorFrameworkPP3
%   Documentation updated: 2026-10-07
%   Status:                Custom PLS core; fold-isolation checks passed
%
% See also customPLS, predictCustomPLS.
%--------------------------------------------------------------------------
% </header> %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%MACA%%
    arguments
        X (:,:) double {mustBeReal,mustBeFinite}
        y (:,1) double {mustBeReal,mustBeFinite}
        foldId (:,1) double {mustBeInteger,mustBePositive,mustBeFinite}
        options.MaxComponents (1,1) double {mustBeInteger,mustBePositive,mustBeFinite} = 20
        options.Reduce (1,1) logical = false
        options.Scale (1,1) logical = false
        options.Blocks (1,:) double {mustBeInteger,mustBeNonnegative,mustBeFinite} = []
    end
    if size(X,1) ~= numel(y) || numel(foldId) ~= numel(y)
        error("CustomPLS:FoldShape","X, y and foldId must have the same row count.");
    end
    if ~isempty(options.Blocks) && numel(options.Blocks) ~= size(X,2)
        error("CustomPLS:Blocks","Provide one block number per predictor.");
    end
    folds = unique(foldId);
    if numel(folds)<2
        error("CustomPLS:Folds","Use at least two validation folds.");
    end
    maximum = min([options.MaxComponents,size(X,2),size(X,1)-1]);
    if maximum < 1
        error("CustomPLS:Shape","PLS tuning requires observations and predictors.");
    end
    predictions = NaN(numel(y),maximum);
    for fold = folds'
        validation = foldId==fold;
        if sum(~validation)<2
            error("CustomPLS:FoldTraining","Every validation fold needs at least two training rows.");
        end
        % Without variable selection, one fit supplies every prefix of betaPath.
        if ~options.Reduce
            fitted = fitSelected(X(~validation,:),y(~validation),maximum,options);
        end
        for component = 1:maximum
            if options.Reduce
                fitted = fitSelected(X(~validation,:),y(~validation),component,options);
            end
            predictions(validation,component) = softsensor.predictCustomPLS(fitted,X(validation,:),component);
        end
    end
    rmse = sqrt(mean((y-predictions).^2,1));
    [~,best] = min(rmse);
    model = fitSelected(X,y,best,options);
    prediction = predictions(:,best);
    denominator = sum((y-mean(y)).^2);
    r2 = NaN(size(rmse));
    if denominator > 0
        r2 = 1-sum((y-predictions).^2,1)/denominator;
    end
    model.cv = struct("foldId",foldId,"predictions",prediction,"predictionPath",predictions,"rmse",rmse, ...
        "r2",r2,"selectedRequestedComponents",best,"type","component-tuning cross-validation");
    model.rmsecv = rmse;
    model.r2cv = r2;
end

function model = fitSelected(X,y,nComponents,options)
% FITSELECTED Fit optional coefficient/VIP/target-projection selection.
    keep = true(size(X,2),1);
    if options.Reduce
        initial = softsensor.customPLS(X,y,nComponents,Scale=options.Scale,Blocks=options.Blocks);
        smallCoefficient = abs(initial.reg)<quantile(abs(initial.reg),0.32);
        lowVip = initial.vip<0.999;
        lowSelectivity = initial.selratio'<quantile(abs(initial.selratio),0.32);
        keep = ~(smallCoefficient | (lowVip & lowSelectivity));
        if ~any(keep)
            [~,strongest] = max(abs(initial.reg));
            keep(strongest) = true;
        end
    end
    blocks = options.Blocks;
    if ~isempty(blocks)
        blocks = blocks(keep);
    end
    model = softsensor.customPLS(X(:,keep),y,nComponents,Scale=options.Scale,Blocks=blocks);
    if all(keep)
        return
    end
    fullCount = size(X,2);
    betaPath = zeros(fullCount+1,size(model.betaPath,2));
    betaPath([true;keep],:) = model.betaPath;
    model.betaPath = betaPath;
    model.beta = betaPath(:,end);
    model.xMean = mean(X,1);
    scale = ones(1,fullCount);
    scale(keep) = model.xScale;
    model.xScale = scale;
    model.xWeights = expandPredictors(model.xWeights,keep);
    model.xLoadings = expandPredictors(model.xLoadings,keep);
    model.rotations = expandPredictors(model.rotations,keep);
    model.vip = expandPredictors(model.vip,keep);
    model.reg = expandPredictors(model.reg,keep);
    model.selratio = expandPredictors(model.selratio',keep)';
    model.wts = model.xWeights;
    model.loads{2,1} = model.xLoadings;
    model.predictorMask = keep;
    model.nPredictors = fullCount;
end

function full = expandPredictors(reduced,keep)
% EXPANDPREDICTORS Retain original predictor indexing for saved models/plots.
    full = zeros(numel(keep),size(reduced,2));
    full(keep,:) = reduced;
end
