function metricsTable = evaluateMetrics(yTrue,yPred)
% <header>  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%MACA%%
% EVALUATEMETRICS Compute finite-pair regression errors and centered R2.
%--------------------------------------------------------------------------
% Syntax
%   metricsTable = softsensor.evaluateMetrics(yTrue, yPred)
%
% Description
%   Convert equal-length response vectors to columns, remove nonfinite
%   pairs, and calculate RMSE, MAE, fractional WAPE and centered R2Y.
%
% Inputs
%   yTrue - Real numeric measured-response vector, either orientation or empty.
%   yPred - Real numeric prediction vector with the same number of elements.
%
% Outputs
%   metricsTable - One-row table named 'Value' with RMSE, MAE, WAPE and R2Y.
%                  RMSE/MAE use response units; WAPE is a fraction.
%
% Notes
%   - NaN/Inf pairs are excluded together. Invalid vector shapes or unequal
%     lengths raise Metrics:Shape; all nonfinite/empty pairs produce NaN metrics.
%   - All-zero measured responses give NaN WAPE; constant responses give NaN R2Y.
%   - R2Y = 1 - sum((yTrue-yPred).^2)/sum((yTrue-mean(yTrue)).^2).
%   - No file writes or changes to RNG, path or working folder.
%
% Example
%   metricsTable = softsensor.evaluateMetrics([100; 110; 120], [100, 110, 110]);
%
% Dependencies
%   Base MATLAB; numeric validation, finite-pair filtering and table construction.
%
% Metadata
%   Project:               SoftSensorFrameworkPP3
%   Documentation updated: 2026-10-07
%   Status:                Public metrics helper; shape/edge-case checks passed
%
% See also softsensor.fitAndEvaluate, softsensor.runSyntheticDemo.
%--------------------------------------------------------------------------
% </header> %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%MACA%%
    arguments
        yTrue {mustBeNumeric,mustBeReal}
        yPred {mustBeNumeric,mustBeReal}
    end
    if (~isvector(yTrue) && ~isempty(yTrue)) || (~isvector(yPred) && ~isempty(yPred)) ...
            || numel(yTrue) ~= numel(yPred)
        error("Metrics:Shape", "Provide numeric vectors with the same number of elements.");
    end
    yTrue = double(yTrue(:));
    yPred = double(yPred(:));
    finitePairs = isfinite(yTrue) & isfinite(yPred);
    yTrue = yTrue(finitePairs);
    yPred = yPred(finitePairs);
    rmse = NaN;
    mae = NaN;
    wape = NaN;
    r2 = NaN;
    if ~isempty(yTrue)
        residual = yTrue-yPred;
        rmse = sqrt(mean(residual.^2));
        mae = mean(abs(residual));
        denominator = sum(abs(yTrue));
        if denominator > 0
            wape = sum(abs(residual))/denominator;
        end
        totalSquares = sum((yTrue-mean(yTrue)).^2);
        if totalSquares > 0
            r2 = 1-sum(residual.^2)/totalSquares;
        end
    end
    metricsTable = table(rmse,mae,wape,r2,RowNames={'Value'}, ...
        VariableNames={'RMSE','MAE','WAPE','R2Y'});
end
