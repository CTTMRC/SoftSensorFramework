function model = customPLS(X,y,nComponents,options)
% <header>  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%MACA%%
% CUSTOMPLS Fit a single-response NIPALS partial least-squares model.
%--------------------------------------------------------------------------
% Syntax
%   model = softsensor.customPLS(X, y, nComponents)
%   model = softsensor.customPLS(X, y, nComponents, Name=Value, ...)
%
% Description
%   Center calibration data, optionally scale predictors, extract PLS1
%   components, and store raw-space regression coefficients and diagnostics.
%   All preprocessing parameters are estimated from the supplied calibration rows.
%
% Inputs
%   X           - Finite real double N-by-P predictor matrix; N >= 2, P >= 1.
%   y           - Finite real double N-by-1 response column, aligned with X.
%   nComponents - Nonnegative integer upper bound on components. Extraction is
%                 limited to min(nComponents, N-1, P) and may stop earlier.
%
% Name-Value Arguments
%   Scale  - Logical flag for predictor standard-deviation scaling. Default: false.
%   Blocks - Row vector of nonnegative integer block labels, one per predictor.
%            Default: []. Each nonzero-norm block is divided by its centered
%            spectral norm after any standard-deviation scaling.
%
% Outputs
%   model - Struct containing nComponents, nPredictors, xMean, xScale, yMean;
%           beta (P+1)-by-1 and betaPath (P+1)-by-(nComponents+1), with intercepts;
%           xWeights, xLoadings, rotations (P-by-nComponents), xScores and yScores
%           (N-by-nComponents), yLoadings (1-by-nComponents), explainedVariance;
%           vip (P-by-1), selratio (1-by-P), reg and predictorMask (P-by-1).
%           Compatibility fields: algorithm, wts, loads and pred. See Notes.
%
% Notes
%   - Centering always applies; predictions retain the units supplied in y.
%   - Autoscaling uses scale one for constant columns; block balancing can
%     multiply that scale afterwards. Zero-norm blocks are left unchanged.
%   - Rank exhaustion or zero response covariance can give zero components;
%     the resulting model predicts the calibration response mean.
%   - beta uses raw X units; reg uses preprocessed predictor units. betaPath's
%     first column is the mean-only model. Scores use rotations, not loadings.
%   - VIP weights component-wise explained response energy. Selectivity ratios
%     use target-projection explained/residual energies with a denominator guard.
%   - loads = {Xscores, Yscores; Xloadings, Yloadings}; pred{2} is fitted y.
%   - No PLSToolbox calls, file writes, RNG changes or global-state changes.
%
% Example
%   X = (1:10)';
%   y = 3*X + 4;
%   model = softsensor.customPLS(X, y, 1);
%   yhat = softsensor.predictCustomPLS(model, [11; 12]);
%
% Dependencies
%   Base MATLAB. No external PLS or Statistics toolbox functions.
%
% References
%   Algorithm documentation and numerical checks: ALGORITHMS.md.
%
% Metadata
%   Project:               SoftSensorFrameworkPP3
%   Documentation updated: 2026-10-07
%   Status:                Custom PLS core; artificial-data checks passed
%
% See also predictCustomPLS, customPLSCV, spectralSNV.
%--------------------------------------------------------------------------
% </header> %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%MACA%%
    arguments
        X (:,:) double {mustBeReal,mustBeFinite}
        y (:,1) double {mustBeReal,mustBeFinite}
        nComponents (1,1) double {mustBeInteger,mustBeNonnegative,mustBeFinite}
        options.Scale (1,1) logical = false
        options.Blocks (1,:) double {mustBeInteger,mustBeNonnegative,mustBeFinite} = []
    end
    [n,p] = size(X);
    if size(y,1) ~= n || n < 2 || p < 1
        error("CustomPLS:Shape","X/y must share at least two rows and X must have a predictor.");
    end
    if ~isempty(options.Blocks) && numel(options.Blocks) ~= p
        error("CustomPLS:Blocks","Provide one block number per predictor.");
    end
    xMean = mean(X,1);
    yMean = mean(y);
    xScale = ones(1,p);
    if options.Scale
        xScale = std(X,0,1);
        xScale(xScale==0) = 1;
    end
    centered = (X-xMean)./xScale;
    for block = unique(options.Blocks)
        columns = options.Blocks==block;
        blockScale = norm(centered(:,columns),2);
        if blockScale > 0
            xScale(columns) = xScale(columns)*blockScale;
        end
    end
    X0 = (X-xMean)./xScale;
    y0 = y-yMean;
    E = X0;
    F = y0;
    maximum = min([nComponents,n-1,p]);
    W = zeros(p,maximum);
    P = zeros(p,maximum);
    T = zeros(n,maximum);
    U = zeros(n,maximum);
    Q = zeros(1,maximum);
    ssY = zeros(1,maximum);
    ssX = zeros(1,maximum);
    betaPath = zeros(p+1,maximum+1);
    betaPath(1,:) = yMean;
    extracted = 0;
    covarianceTolerance = 100*eps*norm(X0,"fro")*norm(y0);
    scoreTolerance = (100*eps*norm(X0,"fro"))^2;
    for component = 1:maximum
        direction = E'*F;
        directionNorm = norm(direction);
        if directionNorm <= covarianceTolerance || directionNorm == 0
            break
        end
        w = direction/directionNorm;
        t = E*w;
        scoreEnergy = t'*t;
        if scoreEnergy <= scoreTolerance || scoreEnergy == 0
            break
        end
        loading = E'*t/scoreEnergy;
        responseLoading = F'*t/scoreEnergy;
        W(:,component) = w;
        P(:,component) = loading;
        T(:,component) = t;
        Q(component) = responseLoading;
        if responseLoading ~= 0
            U(:,component) = F/responseLoading;
        end
        ssY(component) = scoreEnergy*responseLoading^2;
        ssX(component) = scoreEnergy*sum(loading.^2);
        E = E-t*loading';
        F = F-t*responseLoading;
        rotation = W(:,1:component)/(P(:,1:component)'*W(:,1:component));
        scaledBeta = rotation*Q(1:component)';
        rawBeta = scaledBeta./xScale';
        betaPath(:,component+1) = [yMean-xMean*rawBeta;rawBeta];
        extracted = component;
    end
    W = W(:,1:extracted);
    P = P(:,1:extracted);
    T = T(:,1:extracted);
    U = U(:,1:extracted);
    Q = Q(:,1:extracted);
    ssY = ssY(1:extracted);
    ssX = ssX(1:extracted);
    betaPath = betaPath(:,1:extracted+1);
    rotations = zeros(p,0);
    if extracted > 0
        rotations = W/(P'*W);
    end
    beta = betaPath(:,end);
    importance = zeros(p,1);
    if sum(ssY) > 0
        importance = sqrt(p*((W.^2)*ssY')/sum(ssY));
    end
    scaledBeta = beta(2:end).*xScale';
    selectivity = zeros(p,1);
    if norm(scaledBeta) > 0
        targetScores = X0*(scaledBeta/norm(scaledBeta));
        targetEnergy = targetScores'*targetScores;
        if targetEnergy > 0
            targetLoadings = targetScores'*X0/targetEnergy;
            projected = targetScores*targetLoadings;
            explainedEnergy = sum(projected.^2,1);
            residualEnergy = sum((X0-projected).^2,1);
            varying = sum(X0.^2,1)>0;
            selectivity(varying) = (explainedEnergy(varying)./ ...
                max(residualEnergy(varying),eps*sum(X0(:,varying).^2,1)))';
        end
    end
    explained = zeros(2,extracted);
    if sum(X0.^2,"all") > 0
        explained(1,:) = ssX/sum(X0.^2,"all");
    end
    if sum(y0.^2) > 0
        explained(2,:) = ssY/sum(y0.^2);
    end
    model = struct("algorithm","NIPALS PLS1","nComponents",extracted, ...
        "nPredictors",p,"xMean",xMean,"xScale",xScale,"yMean",yMean, ...
        "beta",beta,"betaPath",betaPath,"xWeights",W,"xLoadings",P, ...
        "yLoadings",Q,"xScores",T,"yScores",U,"rotations",rotations, ...
        "explainedVariance",explained,"vip",importance,"selratio",selectivity', ...
        "reg",scaledBeta,"wts",W,"loads",{{T,U;P,Q}}, ...
        "predictorMask",true(p,1),"pred",{{[],[ones(n,1),X]*beta}});
end
