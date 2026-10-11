function verification = verifyCustomPLS(options)
% <header>  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%MACA%%
% VERIFYCUSTOMPLS Check custom PLS equations, edge cases and fold isolation.
%--------------------------------------------------------------------------
% Syntax
%   verification = softsensor.verifyCustomPLS()
%   verification = softsensor.verifyCustomPLS(Reference=true)
%
% Description
%   Verify scores, VIP normalization, exact linear/OLS cases, rank/constant
%   handling, prediction boundaries, row SNV and validation-fold label isolation.
%   Optionally compare held-out predictions with MathWorks plsregress.
%
% Inputs
%   No positional arguments; one name-value argument is supported.
%
% Name-Value Arguments
%   Reference - Logical flag enabling the independent plsregress comparison.
%               Default: false. Reference=true requires Statistics toolbox.
%
% Outputs
%   verification - Struct containing passed, algorithm, referenceChecked,
%                  referenceMaxDifference, foldIsolation and baseMatlabImplementation.
%                  referenceMaxDifference is NaN when Reference=false.
%
% Notes
%   - Artificial numeric data only; no original research data is read.
%   - Label mutation compares every candidate prediction for the held-out fold,
%     including variable selection and block scaling. It does not claim an
%     independent error estimate for the selected component count.
%   - Assertions/errors stop execution; RNG seed 17 is temporary and the caller
%     state is restored even on failure. Prints results; writes no files.
%
% Example
%   verification = softsensor.verifyCustomPLS();
%
% Dependencies
%   Base MATLAB; customPLS, customPLSCV, predictCustomPLS, spectralSNV.
%   Optional Statistics and Machine Learning Toolbox plsregress for Reference=true.
%
% References
%   Independent reference/algorithm notes: ALGORITHMS.md.
%
% Metadata
%   Project:               SoftSensorFrameworkPP3
%   Documentation updated: 2026-10-07
%   Status:                Executable custom PLS verification entry point
%
% See also customPLS, customPLSCV, predictCustomPLS, spectralSNV.
%--------------------------------------------------------------------------
% </header> %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%MACA%%
    arguments
        options.Reference (1,1) logical = false
    end
    priorRng = rng();
    restoreRng = onCleanup(@() rng(priorRng));
    rng(17,"twister");
    X = randn(90,12);
    y = X(:,1:4)*[2;-1;0.8;3] + 0.15*randn(90,1) + 10;
    newX = randn(25,12);
    maximumDifference = NaN;
    if options.Reference
        maximumDifference = 0;
    end
    for component = 1:6
        model = softsensor.customPLS(X,y,component);
        [~,scores] = softsensor.predictCustomPLS(model,X);
        assert(norm(scores-model.xScores,"fro")<1e-9);
        assert(abs(sum(model.vip.^2)-size(X,2))<1e-9);
        residual = (X-model.xMean)./model.xScale-model.xScores*model.xLoadings';
        assert(norm(model.xScores'*residual,"fro")<1e-8);
        if options.Reference
            [~,~,~,~,reference] = plsregress(X,y,component);
            difference = max(abs(softsensor.predictCustomPLS(model,newX)-[ones(25,1),newX]*reference));
            assert(difference<1e-9);
            maximumDifference = max(maximumDifference,difference);
        end
    end
    exactX = (-5:5)';
    exact = softsensor.customPLS(exactX,3*exactX+4,1);
    assert(max(abs(exact.beta-[4;3]))<1e-12);
    complete = softsensor.customPLS(X(:,1:6),y,6);
    ols = [ones(90,1),X(:,1:6)]\y;
    assert(norm(complete.beta-ols)<1e-9);
    constant = softsensor.customPLS(ones(12,3),7*ones(12,1),5);
    assert(constant.nComponents==0);
    assert(all(softsensor.predictCustomPLS(constant,ones(2,3))==7));
    fixedResponse = softsensor.customPLS(X,5*ones(90,1),4);
    assert(fixedResponse.nComponents==0 && all(fixedResponse.beta(2:end)==0));
    collinear = softsensor.customPLS([X(:,1),X(:,1),ones(90,1)],y,3);
    assert(collinear.nComponents==1 && all(isfinite(collinear.beta)));
    assert(abs(collinear.vip(1)-collinear.vip(2))<1e-12 && collinear.vip(3)==0);
    scaled = softsensor.customPLS(X,y,4,Scale=true);
    if options.Reference
        standardizedX = (X-mean(X,1))./std(X,0,1);
        [~,~,~,~,reference] = plsregress(standardizedX,y,4);
        standardizedNew = (newX-mean(X,1))./std(X,0,1);
        difference = max(abs(softsensor.predictCustomPLS(scaled,newX)-[ones(25,1),standardizedNew]*reference));
        assert(difference<1e-9);
        maximumDifference = max(maximumDifference,difference);
    end
    incomplete = newX;
    incomplete(1,1) = NaN;
    incomplete(2,1) = Inf;
    missingPrediction = softsensor.predictCustomPLS(scaled,incomplete);
    assert(all(isnan(missingPrediction(1:2))) && all(isfinite(missingPrediction(3:end))));
    expectError(@() softsensor.predictCustomPLS(scaled,newX(:,1:3)),"CustomPLS:Predictors");
    expectError(@() softsensor.customPLS(X,y(1:2),2),"CustomPLS:Shape");
    expectError(@() softsensor.customPLS(X,y,2,Blocks=1),"CustomPLS:Blocks");
    normalized = softsensor.spectralSNV([1,2,3;4,4,4;NaN,1,2]);
    assert(abs(mean(normalized(1,:)))<1e-12 && abs(std(normalized(1,:))-1)<1e-12);
    assert(all(normalized(2,:)==0) && all(isnan(normalized(3,:))));

    folds = repelem((1:5)',18);
    [selected,oof] = softsensor.customPLSCV(X,y,folds,MaxComponents=4,Reduce=true, ...
        Blocks=[ones(1,7),2*ones(1,5)]);
    changedY = y;
    changedY(folds==1) = changedY(folds==1)+1000;
    changed = softsensor.customPLSCV(X,changedY,folds,MaxComponents=4,Reduce=true, ...
        Blocks=[ones(1,7),2*ones(1,5)]);
    assert(isequal(selected.cv.predictionPath(folds==1,:),changed.cv.predictionPath(folds==1,:)), ...
        "A validation fold's labels leaked into its fitted predictions.");
    assert(all(isfinite(oof)) && all(isfinite(softsensor.predictCustomPLS(selected,newX))));
    assert(selected.nPredictors==12 && numel(selected.predictorMask)==12);
    assert(all(selected.beta([false;~selected.predictorMask])==0));
    expectError(@() softsensor.customPLSCV(X,y,ones(90,1)),"CustomPLS:Folds");
    verification = struct("passed",true,"algorithm","NIPALS PLS1", ...
        "referenceChecked",options.Reference,"referenceMaxDifference",maximumDifference, ...
        "foldIsolation",true,"baseMatlabImplementation",true);
    disp(verification);
end

function expectError(action,identifier)
% EXPECTERROR Check the public error identifier of an invalid input.
    observed = "";
    try
        action();
    catch exception
        observed = string(exception.identifier);
    end
    assert(observed==identifier,"Expected %s; received %s.",identifier,observed);
end
