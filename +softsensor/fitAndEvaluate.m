function result = fitAndEvaluate(X,y,batchId,heldOut,options)
% <header> %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%MACA%%
% FITANDEVALUATE Tune regressors on calibration batches and evaluate test batches.
% X is N-by-P numeric predictors; y, batchId and heldOut are aligned N-by-1
% columns. NaN/Inf predictor/response rows are dropped together. Batch IDs
% must be nonmissing and a batch cannot belong to both calibration and test.
% heldOut accepts logical values or numeric zero/one only.
% Group is required for Groupwise PLS; arbitrary nonempty string labels work.
% Models defaults to all four families; select PLS alone for base MATLAB.
% For the default families, Group is required by Groupwise PLS.
% NumFolds=4, MaxComponents=5 and Seed=42 control tuning and reproducibility.
% Scaling is learned within each fold; CV RMSE is a tuning statistic.
% Predictors must be available at prediction time. Any caller-side learned
% preprocessing must also respect the calibration/validation/test boundary.
% No files are read/written. Returned models support predictFittedModel.
% Example: r = softsensor.fitAndEvaluate(X,y,batches,test,Models="PLS");
% </header> %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%MACA%%
    arguments
        X (:,:) double {mustBeReal}
        y (:,1) double {mustBeReal}
        batchId (:,1) string
        heldOut (:,1) {mustBeNumericOrLogical,mustBeReal,mustBeFinite}
        options.Group (:,1) string = strings(0,1)
        options.Models (1,:) string {mustBeMember(options.Models, ...
            ["PLS","Groupwise PLS","ElasticNet","Tree Ensemble"])} = ...
            ["PLS","Groupwise PLS","ElasticNet","Tree Ensemble"]
        options.NumFolds (1,1) double {mustBeInteger,mustBeFinite, ...
            mustBeGreaterThanOrEqual(options.NumFolds,2)} = 4
        options.MaxComponents (1,1) double {mustBeInteger,mustBePositive,mustBeFinite} = 5
        options.Seed (1,1) double {mustBeInteger,mustBeNonnegative,mustBeFinite, ...
            mustBeLessThanOrEqual(options.Seed,4294967295)} = 42
    end
    n = size(X,1);
    if size(X,2)<1 || numel(y)~=n || numel(batchId)~=n || numel(heldOut)~=n
        error("SoftSensor:Shape","X, y, batchId and heldOut must share rows; X needs predictors.");
    end
    if any(heldOut~=0 & heldOut~=1)
        error("SoftSensor:SplitMask","heldOut must contain only logical values or numeric zero/one.");
    end
    heldOut = logical(heldOut);
    if isempty(options.Models) || numel(unique(options.Models))~=numel(options.Models)
        error("SoftSensor:Models","Choose at least one model, without duplicates.");
    end
    if any(ismissing(batchId) | strlength(batchId)==0)
        error("SoftSensor:BatchID","Batch identifiers must be nonempty and nonmissing.");
    end
    if ~isempty(intersect(batchId(heldOut),batchId(~heldOut)))
        error("SoftSensor:BatchOverlap","A batch cannot occur in both calibration and test.");
    end
    groups = options.Group;
    if isempty(groups)
        if any(options.Models=="Groupwise PLS")
            error("SoftSensor:Groups","Supply Group for groupwise fitting.");
        end
        groups = repmat("Global",n,1);
    end
    if numel(groups)~=n || any(ismissing(groups) | strlength(groups)==0)
        error("SoftSensor:Groups","Group must contain one nonmissing, nonempty label per row.");
    end
    valid = isfinite(y) & all(isfinite(X),2);
    rowIndex = find(valid);
    X = X(valid,:);
    y = y(valid);
    batchId = batchId(valid);
    heldOut = heldOut(valid);
    groups = groups(valid);
    training = ~heldOut;
    if sum(training)<3 || ~any(heldOut)
        error("SoftSensor:Split","Need at least three complete calibration rows and one test row.");
    end
    calibrationBatches = unique(batchId(training),"stable");
    if numel(calibrationBatches)<2
        error("SoftSensor:Folds","Need at least two complete calibration batches.");
    end
    if any(options.Models=="Groupwise PLS") && any(~ismember(groups(heldOut),unique(groups(training))))
        error("SoftSensor:UnseenGroup","A test group has no complete calibration observations.");
    end
    previousRng = rng();
    restoreRng = onCleanup(@() rng(previousRng));
    rng(options.Seed,"twister");
    predictions = NaN(sum(heldOut),numel(options.Models));
    models = cell(numel(options.Models),1);
    tuning = table();
    metrics = table();
    for modelIndex = 1:numel(options.Models)
        modelName = options.Models(modelIndex);
        if modelName=="Groupwise PLS"
            groupNames = unique(groups(training),"stable");
            estimators = cell(numel(groupNames),1);
            for groupIndex = 1:numel(groupNames)
                rows = training & groups==groupNames(groupIndex);
                folds = makeFolds(batchId(rows),options.NumFolds);
                [estimators{groupIndex},localTuning] = selectModel("PLS",X(rows,:),y(rows),folds,options);
                localTuning.Model(:) = modelName;
                localTuning.Group(:) = groupNames(groupIndex);
                tuning = [tuning;localTuning]; %#ok<AGROW>
            end
            model = struct("name",modelName,"nPredictors",size(X,2), ...
                "groupNames",groupNames,"estimators",{estimators});
        else
            folds = makeFolds(batchId(training),options.NumFolds);
            [model,localTuning] = selectModel(modelName,X(training,:),y(training),folds,options);
            tuning = [tuning;localTuning]; %#ok<AGROW>
        end
        predictions(:,modelIndex) = softsensor.predictFittedModel(model,X(heldOut,:),Group=groups(heldOut));
        models{modelIndex} = model;
        if any(~isfinite(predictions(:,modelIndex)))
            error("SoftSensor:Prediction","A fitted model produced nonfinite complete-row predictions.");
        end
        row = softsensor.evaluateMetrics(y(heldOut),predictions(:,modelIndex));
        row.Properties.RowNames = {};
        row = addvars(row,modelName,sum(heldOut),Before=1,NewVariableNames={'Model','N'});
        metrics = [metrics;row]; %#ok<AGROW>
    end
    predictionTable = table(rowIndex(heldOut),batchId(heldOut),groups(heldOut),y(heldOut), ...
        VariableNames={'RowIndex','BatchID','Group','Measured'});
    predictionNames = matlab.lang.makeValidName(cellstr(options.Models));
    predictionTable = [predictionTable,array2table(predictions,VariableNames=predictionNames)];
    result = struct("schemaVersion",1,"seed",options.Seed,"modelNames",options.Models, ...
        "metrics",metrics,"predictions",predictionTable,"tuning",tuning,"models",{models}, ...
        "calibrationBatches",calibrationBatches,"testBatches",unique(batchId(heldOut),"stable"), ...
        "inputRows",n,"completeRows",sum(valid),"droppedRows",find(~valid), ...
        "settings",options,"matlabVersion",version());
    % Group labels are data, not configuration; do not duplicate the full input.
    result.settings = rmfield(result.settings,"Group");
end

function folds = makeFolds(batches,requested)
    names = unique(batches,"stable");
    if numel(names)<2
        error("SoftSensor:Folds","Every fitted group needs at least two complete calibration batches.");
    end
    folds = zeros(numel(batches),1);
    count = min(requested,numel(names));
    for index = 1:numel(names)
        folds(batches==names(index)) = mod(index-1,count)+1;
    end
    if any(arrayfun(@(fold) sum(folds~=fold)<2,unique(folds)))
        error("SoftSensor:FoldTraining","Every validation fold must leave at least two training rows.");
    end
end

function [model,tuning] = selectModel(name,X,y,folds,options)
    if name=="PLS"
        [core,~] = softsensor.customPLSCV(X,y,folds,MaxComponents=options.MaxComponents,Scale=true);
        parameters = (1:numel(core.cv.rmse))';
        errors = core.cv.rmse';
        selected = parameters==core.cv.selectedRequestedComponents;
        candidates = table(parameters,zeros(size(parameters)),VariableNames={'Parameter','Alpha'});
        model = struct("name","PLS","nPredictors",size(X,2),"estimator",compactPLS(core));
    else
        switch name
            case "ElasticNet"
                [lambda,alpha] = ndgrid([0.005,0.02,0.08,0.3],[0.2,0.5,0.8]);
                candidates = table(lambda(:),alpha(:),VariableNames={'Parameter','Alpha'});
            case "Tree Ensemble"
                candidates = table([2;5;10],zeros(3,1),VariableNames={'Parameter','Alpha'});
            otherwise
                error("SoftSensor:Model","Unsupported model family.");
        end
        errors = NaN(height(candidates),1);
        for index = 1:height(candidates)
            predicted = NaN(size(y));
            for fold = unique(folds)'
                validation = folds==fold;
                fitted = fitModel(name,X(~validation,:),y(~validation), ...
                    candidates.Parameter(index),candidates.Alpha(index));
                predicted(validation) = softsensor.predictFittedModel(fitted,X(validation,:));
            end
            errors(index) = sqrt(mean((y-predicted).^2));
        end
        if any(~isfinite(errors))
            error("SoftSensor:Tuning","Candidate validation errors must be finite.");
        end
        [~,best] = min(errors);
        selected = (1:height(candidates))'==best;
        model = fitModel(name,X,y,candidates.Parameter(best),candidates.Alpha(best));
    end
    tuning = addvars(candidates,repmat(name,height(candidates),1), ...
        repmat("Global",height(candidates),1),errors,selected,Before=1, ...
        NewVariableNames={'Model','Group','BatchCV_RMSE','Selected'});
end

function model = fitModel(name,X,y,parameter,alpha)
    scales = std(X,0,1);
    model = struct("name",name,"nPredictors",size(X,2),"mean",mean(X,1), ...
        "scale",scales,"keep",scales>0,"yMean",mean(y),"yScale",std(y), ...
        "parameter",parameter,"alpha",alpha,"estimator",[]);
    model.scale(~model.keep) = 1;
    if model.yScale==0 || ~any(model.keep)
        model.name = "Mean";
        return
    end
    scaled = (X(:,model.keep)-model.mean(model.keep))./model.scale(model.keep);
    scaledY = (y-model.yMean)/model.yScale;
    switch name
        case "ElasticNet"
            [beta,info] = lasso(scaled,scaledY,Lambda=parameter,Alpha=alpha,Standardize=false);
            model.estimator = [info.Intercept;beta];
        case "Tree Ensemble"
            learner = templateTree(MinLeafSize=parameter);
            model.estimator = compact(fitrensemble(scaled,scaledY,Method="Bag", ...
                NumLearningCycles=40,Learners=learner));
        otherwise
            error("SoftSensor:Model","Unsupported fitted family.");
    end
end

function compact = compactPLS(model)
    fields = {'nComponents','nPredictors','xMean','xScale','yMean','betaPath','rotations'};
    compact = struct();
    for index = 1:numel(fields)
        compact.(fields{index}) = model.(fields{index});
    end
end
