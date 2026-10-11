function prediction = predictFittedModel(model,X,options)
% <header> %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%MACA%%
% PREDICTFITTEDMODEL Reuse a fitted regressor without refitting preprocessing.
% model is a model from fitAndEvaluate; X is M-by-P in the fitted column order.
% Group is an M-by-1 string column required for Groupwise PLS. Unseen groups
% raise SoftSensor:UnseenGroup. Incomplete X rows return NaN predictions.
% Empty M-by-P input returns an empty column. No file reads/writes or RNG changes.
% </header> %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%MACA%%
    arguments
        model (1,1) struct
        X (:,:) double {mustBeReal}
        options.Group (:,1) string = strings(0,1)
    end
    if size(X,2)~=model.nPredictors
        error("SoftSensor:Predictors","X must have the fitted predictor count and column order.");
    end
    if ~ismember(model.name,["PLS","Groupwise PLS","Mean","ElasticNet","Tree Ensemble"])
        error("SoftSensor:Model","Unsupported saved model family.");
    end
    finite = all(isfinite(X),2);
    prediction = NaN(size(X,1),1);
    if model.name=="Groupwise PLS"
        groups = options.Group;
        if numel(groups)~=size(X,1) || any(ismissing(groups) | strlength(groups)==0)
            error("SoftSensor:Groups","Provide one defined group label per prediction row.");
        end
        if any(~ismember(groups,model.groupNames))
            error("SoftSensor:UnseenGroup","A prediction group has no fitted model.");
        end
        for index = 1:numel(model.groupNames)
            rows = groups==model.groupNames(index);
            prediction(rows) = softsensor.predictFittedModel(model.estimators{index},X(rows,:));
        end
    elseif model.name=="PLS"
        prediction = softsensor.predictCustomPLS(model.estimator,X);
    elseif model.name=="Mean"
        prediction(finite) = model.yMean;
    elseif any(finite)
        scaled = (X(finite,model.keep)-model.mean(model.keep))./model.scale(model.keep);
        switch model.name
            case "Tree Ensemble"
                standardized = predict(model.estimator,scaled);
            case "ElasticNet"
                standardized = [ones(sum(finite),1),scaled]*model.estimator;
            otherwise
                error("SoftSensor:Model","Unsupported saved model family.");
        end
        prediction(finite) = standardized*model.yScale+model.yMean;
    end
end
