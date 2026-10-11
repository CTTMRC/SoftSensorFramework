function result = runSyntheticDemo(options)
% <header> %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%MACA%%
% RUNSYNTHETICDEMO Run four model families on artificial held-out batches.
% Target is "Mw" (default) or "SmallMolecules"; Seed defaults to 42.
% Models selects any subset of fitAndEvaluate's four supported families.
% OutputFolder="" writes nothing. A supplied parent gets a new, complete run
% folder containing compact reusable models and small result CSVs. The raw
% synthetic dataset is regenerated on demand and is not duplicated per run.
% result.outputFolder identifies the committed run directory, if requested.
% Repeated calls restore caller RNG and never overwrite prior runs.
% Requires Statistics toolbox for elastic net and trees; PLS uses base MATLAB.
% These artificial results do not reproduce the original scientific study.
% Example: r = softsensor.runSyntheticDemo(OutputFolder="results");
% </header> %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%MACA%%
    arguments
        options.Target (1,1) string {mustBeMember(options.Target,["Mw","SmallMolecules"])} = "Mw"
        options.Seed (1,1) double {mustBeInteger,mustBeNonnegative,mustBeFinite, ...
            mustBeLessThanOrEqual(options.Seed,4294967295)} = 42
        options.Models (1,:) string {mustBeMember(options.Models, ...
            ["PLS","Groupwise PLS","ElasticNet","Tree Ensemble"])} = ...
            ["PLS","Groupwise PLS","ElasticNet","Tree Ensemble"]
        options.OutputFolder (1,1) string = ""
    end
    [Xp,Xs,Y,metadata] = softsensor.generateSyntheticData(Seed=options.Seed);
    bands = zeros(height(Xs),7);
    for index = 1:7
        columns = (index-1)*20+1:min(index*20,width(Xs));
        if index==7
            columns = (index-1)*20+1:width(Xs);
        end
        bands(:,index) = mean(Xs{:,columns},2);
    end
    rowIndex = (1:height(Y))';
    batchStart = floor((rowIndex-1)/240)*240+1;
    source = max(rowIndex-metadata.residenceTimeMinutes,batchStart);
    % Explicit fixed artificial category vocabulary, with no data-fitted encoding.
    enzymeNames = ["A1","A2","B","C","D","E"];
    mixNames = ["Mix1","Mix2","Mix3"];
    enzyme = double(string(Y.EnzymeCode)==enzymeNames);
    mix = double(string(Y.RawMaterialMix)==mixNames);
    predictors = [bands(source,:),Xp{source,:},enzyme,mix];
    result = softsensor.fitAndEvaluate(predictors,Y{:,options.Target},metadata.batchId, ...
        string(Y.testset)=="1",Group=string(Y.EnzymeCode),Models=options.Models,Seed=options.Seed);
    result.isSynthetic = true;
    result.target = options.Target;
    result.generator = struct("Seed",options.Seed,"Batches",18,"MinutesPerBatch",240);
    result.featureDescription = "Seven fixed spectral bands, 17 process variables and fixed category indicators.";
    result.predictions = addvars(result.predictions,Y.Time(result.predictions.RowIndex), ...
        Before=1,NewVariableNames="Time");
    result.outputFolder = "";
    if options.OutputFolder~=""
        result.outputFolder = softsensor.saveResult(result,options.OutputFolder);
    end
    disp(result.metrics);
end
