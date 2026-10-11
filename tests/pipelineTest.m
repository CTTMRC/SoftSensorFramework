classdef pipelineTest < matlab.unittest.TestCase
    % Public fitting, metrics, input contracts and saved inference.
    properties (TestParameter)
        family = {"PLS","Groupwise PLS","ElasticNet","Tree Ensemble"}
        target = {"Mw","SmallMolecules"}
    end

    properties
        X
        Y
        Batches
        HeldOut
        Groups
    end

    methods (TestClassSetup)
        function addRelease(testCase)
            root = fileparts(fileparts(mfilename('fullpath')));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(root));
        end
    end

    methods (TestMethodSetup)
        function artificialFixture(testCase)
            previous = rng();
            testCase.addTeardown(@() rng(previous));
            rng(31,"twister");
            testCase.X = [randn(60,3),ones(60,1)];
            testCase.Y = testCase.X*[3;-2;0.5;0]+2+0.1*randn(60,1);
            testCase.Batches = repelem("batch "+string((1:10)'),6);
            testCase.HeldOut = ismember(testCase.Batches,["batch 5","batch 10"]);
            testCase.Groups = repelem(["first group";"group / two"],30);
        end
    end

    methods (Test)
        function plsAndArbitraryGroups(testCase)
            r = softsensor.fitAndEvaluate(testCase.X,testCase.Y,testCase.Batches,testCase.HeldOut, ...
                Group=testCase.Groups,Models=["Groupwise PLS","PLS"]);
            testCase.verifyEqual(r.metrics.Model,["Groupwise PLS";"PLS"]);
            testCase.verifyEqual(r.predictions.Properties.VariableNames(end-1:end),{'GroupwisePLS','PLS'});
            testCase.verifyEmpty(intersect(r.calibrationBatches,r.testBatches));
            predicted = softsensor.predictFittedModel(r.models{1},testCase.X(testCase.HeldOut,:), ...
                Group=testCase.Groups(testCase.HeldOut));
            testCase.verifyEqual(predicted,r.predictions.GroupwisePLS,AbsTol=1e-12);
            testCase.verifyFalse(isfield(r.settings,'Group'));
            testCase.verifyFalse(isfield(r.models{2}.estimator,'xScores'));
        end

        function testLabelsAreIsolated(testCase)
            changedY = testCase.Y;
            changedY(testCase.HeldOut) = changedY(testCase.HeldOut)+1000;
            original = softsensor.fitAndEvaluate(testCase.X,testCase.Y,testCase.Batches, ...
                testCase.HeldOut,Models="PLS");
            changed = softsensor.fitAndEvaluate(testCase.X,changedY,testCase.Batches, ...
                testCase.HeldOut,Models="PLS");
            testCase.verifyEqual(original.predictions.PLS,changed.predictions.PLS,AbsTol=1e-12);
            testCase.verifyEqual(original.tuning.BatchCV_RMSE,changed.tuning.BatchCV_RMSE,AbsTol=1e-12);
        end

        function droppedRowsRetainAlignment(testCase)
            predictors = testCase.X;
            y = testCase.Y;
            predictors(2,1) = NaN;
            y(25) = Inf;
            predictors(29,2) = Inf;
            r = softsensor.fitAndEvaluate(predictors,y,testCase.Batches,testCase.HeldOut,Models="PLS");
            testCase.verifyEqual(r.droppedRows,[2;25;29]);
            testCase.verifyEqual(r.completeRows,57);
            expectedRows = find(testCase.HeldOut & isfinite(y) & all(isfinite(predictors),2));
            testCase.verifyEqual(r.predictions.RowIndex,expectedRows);
            testCase.verifyEqual(r.predictions.Measured,testCase.Y(expectedRows),AbsTol=1e-12);
        end

        function rejectsOverlapBeforeDropping(testCase)
            heldOut = testCase.HeldOut;
            heldOut(1) = true;
            predictors = testCase.X;
            predictors(1,:) = NaN;
            testCase.verifyError(@() softsensor.fitAndEvaluate(predictors,testCase.Y,testCase.Batches, ...
                heldOut,Models="PLS"),"SoftSensor:BatchOverlap");
        end

        function rejectsUnseenGroups(testCase)
            groups = testCase.Groups;
            groups(testCase.HeldOut) = "unseen";
            testCase.verifyError(@() softsensor.fitAndEvaluate(testCase.X,testCase.Y,testCase.Batches, ...
                testCase.HeldOut,Group=groups,Models="Groupwise PLS"),"SoftSensor:UnseenGroup");
            r = softsensor.fitAndEvaluate(testCase.X,testCase.Y,testCase.Batches,testCase.HeldOut, ...
                Group=testCase.Groups,Models="Groupwise PLS");
            testCase.verifyError(@() softsensor.predictFittedModel(r.models{1},testCase.X(1,:), ...
                Group="unseen"),"SoftSensor:UnseenGroup");
            testCase.verifyError(@() softsensor.predictFittedModel(r.models{1},testCase.X(1,:)), ...
                "SoftSensor:Groups");
            testCase.verifySize(softsensor.predictFittedModel(r.models{1},zeros(0,4)),[0,1]);
        end

        function rejectsInvalidInputs(testCase)
            testCase.verifyError(@() softsensor.fitAndEvaluate(testCase.X,testCase.Y(1:2), ...
                testCase.Batches,testCase.HeldOut,Models="PLS"),"SoftSensor:Shape");
            testCase.verifyError(@() softsensor.fitAndEvaluate(testCase.X,testCase.Y, ...
                testCase.Batches,2*double(testCase.HeldOut),Models="PLS"),"SoftSensor:SplitMask");
            testCase.verifyError(@() softsensor.fitAndEvaluate(testCase.X,testCase.Y, ...
                testCase.Batches,testCase.HeldOut,Models=["PLS","PLS"]),"SoftSensor:Models");
            testCase.verifyError(@() softsensor.fitAndEvaluate(testCase.X,testCase.Y, ...
                testCase.Batches,testCase.HeldOut,Models="Groupwise PLS"),"SoftSensor:Groups");
            testCase.verifyError(@() softsensor.fitAndEvaluate(testCase.X,testCase.Y, ...
                testCase.Batches,false(60,1),Models="PLS"),"SoftSensor:Split");
        end

        function missingIdentifiers(testCase)
            batches = testCase.Batches;
            batches(1) = missing;
            groups = testCase.Groups;
            groups(1) = "";
            testCase.verifyError(@() softsensor.fitAndEvaluate(testCase.X,testCase.Y,batches, ...
                testCase.HeldOut,Models="PLS"),"SoftSensor:BatchID");
            testCase.verifyError(@() softsensor.fitAndEvaluate(testCase.X,testCase.Y,testCase.Batches, ...
                testCase.HeldOut,Group=groups,Models="PLS"),"SoftSensor:Groups");
        end

        function insufficientGroupBatches(testCase)
            groups = testCase.Groups;
            groups(testCase.Batches=="batch 1") = "singleton";
            testCase.verifyError(@() softsensor.fitAndEvaluate(testCase.X,testCase.Y,testCase.Batches, ...
                testCase.HeldOut,Group=groups,Models="Groupwise PLS"),"SoftSensor:Folds");
            testCase.verifyError(@() softsensor.fitAndEvaluate([1;2;3;4],[1;2;3;4], ...
                ["a";"a";"b";"c"],[false;false;false;true],Models="PLS"),"SoftSensor:FoldTraining");
        end

        function metricsFinitePairsAndEdges(testCase)
            metrics = softsensor.evaluateMetrics([1,2,NaN,Inf],[1;0;7;3]);
            testCase.verifyEqual(metrics{1,:},[sqrt(2),1,2/3,-7],AbsTol=1e-12);
            constant = softsensor.evaluateMetrics([0;0],[1;1]);
            testCase.verifyTrue(isnan(constant.WAPE) && isnan(constant.R2Y));
            emptyMetrics = softsensor.evaluateMetrics([],[]);
            testCase.verifyTrue(all(isnan(emptyMetrics{1,:})));
            testCase.verifyError(@() softsensor.evaluateMetrics(ones(2),ones(2)),"Metrics:Shape");
        end

        function rejectsUnknownSavedModel(testCase)
            model = struct('name',"unknown",'nPredictors',4);
            testCase.verifyError(@() softsensor.predictFittedModel(model,zeros(0,4)),"SoftSensor:Model");
        end
    end

    methods (Test, TestTags = {'Toolbox'})
        function savedFamilyInference(testCase,family)
            fixture = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture());
            parent = fixture.Folder;
            r = softsensor.fitAndEvaluate(testCase.X,testCase.Y,testCase.Batches,testCase.HeldOut, ...
                Group=testCase.Groups,Models=family);
            folder = softsensor.saveResult(r,parent);
            loaded = load(fullfile(folder,"result.mat"),"result");
            predicted = softsensor.predictFittedModel(loaded.result.models{1}, ...
                testCase.X(testCase.HeldOut,:),Group=testCase.Groups(testCase.HeldOut));
            testCase.verifyEqual(predicted,r.predictions{:,end},AbsTol=1e-12);
            testCase.verifySize(softsensor.predictFittedModel(loaded.result.models{1},zeros(0,4)),[0,1]);
            incomplete = testCase.X(1,:);
            incomplete(1) = NaN;
            testCase.verifyTrue(isnan(softsensor.predictFittedModel(loaded.result.models{1},incomplete, ...
                Group=testCase.Groups(1))));
        end

        function constantsAllFamilies(testCase)
            r = softsensor.fitAndEvaluate(ones(60,4),testCase.Y,testCase.Batches,testCase.HeldOut, ...
                Group=testCase.Groups);
            expected = mean(testCase.Y(~testCase.HeldOut));
            testCase.verifyEqual(r.predictions{:,{'PLS','ElasticNet','TreeEnsemble'}}, ...
                repmat(expected,12,3),AbsTol=1e-12);
            constantY = softsensor.fitAndEvaluate(testCase.X,7*ones(60,1), ...
                testCase.Batches,testCase.HeldOut,Group=testCase.Groups);
            testCase.verifyEqual(constantY.predictions{:,end-3:end},7*ones(12,4),AbsTol=1e-12);
        end

        function syntheticTargetsDeterministic(testCase,target)
            previous = rng();
            first = softsensor.runSyntheticDemo(Target=target);
            after = rng();
            second = softsensor.runSyntheticDemo(Target=target);
            testCase.verifyEqual(after,previous);
            testCase.verifyEqual(first.predictions,second.predictions);
            testCase.verifyEqual(first.tuning,second.tuning);
            testCase.verifyTrue(all(isfinite(first.predictions{:,end-3:end}),"all"));
            testCase.verifyEqual(first.outputFolder,"");
            testCase.verifyEqual(first.modelNames,["PLS","Groupwise PLS","ElasticNet","Tree Ensemble"]);
        end
    end
end
