classdef customPLSTest < matlab.unittest.TestCase
    % Numerical regressions and isolation checks using artificial data only.
    properties (TestParameter)
        component = {1,2,3,4}
        scale = {false,true}
        reduce = {false,true}
    end

    methods (TestClassSetup)
        function addRelease(testCase)
            root = fileparts(fileparts(mfilename('fullpath')));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(root));
        end
    end

    methods (TestMethodSetup)
        function seedRandom(testCase)
            previous = rng();
            testCase.addTeardown(@() rng(previous));
            rng(17,"twister");
        end
    end

    methods (Test)
        function exactRegression(testCase)
            X = (-5:5)';
            model = softsensor.customPLS(X,3*X+4,1);
            testCase.verifyEqual(model.beta,[4;3],AbsTol=1e-12);
            testCase.verifyEqual(softsensor.predictCustomPLS(model,[11;12]),[37;40],AbsTol=1e-12);
        end

        function fullRankMatchesOLS(testCase)
            X = randn(50,6);
            y = randn(50,1);
            model = softsensor.customPLS(X,y,6);
            testCase.verifyEqual(model.beta,[ones(50,1),X]\y,AbsTol=1e-10);
        end

        function rankExhaustion(testCase)
            x = (-5:5)';
            model = softsensor.customPLS([x,x,ones(11,1)],3*x+4,10,Scale=true);
            testCase.verifyEqual(model.nComponents,1);
            testCase.verifyEqual(model.beta,[4;1.5;1.5;0],AbsTol=1e-12);
            testCase.verifyEqual(model.vip(3),0,AbsTol=1e-12);
        end

        function constantResponse(testCase)
            model = softsensor.customPLS(randn(20,3),7*ones(20,1),3);
            testCase.verifyEqual(model.nComponents,0);
            testCase.verifyEqual(softsensor.predictCustomPLS(model,randn(2,3)),[7;7],AbsTol=1e-12);
        end

        function constantPredictors(testCase)
            model = softsensor.customPLS(ones(10,3),(1:10)',3,Scale=true,Blocks=[1,1,2]);
            testCase.verifyEqual(model.nComponents,0);
            testCase.verifyEqual(softsensor.predictCustomPLS(model,ones(2,3)),[5.5;5.5],AbsTol=1e-12);
        end

        function scoreProjection(testCase,component,scale)
            X = randn(50,6);
            y = X(:,1)+randn(50,1);
            model = softsensor.customPLS(X,y,component,Scale=scale,Blocks=[1,1,1,2,2,2]);
            [~,scores] = softsensor.predictCustomPLS(model,X);
            testCase.verifyEqual(scores,model.xScores,AbsTol=1e-10);
            testCase.verifyEqual(sum(model.vip.^2),6,AbsTol=1e-10);
        end

        function cvPrefixMatchesSeparateFit(testCase,component,scale)
            X = randn(40,6);
            y = X(:,1)-2*X(:,2)+randn(40,1);
            folds = repelem((1:4)',10);
            tuned = softsensor.customPLSCV(X,y,folds,MaxComponents=4,Scale=scale,Blocks=[1,1,1,2,2,2]);
            independent = softsensor.customPLS(X(folds~=1,:),y(folds~=1),component, ...
                Scale=scale,Blocks=[1,1,1,2,2,2]);
            testCase.verifyEqual(tuned.cv.predictionPath(folds==1,component), ...
                softsensor.predictCustomPLS(independent,X(folds==1,:)),AbsTol=1e-10);
        end

        function validationLabelsAreIsolated(testCase,reduce,scale)
            X = randn(40,6);
            y = X(:,1)-2*X(:,2)+randn(40,1);
            folds = repelem((1:4)',10);
            changedY = y;
            changedY(folds==1) = 1000+changedY(folds==1);
            original = softsensor.customPLSCV(X,y,folds,MaxComponents=4,Scale=scale, ...
                Reduce=reduce,Blocks=[1,1,1,2,2,2]);
            changed = softsensor.customPLSCV(X,changedY,folds,MaxComponents=4,Scale=scale, ...
                Reduce=reduce,Blocks=[1,1,1,2,2,2]);
            testCase.verifyEqual(original.cv.predictionPath(folds==1,:), ...
                changed.cv.predictionPath(folds==1,:),AbsTol=1e-12);
            testCase.verifySize(original.predictorMask,[6,1]);
            testCase.verifyTrue(all(isfinite(softsensor.predictCustomPLS(original,X))));
        end

        function reducedPredictorsStayAligned(testCase)
            X = randn(40,8);
            X(:,8) = 1;
            y = 3*X(:,1)+0.01*randn(40,1);
            model = softsensor.customPLSCV(X,y,repelem((1:4)',10),MaxComponents=3,Reduce=true);
            testCase.verifyEqual(model.nPredictors,8);
            testCase.verifyFalse(model.predictorMask(8));
            testCase.verifyEqual(model.betaPath([false;~model.predictorMask],:), ...
                zeros(sum(~model.predictorMask),size(model.betaPath,2)),AbsTol=1e-12);
        end

        function incompleteAndEmptyInference(testCase)
            model = softsensor.customPLS((1:5)',(2:6)',1);
            testCase.verifyEqual(softsensor.predictCustomPLS(model,[NaN;Inf;6]),[NaN;NaN;7],AbsTol=1e-12);
            [prediction,scores] = softsensor.predictCustomPLS(model,zeros(0,1));
            testCase.verifySize(prediction,[0,1]);
            testCase.verifySize(scores,[0,1]);
            testCase.verifyEqual(softsensor.predictCustomPLS(model,6,0),4,AbsTol=1e-12);
            testCase.verifyEqual(softsensor.predictCustomPLS(model,6,100),7,AbsTol=1e-12);
        end

        function invalidShapesAndFolds(testCase)
            X = randn(4,3);
            y = (1:4)';
            testCase.verifyError(@() softsensor.customPLS(X,y(1:2),2),"CustomPLS:Shape");
            testCase.verifyError(@() softsensor.customPLSCV(X,y,[1;1;1;2]),"CustomPLS:FoldTraining");
            testCase.verifyError(@() softsensor.customPLSCV(X,y,ones(4,1)),"CustomPLS:Folds");
            testCase.verifyError(@() softsensor.customPLSCV(X,y,[1;1;2;2],Blocks=1),"CustomPLS:Blocks");
            model = softsensor.customPLS(X,y,2);
            testCase.verifyError(@() softsensor.predictCustomPLS(model,X(:,1:2)),"CustomPLS:Predictors");
        end

        function rowSNV(testCase)
            transformed = softsensor.spectralSNV([1,2,3;4,4,4;NaN,1,2]);
            testCase.verifyEqual(transformed(1,:),[-1,0,1],AbsTol=1e-12);
            testCase.verifyEqual(transformed(2,:),[0,0,0],AbsTol=1e-12);
            testCase.verifyTrue(all(isnan(transformed(3,:))));
            testCase.verifySize(softsensor.spectralSNV(zeros(0,3)),[0,3]);
        end

        function publicVerification(testCase)
            verification = softsensor.verifyCustomPLS();
            testCase.verifyTrue(verification.passed);
        end
    end

    methods (Test, TestTags = {'Toolbox'})
        function matchesPLSReference(testCase,component,scale)
            X = randn(50,6);
            y = X(:,1)-2*X(:,2)+0.2*randn(50,1);
            newX = randn(10,6);
            model = softsensor.customPLS(X,y,component,Scale=scale);
            referenceX = (X-model.xMean)./model.xScale;
            [~,~,~,~,beta] = plsregress(referenceX,y,component);
            expected = [ones(10,1),(newX-model.xMean)./model.xScale]*beta;
            testCase.verifyEqual(softsensor.predictCustomPLS(model,newX),expected,AbsTol=1e-9);
        end
    end
end
