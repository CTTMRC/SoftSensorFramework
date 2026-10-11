classdef exportTest < matlab.unittest.TestCase
    % Reproducibility, opt-in exports and failure isolation.
    properties (TestParameter)
        format = {"mat","csv","both"}
    end

    properties
        Parent
    end

    methods (TestClassSetup)
        function addRelease(testCase)
            root = fileparts(fileparts(mfilename('fullpath')));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(root));
        end
    end

    methods (TestMethodSetup)
        function temporaryParent(testCase)
            fixture = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture());
            testCase.Parent = string(fixture.Folder);
        end
    end

    methods (Test)
        function generatorDeterminismAndRNG(testCase)
            before = rng();
            [xp,xs,y,m] = softsensor.generateSyntheticData(MinutesPerBatch=72);
            [xp2,xs2,y2,m2] = softsensor.generateSyntheticData(MinutesPerBatch=72);
            after = rng();
            testCase.verifyEqual(before,after);
            testCase.verifyEqual(xp,xp2);
            testCase.verifyEqual(xs,xs2);
            testCase.verifyEqual(y,y2);
            testCase.verifyEqual(m,m2);
            testCase.verifySize(xp,[1296,17]);
            testCase.verifySize(xs,[1296,141]);
            testCase.verifyEqual(xp.Time,y.Time);
            testCase.verifyFalse(isfield(m,'outputFolder'));
            testCase.verifyTrue(m.isSynthetic);
        end

        function defaultDoesNotWrite(testCase)
            fixture = testCase.applyFixture(matlab.unittest.fixtures.WorkingFolderFixture());
            softsensor.generateSyntheticData(MinutesPerBatch=72);
            result = softsensor.runSyntheticDemo(Models="PLS");
            files = dir(fixture.Folder);
            testCase.verifyEqual(string({files.name}),[".",".."]);
            testCase.verifyEqual(result.outputFolder,"");
        end

        function exportFormats(testCase,format)
            [~,~,~,metadata] = softsensor.generateSyntheticData(MinutesPerBatch=72, ...
                OutputFolder=testCase.Parent,ExportFormat=format);
            files = dir(metadata.outputFolder);
            names = sort(string({files(~[files.isdir]).name}));
            expected = exportTest.expectedFiles(format);
            testCase.verifyEqual(names,expected);
            testCase.verifyEmpty(dir(fullfile(testCase.Parent,"*.partial")));
        end

        function repeatedExportsStayIntact(testCase)
            [~,~,~,first] = softsensor.generateSyntheticData(MinutesPerBatch=72,OutputFolder=testCase.Parent);
            file = fullfile(first.outputFolder,"syntheticData.mat");
            firstBytes = exportTest.readBytes(file);
            [~,~,~,second] = softsensor.generateSyntheticData(MinutesPerBatch=72,OutputFolder=testCase.Parent);
            testCase.verifyNotEqual(first.outputFolder,second.outputFolder);
            testCase.verifyEqual(firstBytes,exportTest.readBytes(file));
            loaded = load(file,"metadata","Xp");
            testCase.verifyTrue(loaded.metadata.isSynthetic);
            testCase.verifySize(loaded.Xp,[1296,17]);
        end

        function rejectsNonSyntheticMetadata(testCase)
            [xp,xs,y,m] = softsensor.generateSyntheticData(MinutesPerBatch=72);
            m.isSynthetic = false;
            testCase.verifyError(@() softsensor.exportSyntheticData(xp,xs,y,m,testCase.Parent), ...
                "SoftSensor:SyntheticOnly");
        end

        function failedExportIsPartial(testCase)
            [xp,xs,y,m] = softsensor.generateSyntheticData(MinutesPerBatch=72);
            m.batchId = m.batchId(1:end-1);
            testCase.verifyError(@() softsensor.exportSyntheticData(xp,xs,y,m,testCase.Parent, ...
                Format="csv"),"MATLAB:table:UnequalVarLengths");
            folders = dir(testCase.Parent);
            testCase.verifyEqual(numel(folders),3);
            testCase.verifyEqual(numel(dir(fullfile(testCase.Parent,"*.partial"))),1);
            testCase.verifyEmpty(dir(fullfile(testCase.Parent,"**","COMPLETE.json")));
        end

        function repeatedResultsStayIntact(testCase)
            result = exportTest.smallResult();
            first = softsensor.saveResult(result,testCase.Parent);
            bytes = exportTest.readBytes(fullfile(first,"result.mat"));
            second = softsensor.saveResult(result,testCase.Parent);
            testCase.verifyNotEqual(first,second);
            testCase.verifyEqual(bytes,exportTest.readBytes(fullfile(first,"result.mat")));
            loaded = load(fullfile(first,"result.mat"),"result");
            testCase.verifyEqual(loaded.result.outputFolder,first);
            record = jsondecode(fileread(fullfile(first,"COMPLETE.json")));
            testCase.verifyTrue(record.complete);
            testCase.verifyFalse(isfile(fullfile(first,"syntheticData.mat")));
            testCase.verifyEmpty(dir(fullfile(testCase.Parent,"*.partial")));
        end

        function failedResultHasNoCompletion(testCase)
            result = exportTest.smallResult();
            result.tuning = 7;
            testCase.verifyError(@() softsensor.saveResult(result,testCase.Parent), ...
                "MATLAB:table:write:UnsupportedInputType");
            testCase.verifyEqual(numel(dir(fullfile(testCase.Parent,"*.partial"))),1);
            testCase.verifyEmpty(dir(fullfile(testCase.Parent,"**","COMPLETE.json")));
        end

        function invalidOutputParents(testCase)
            result = exportTest.smallResult();
            testCase.verifyError(@() softsensor.saveResult(result,""),"SoftSensor:OutputFolder");
            testCase.verifyError(@() softsensor.saveResult(struct(),testCase.Parent),"SoftSensor:Result");
            occupied = fullfile(testCase.Parent,"file.txt");
            writelines("occupied",occupied);
            testCase.verifyError(@() softsensor.saveResult(result,occupied),"SoftSensor:OutputFolder");
        end

        function invalidGeneratorArguments(testCase)
            testCase.verifyError(@() softsensor.generateSyntheticData(Batches=19),"SyntheticData:BatchCount");
            testCase.verifyError(@() softsensor.generateSyntheticData(MinutesPerBatch=60), ...
                "MATLAB:validators:mustBeGreaterThanOrEqual");
            testCase.verifyError(@() softsensor.generateSyntheticData(Seed=4294967296), ...
                "MATLAB:validators:mustBeLessThanOrEqual");
        end
    end

    methods (Static, Access = private)
        function names = expectedFiles(format)
            switch format
                case "mat"
                    names = "syntheticData.mat";
                case "csv"
                    names = ["batches.csv","process.csv","responses.csv","spectra.csv"];
                case "both"
                    names = ["batches.csv","process.csv","responses.csv","spectra.csv","syntheticData.mat"];
                otherwise
                    error("Test:Format","Unsupported fixture format.");
            end
        end

        function result = smallResult()
            result = softsensor.fitAndEvaluate((1:12)',2*(1:12)'+1, ...
                repelem(["a";"b";"c";"test"],3),[false(9,1);true(3,1)],Models="PLS");
        end

        function bytes = readBytes(file)
            handle = fopen(file,"rb");
            closer = onCleanup(@() fclose(handle));
            bytes = fread(handle,Inf,"*uint8");
        end
    end
end
