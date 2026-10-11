function verification = verifyRelease(options)
% <header> %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%MACA%%
% VERIFYRELEASE Run the release's class-based tests and source Code Analyzer.
% IncludeToolboxTests=true (default) includes Statistics toolbox reference,
% tree/elastic-net, serialization and both-target integration checks. Set false
% for the base-MATLAB suite. Tests use artificial data and temporary directories.
% Returns counts and per-test summaries, without absolute paths or datasets.
% Errors if any test fails/is incomplete or any source analyzer message remains.
% Writes no persistent report; callers can save the returned struct explicitly.
% Example: verification = softsensor.verifyRelease();
% </header> %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%MACA%%
    arguments
        options.IncludeToolboxTests (1,1) logical = true
    end
    root = fileparts(fileparts(mfilename('fullpath')));
    suite = matlab.unittest.TestSuite.fromFolder(fullfile(root,"tests"));
    if ~options.IncludeToolboxTests
        suite = suite.selectIf(~matlab.unittest.selectors.HasTag("Toolbox"));
    end
    results = run(suite);
    summaries = repmat(struct("name","","passed",false,"failed",false, ...
        "incomplete",false,"durationSeconds",0),numel(results),1);
    for index = 1:numel(results)
        summaries(index) = struct("name",string(results(index).Name), ...
            "passed",results(index).Passed,"failed",results(index).Failed, ...
            "incomplete",results(index).Incomplete,"durationSeconds",results(index).Duration);
    end
    sourceFiles = dir(fullfile(root,"+softsensor","*.m"));
    testFiles = dir(fullfile(root,"tests","*.m"));
    files = [sourceFiles;testFiles];
    messages = struct("file",{},"line",{},"message",{});
    for index = 1:numel(files)
        issues = checkcode(fullfile(files(index).folder,files(index).name),"-id");
        for issue = 1:numel(issues)
            messages(end+1) = struct("file",string(files(index).name), ...
                "line",issues(issue).line,"message",string(issues(issue).message)); %#ok<AGROW>
        end
    end
    verification = struct("passed",all([results.Passed]) && isempty(messages), ...
        "matlabVersion",string(version()),"includeToolboxTests",options.IncludeToolboxTests, ...
        "testCount",numel(results),"passedCount",sum([results.Passed]), ...
        "failedCount",sum([results.Failed]),"incompleteCount",sum([results.Incomplete]), ...
        "durationSeconds",sum([results.Duration]),"tests",summaries, ...
        "analyzerFileCount",numel(files),"analyzerMessages",messages);
    disp(verification);
    assert(all([results.Passed]),"SoftSensor:Verification","Release tests failed or were incomplete.");
    assert(isempty(messages),"SoftSensor:Analyzer","Release has Code Analyzer messages.");
end
