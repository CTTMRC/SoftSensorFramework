function outputFolder = saveResult(result,parentFolder)
% <header> %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%MACA%%
% SAVERESULT Save a fitted result once in a new run folder under parentFolder.
% Writes result.mat plus metrics/predictions/tuning CSVs and COMPLETE.json.
% A same-parent staging directory is renamed only after all writes succeed.
% Failure leaves an uncommitted .partial folder for inspection. Never overwrites
% a completed run. Caller-supplied results can contain private data: save them
% outside the source release. This saves the supplied struct in full;
% fitAndEvaluate returns compact models without raw calibration matrices.
% </header> %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%MACA%%
    arguments
        result (1,1) struct
        parentFolder (1,1) string
    end
    if ismissing(parentFolder) || strlength(parentFolder)==0
        error("SoftSensor:OutputFolder","Provide a nonempty output parent.");
    end
    required = ["models","metrics","predictions","tuning"];
    if ~all(isfield(result,required))
        error("SoftSensor:Result","Result must contain models, metrics, predictions and tuning.");
    end
    if ~isfolder(parentFolder)
        [ok,message] = mkdir(parentFolder);
        if ~ok
            error("SoftSensor:OutputFolder","%s",message);
        end
    end
    [ok,attributes] = fileattrib(parentFolder);
    if ~ok
        error("SoftSensor:OutputFolder","Cannot resolve output parent.");
    end
    base = string(tempname(attributes.Name));
    staging = base+".partial";
    outputFolder = base;
    [ok,message] = mkdir(staging);
    if ~ok
        error("SoftSensor:OutputFolder","%s",message);
    end
    result.outputFolder = outputFolder;
    save(fullfile(staging,"result.mat"),"result","-v7");
    writetable(result.metrics,fullfile(staging,"metrics.csv"));
    writetable(result.predictions,fullfile(staging,"predictions.csv"));
    writetable(result.tuning,fullfile(staging,"tuning.csv"));
    summary = struct("complete",true,"schemaVersion",1,"matlabVersion",version(), ...
        "files",["result.mat","metrics.csv","predictions.csv","tuning.csv"]);
    writeText(fullfile(staging,"COMPLETE.json"),jsonencode(summary,PrettyPrint=true));
    if isfolder(outputFolder) || isfile(outputFolder)
        error("SoftSensor:OutputExists","Refusing to overwrite an existing run.");
    end
    [ok,message] = movefile(staging,outputFolder);
    if ~ok
        error("SoftSensor:Commit","%s",message);
    end
end

function writeText(file,text)
    [handle,message] = fopen(file,"w","n","UTF-8");
    if handle<0
        error("SoftSensor:Write","%s",message);
    end
    closer = onCleanup(@() fclose(handle));
    count = fprintf(handle,"%s\n",text);
    if count<strlength(text)
        error("SoftSensor:Write","Incomplete completion-record write.");
    end
end
