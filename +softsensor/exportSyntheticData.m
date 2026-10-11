function outputFolder = exportSyntheticData(Xp,Xs,Y,metadata,parentFolder,options)
% <header> %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%MACA%%
% EXPORTSYNTHETICDATA Explicitly export generated artificial data into a new folder.
% Used by generateSyntheticData; Format="mat" (default), "csv" or "both".
% Never overwrites completed exports. Failure leaves a .partial directory.
% Only pass artificial data here; an isSynthetic flag does not anonymize data.
% </header> %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%MACA%%
    arguments
        Xp timetable
        Xs timetable
        Y timetable
        metadata (1,1) struct
        parentFolder (1,1) string
        options.Format (1,1) string {mustBeMember(options.Format,["mat","csv","both"])} = "mat"
    end
    if ismissing(parentFolder) || strlength(parentFolder)==0
        error("SoftSensor:OutputFolder","Provide a nonempty output parent.");
    end
    if ~isfield(metadata,"isSynthetic") || ~isequal(metadata.isSynthetic,true)
        error("SoftSensor:SyntheticOnly","Synthetic metadata is required.");
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
    outputFolder = string(tempname(attributes.Name));
    staging = outputFolder+".partial";
    [ok,message] = mkdir(staging);
    if ~ok
        error("SoftSensor:OutputFolder","%s",message);
    end
    if options.Format=="mat" || options.Format=="both"
        save(fullfile(staging,"syntheticData.mat"),"Xp","Xs","Y","metadata","-v7");
    end
    if options.Format=="csv" || options.Format=="both"
        writetimetable(Xp,fullfile(staging,"process.csv"));
        writetimetable(Xs,fullfile(staging,"spectra.csv"));
        writetimetable(Y,fullfile(staging,"responses.csv"));
        batches = table(Y.Time,metadata.batchId,string(Y.testset)=="1", ...
            VariableNames={'Time','BatchID','HeldOut'});
        writetable(batches,fullfile(staging,"batches.csv"));
    end
    if isfolder(outputFolder) || isfile(outputFolder)
        error("SoftSensor:OutputExists","Refusing to overwrite an existing export.");
    end
    [ok,message] = movefile(staging,outputFolder);
    if ~ok
        error("SoftSensor:Commit","%s",message);
    end
end
