function [Xp, Xs, Y, metadata] = generateSyntheticData(options)
% <header>  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%MACA%%
% GENERATESYNTHETICDATA Create artificial process, NIR and sparse laboratory timetables.
%--------------------------------------------------------------------------
% Syntax
%   [Xp, Xs, Y, metadata] = softsensor.generateSyntheticData()
%   [Xp, Xs, Y, metadata] = softsensor.generateSyntheticData(Name=Value, ...)
%
% Description
%   Generate hand-written artificial process, spectral and response signals
%   in repeated enzyme/material batches. No original rows, private fitted models
%   or private files are read. Optionally export the synthetic dataset.
%
% Inputs
%   No positional arguments; all settings are name-value arguments.
%
% Name-Value Arguments
%   Seed            - Nonnegative integer RNG seed. Default: 42.
%   Batches         - Integer multiple of six, at least 18. Default: 18.
%   MinutesPerBatch - Integer >= 72; one-minute observations. Default: 240.
%   OutputFolder    - Scalar string export directory. Default: "" (no files).
%   ExportFormat    - "mat" (default), "csv" or "both" for explicit exports.
%
% Outputs
%   Xp       - N-by-17 process timetable, N = Batches*MinutesPerBatch. Artificial
%              flow is L/h, temperature deg C and pressure bar.
%   Xs       - N-by-141 spectral timetable with NIR wavelength column names.
%   Y        - N-by-11 category/laboratory timetable; Mw is in Da. Laboratory
%              responses are sparsely observed; missing measurements are NaN.
%   metadata - Struct with isSynthetic, seed, batchId, wavelengthNm, composition
%              timetable (percent), residenceTimeMinutes=45 and description.
%
% Notes
%   - Row times are shared one-minute observations starting 2020-01-01.
%   - Six enzyme labels and three artificial mixes repeat across batches; final
%     six batches form the marked test set, one for each enzyme label.
%   - Responses use an artificial 45-minute within-batch residence time. Sparse
%     lab sampling starts at minute offset 60; measurements can still be missing.
%   - Includes missing sensor rows/outages. Structural similarity to historical
%     interfaces does not establish distribution or scientific-result equivalence.
%   - Caller RNG is restored on return/error. OutputFolder receives a unique
%     completed export directory; prior exports are never overwritten. MAT is
%     the default format. With the default folder, no files are written.
%
% Example
%   [Xp, Xs, Y, metadata] = softsensor.generateSyntheticData(Seed=42);
%   % Optional export of synthetic files only:
%   softsensor.generateSyntheticData(OutputFolder="my_synthetic_data");
%
% Dependencies
%   Base MATLAB; timetable/categorical, RNG, filtering and optional MAT/CSV I/O.
%
% Metadata
%   Project:               SoftSensorFrameworkPP3
%   Documentation updated: 2026-10-07
%   Status:                Synthetic public entry point; verification passed
%
% See also softsensor.runSyntheticDemo, softsensor.exportSyntheticData.
%--------------------------------------------------------------------------
% </header> %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%MACA%%
    arguments
        options.Seed (1,1) double {mustBeInteger,mustBeNonnegative,mustBeFinite, ...
            mustBeLessThanOrEqual(options.Seed,4294967295)} = 42
        options.Batches (1,1) double {mustBeInteger,mustBePositive,mustBeFinite} = 18
        options.MinutesPerBatch (1,1) double {mustBeInteger,mustBeFinite, ...
            mustBeGreaterThanOrEqual(options.MinutesPerBatch,72)} = 240
        options.OutputFolder (1,1) string = ""
        options.ExportFormat (1,1) string {mustBeMember(options.ExportFormat,["mat","csv","both"])} = "mat"
    end
    if mod(options.Batches,6) ~= 0 || options.Batches < 18
        error("SyntheticData:BatchCount", "Use a multiple of six batches, at least 18.");
    end
    previousRng = rng();
    restoreRng = onCleanup(@() rng(previousRng));
    rng(options.Seed,"twister");

    n = options.Batches*options.MinutesPerBatch;
    time = datetime(2020,1,1) + minutes((0:n-1)');
    batch = repelem((1:options.Batches)',options.MinutesPerBatch);
    group = mod(batch-1,6)+1;
    mix = mod(floor((batch-1)/6) + group-1,3)+1;
    enzymeNames = ["A1","A2","B","C","D","E"];
    enzyme = categorical(enzymeNames(group)',enzymeNames,enzymeNames);
    mixNames = ["Mix1","Mix2","Mix3"];
    rawMix = categorical(mixNames(mix)',mixNames,mixNames);
    slowNoise = filter(1,[1,-0.94],randn(n,4))*sqrt(1-0.94^2);
    batchNoise = randn(options.Batches,4);
    rawState = 2200 + 180*batchNoise(batch,1) + 100*(mix-2);
    waterState = rawState + 600 + 70*batchNoise(batch,2);
    reactorState = rawState + waterState + 300;
    rawFlow = rawState + 65*slowNoise(:,1);
    waterFlow = waterState + 65*slowNoise(:,2);
    reactorFlow = reactorState + 80*slowNoise(:,3);
    rawPercent = 100*rawFlow./(rawFlow + waterFlow + 300);
    pressure = 3.3 + 0.20*batchNoise(batch,3) + 0.13*slowNoise(:,3);
    tt07 = 62 + 0.35*slowNoise(:,1);
    tt08 = 60.7 + 0.45*slowNoise(:,2);
    tt09 = 60.3 + 0.55*slowNoise(:,3);
    tt20 = 59.7 + 0.55*slowNoise(:,4);
    tt12 = 94 + 1.1*slowNoise(:,4);
    dose1 = (1.6 + 0.18*group).*rawState/2200;
    dose2 = double(ismember(group,[4,5]))*0.8.*rawState/2200;
    ratioState = rawState./(rawState+waterState+300);
    process = [pressure,rawFlow,rawState,rawPercent,reactorFlow,reactorState, ...
        tt07,tt08,tt09,tt12,tt20,waterFlow,waterState,ratioState, ...
        hypot(dose1,dose2),dose1,dose2];
    processNames = ["PT03","RawMatFlow","RawMatFlowState","RawMatPercent", ...
        "ReactorFlow","ReactorFlowState","TT07","TT08","TT09","TT12", ...
        "TT20","WaterFlow","WaterFlowState","rawmatRatioState", ...
        "EnzymeNorm","EnzymeLvL1","EnzymeLvL2"];
    Xp = array2timetable(process,"RowTimes",time,"VariableNames",cellstr(processNames));

    protein = 16 + 0.8*(mix-2) + 0.45*slowNoise(:,1);
    fat = 15 + 1.2*batchNoise(batch,4) + 0.65*slowNoise(:,2);
    ash = 5.4 + 0.3*(mix-2) + 0.2*slowNoise(:,3);
    water = 100 - protein - fat - ash;
    wavelength = linspace(1000,2500,141);
    bands = exp(-((wavelength-[1200;1450;1750;2100])/100).^2);
    spectralValues = -1.4 + [fat,water,protein,ash]*bands/65 ...
        + 0.008*randn(n,141) + 0.025*slowNoise(:,4)*ones(1,141);
    Xs = array2timetable(spectralValues,"RowTimes",time, ...
        "VariableNames",cellstr("NIR_"+string(round(wavelength))));

    % The artificial residence time is 45 minutes; never cross batch boundaries.
    offsetInBatch = mod((0:n-1)',options.MinutesPerBatch);
    source = max((1:n)'-45,(batch-1)*options.MinutesPerBatch+1);
    molecularWeight = 6500 + 750*(group-3.5) + 250*(mix-2) ...
        + 500*(pressure(source)-3.3) - 700*(dose1(source)-2.2) ...
        + 220*(protein(source)-16) + 350*randn(n,1);
    smallMolecules = 11 - (molecularWeight-6500)/1900 + 0.35*randn(n,1);
    collagen = 24 + 1.8*(protein(source)-16) + 2.8*randn(n,1);
    brix = 3.5 + 0.015*(rawPercent(source)-40) + 0.25*randn(n,1);
    brixAdjusted = brix./rawPercent(source);
    yield = 33 + 0.8*(protein(source)-16) + 3*randn(n,1);
    yieldWithoutNir = yield + 1.5*randn(n,1);
    responses = [collagen,molecularWeight,smallMolecules,brix,brixAdjusted,yield,yieldWithoutNir];
    observed = offsetInBatch >= 60 & mod(offsetInBatch,12) == 0 & rand(n,1)>0.03;
    responses(~observed,:) = NaN;
    sampleNames = repmat("",n,1);
    sampleNames(observed) = "SYN_" + compose("%06d",find(observed));
    sampleId = categorical(sampleNames);
    sampleId(sampleNames=="") = categorical(missing);
    % Last complete cycle supplies one held-out batch for each enzyme type.
    heldOut = batch > options.Batches-6;
    testSet = categorical(double(heldOut),[0,1],["0","1"]);
    Y = timetable(enzyme,rawMix,sampleId,testSet, ...
        responses(:,1),responses(:,2),responses(:,3),responses(:,4), ...
        responses(:,5),responses(:,6),responses(:,7),RowTimes=time, ...
        VariableNames={'EnzymeCode','RawMaterialMix','product_sampleID','testset', ...
        'Collagen','Mw','SmallMolecules','Brix','BrixAdjusted','Yield','Yield_withoutNIR'});
    % Missingness is sensor-level, plus short acquisition outages.
    Xp{rand(n,1)<0.005,1:14} = NaN;
    Xs{rand(n,1)<0.005,:} = NaN;
    outageStart = options.MinutesPerBatch + 95;
    Xs{outageStart:outageStart+7,:} = NaN;
    batchId = "SYN_BATCH_" + compose("%03d",batch);
    composition = timetable(fat,protein,ash,water,RowTimes=time, ...
        VariableNames={'NIR_Fat','NIR_Protein','NIR_Ash','NIR_Water'});
    metadata = struct("isSynthetic",true,"seed",options.Seed,"batchId",batchId, ...
        "wavelengthNm",wavelength,"composition",composition,"residenceTimeMinutes",45, ...
        "description","Artificial values; structural similarity only; no original records used.");
    if options.OutputFolder ~= ""
        metadata.outputFolder = softsensor.exportSyntheticData(Xp,Xs,Y,metadata, ...
            options.OutputFolder,Format=options.ExportFormat);
    end
end
