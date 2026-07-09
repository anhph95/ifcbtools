function export_ifcb(inputMatFile, options)
%EXPORT_IFCB Export one IFCB summary MAT product to processing CSV files.
%
% export_ifcb("path/to/product.mat")
%
% The one-argument form scans available table variables, asks you to confirm
% or select data_table, infers data_type when possible, and asks you to
% confirm count/carbon before export.
%
% Scripted carbon example:
%   export_ifcb( ...
%       "\\sosiknas1\IFCB_products\NESLTER_transect\summary\carbon_group_class_withTS.mat", ...
%       data_table="classC_opt_adhoc_merge", ...
%       data_type="carbon")
%
% Scripted count example:
%   export_ifcb( ...
%       "\\sosiknas1\IFCB_products\NESLTER_transect\summary\count_group_class_withTS.mat", ...
%       data_table="classcount_opt_adhoc_merge", ...
%       data_type="count")
%
% Outputs:
%   <output_dir>/ifcb_class.csv
%   <output_dir>/ifcb_count.csv or <output_dir>/ifcb_carbon.csv
%
% Scientific data model:
%   rows    = IFCB samples or sample-derived records
%   columns = metadata fields plus taxon/class biomass/count variables

arguments
    inputMatFile (1, 1) string
    options.data_table (1, 1) string = ""
    options.data_type (1, 1) string = ""
    options.output_dir (1, 1) string = ""
end

%% --- RESOLVE INPUT AND OUTPUT PATHS ---

dataTableName = options.data_table;
dataType = lower(options.data_type);

if strlength(dataType) > 0 && ~ismember(dataType, ["count", "carbon"])
    error("IFCB:InvalidDataType", "Provide data_type=""count"" or data_type=""carbon"".");
end

outputDir = options.output_dir;

if ~isfile(inputMatFile)
    error("IFCB:MissingInput", "Input MAT file does not exist: %s", inputMatFile);
end

%% --- CONFIGURE WORKFLOW LOGGING ---

logDir = fullfile(pwd, 'logs');
if ~exist(logDir, 'dir')
    mkdir(logDir);
end
workflowName = "ifcb_export_MATLAB";
logTimestamp = string(datetime("now", "Format", "yyyyMMdd_HHmmss"));
outLogPath = fullfile(logDir, workflowName + "_" + logTimestamp + ".out.log");
errLogPath = fullfile(logDir, workflowName + "_" + logTimestamp + ".err.log");
fclose(fopen(errLogPath, "a"));
diary(outLogPath);
fprintf("%s | INFO | %s | Logging to: %s\n", ...
    string(datetime("now", "Format", "yyyy-MM-dd HH:mm:ss")), workflowName, outLogPath);
fprintf("%s | INFO | %s | Errors to: %s\n", ...
    string(datetime("now", "Format", "yyyy-MM-dd HH:mm:ss")), workflowName, errLogPath);
fprintf("Run configuration:\n");
fprintf("  working_directory: %s\n", pwd);
fprintf("  input_mat_file: %s\n", inputMatFile);
fprintf("  output_dir: %s\n", outputDir);
fprintf("  data_table: %s\n", dataTableName);
fprintf("  data_type: %s\n", dataType);
fprintf("  log_dir: %s\n", logDir);

try
    fprintf('\nLoading file: %s\n', inputMatFile);
    loadedData = load(inputMatFile);

    if strlength(dataTableName) == 0
        dataTableName = choose_data_table(loadedData);
    end

    if strlength(dataType) == 0
        dataType = choose_data_type(inputMatFile, dataTableName);
    end

    if strlength(outputDir) == 0
        outputDir = choose_output_dir();
    end
    if ~exist(outputDir, 'dir')
        fprintf('Creating output directory: %s\n', outputDir);
        mkdir(outputDir);
    end

    fprintf('Resolved data_table: %s\n', dataTableName);
    fprintf('Resolved data_type: %s\n', dataType);
    fprintf('Resolved output_dir: %s\n', outputDir);

    %% --- REQUIRE AND PREPARE METADATA ---

    if ~isfield(loadedData, 'meta_data') || ~istable(loadedData.meta_data)
        error("IFCB:MissingMetadata", "Input MAT file must contain table variable 'meta_data'.");
    end

    metaData = loadedData.meta_data;
    if ~ismember("pid", string(metaData.Properties.VariableNames))
        error("IFCB:MissingPid", "meta_data must contain a 'pid' column.");
    end
    if numel(unique(metaData.pid)) ~= height(metaData)
        error("IFCB:DuplicateMetadataPid", "meta_data contains duplicate pid values.");
    end

    if ismember("sample_time", string(metaData.Properties.VariableNames))
        rawTimes = metaData.sample_time;
        if iscellstr(rawTimes) || isstring(rawTimes)
            trimmedTimes = erase(string(rawTimes), '+00:00');
            metaData.sample_time = datetime(trimmedTimes, ...
                'InputFormat', 'yyyy-MM-dd HH:mm:ss', ...
                'TimeZone','UTC');
        end
    end

    %% --- REQUIRE AND EXPORT CLASS LABELS ---

    loadedNames = string(fieldnames(loadedData));
    labelNames = loadedNames(endsWith(loadedNames, "_label"));
    if isempty(labelNames)
        error("IFCB:MissingClassLabels", "Input MAT file must contain one or more '*_label' variables.");
    end

    classTable = table();
    fprintf('Extracting class labels...\n');
    for i = 1:numel(labelNames)
        labelName = labelNames(i);
        labelData = loadedData.(char(labelName));

        if isstring(labelData)
            labelData = cellstr(labelData);
        elseif ischar(labelData)
            labelData = cellstr(labelData);
        elseif ~iscell(labelData)
            error("IFCB:InvalidClassLabels", "Class label variable '%s' must be a cell, string, or char array.", labelName);
        end

        labelData = labelData(:);
        groupName = erase(labelName, "_label");
        groupLabels = repmat(cellstr(groupName), size(labelData, 1), 1);
        tmp = table(labelData, groupLabels, 'VariableNames', {'class', 'label'});
        classTable = [classTable; tmp]; %#ok<AGROW>
        fprintf('  Processed label group: %s (%d entries)\n', groupName, size(labelData, 1));
    end

    classOutFile = fullfile(outputDir, 'ifcb_class.csv');
    fprintf('Writing IFCB class table: %s\n', classOutFile);
    writetable(classTable, classOutFile);

    %% --- REQUIRE ONE MEASUREMENT PRODUCT ---

    if ~isfield(loadedData, dataTableName) || ~istable(loadedData.(char(dataTableName)))
        error("IFCB:MissingMeasurementTable", ...
            "Input MAT file must contain table variable '%s'.", dataTableName);
    end

    dataTable = loadedData.(char(dataTableName));
    if dataType == "count"
        outFile = fullfile(outputDir, 'ifcb_count.csv');
    else
        outFile = fullfile(outputDir, 'ifcb_carbon.csv');
    end
    fprintf('Using measurement table: %s (%s)\n', dataTableName, dataType);

    if ~ismember("pid", string(dataTable.Properties.VariableNames))
        error("IFCB:MissingPid", "The %s measurement table must contain a 'pid' column.", dataType);
    end

    [foundPid, metaIndex] = ismember(dataTable.pid, metaData.pid);
    if any(~foundPid)
        missingPid = string(dataTable.pid(find(~foundPid, 1, 'first')));
        error("IFCB:UnmatchedPid", "Measurement pid has no matching meta_data row. First missing pid: %s", missingPid);
    end

    matchedMeta = metaData(metaIndex, :);
    metadataNames = string(matchedMeta.Properties.VariableNames);
    dataNames = string(dataTable.Properties.VariableNames);
    measurementNames = dataNames(~ismember(dataNames, metadataNames));
    mergedTable = [matchedMeta dataTable(:, cellstr(measurementNames))];

    fprintf('Writing merged IFCB %s table: %s\n', dataType, outFile);
    writetable(mergedTable, outFile);

    fprintf('\nIFCB %s product exported successfully.\n', dataType);
    diary off
catch ME
    diary off
    fid = fopen(errLogPath, "a");
    fprintf(fid, "%s | ERROR | %s | %s\n", ...
        string(datetime("now", "Format", "yyyy-MM-dd HH:mm:ss")), ...
        workflowName, getReport(ME, "extended", "hyperlinks", "off"));
    fclose(fid);
    rethrow(ME)
end
end

function dataTableName = choose_data_table(loadedData)
loadedNames = string(fieldnames(loadedData));
tableNames = strings(0);
for i = 1:numel(loadedNames)
    thisName = loadedNames(i);
    if thisName ~= "meta_data" && istable(loadedData.(char(thisName)))
        tableNames(end + 1) = thisName; %#ok<AGROW>
    end
end

if isempty(tableNames)
    error("IFCB:MissingDataTable", ...
        "Input MAT file contains no table variables besides 'meta_data'.");
end

if numel(tableNames) == 1
    thisTable = loadedData.(char(tableNames(1)));
    fprintf('\nDetected measurement table variable:\n');
    fprintf('  %s (%d rows, %d columns)\n', tableNames(1), height(thisTable), width(thisTable));
    answer = lower(strtrim(string(input(sprintf('Use data_table "%s"? [Y/n]: ', tableNames(1)), 's'))));
    if answer == "" || answer == "y" || answer == "yes"
        dataTableName = tableNames(1);
        return
    end
else
    fprintf('\nAvailable measurement table variables:\n');
    for i = 1:numel(tableNames)
        thisTable = loadedData.(char(tableNames(i)));
        fprintf('  %d. %s (%d rows, %d columns)\n', ...
            i, tableNames(i), height(thisTable), width(thisTable));
    end
end

choice = strtrim(string(input('Select data table number or enter data_table name: ', 's')));
choiceNumber = str2double(choice);
if ~isnan(choiceNumber)
    if choiceNumber < 1 || choiceNumber > numel(tableNames) || choiceNumber ~= fix(choiceNumber)
        error("IFCB:InvalidDataTableChoice", "Selected data table number is invalid.");
    end
    dataTableName = tableNames(choiceNumber);
else
    dataTableName = choice;
end
end

function dataType = choose_data_type(inputMatFile, dataTableName)
inferredType = infer_data_type(inputMatFile, dataTableName);
if strlength(inferredType) > 0
    answer = lower(strtrim(string(input(sprintf('Use data_type "%s"? [Y/n]: ', inferredType), 's'))));
    if answer == "" || answer == "y" || answer == "yes"
        dataType = inferredType;
        return
    end
end

fprintf('\nData type options:\n');
fprintf('  1. count\n');
fprintf('  2. carbon\n');
dataType = "";
while strlength(dataType) == 0
    typeChoice = lower(strtrim(string(input('Select data_type [count/carbon or 1/2]: ', 's'))));
    if typeChoice == "1" || typeChoice == "count"
        dataType = "count";
    elseif typeChoice == "2" || typeChoice == "carbon"
        dataType = "carbon";
    else
        fprintf('Please enter count, carbon, 1, or 2.\n');
    end
end
end

function outputDir = choose_output_dir()
answer = strtrim(string(input('Output folder relative to current working directory [data]: ', 's')));
if strlength(answer) == 0
    answer = "data";
end
outputDir = answer;
end

function dataType = infer_data_type(inputMatFile, dataTableName)
tableText = lower(dataTableName);
fileText = lower(inputMatFile);

if contains(tableText, "count")
    dataType = "count";
elseif contains(tableText, "carbon") || contains(tableText, "classc")
    dataType = "carbon";
elseif contains(fileText, "count")
    dataType = "count";
elseif contains(fileText, "carbon")
    dataType = "carbon";
else
    dataType = "";
end
end
