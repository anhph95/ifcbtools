function export_ifcb(dataset, options)
%EXPORT_IFCB Export IFCB summary MAT products to processing CSV files.
%
% export_ifcb()
% export_ifcb("NESLTER_broadscale")
% export_ifcb("NESLTER_broadscale", SummaryDir="path/to/summary")
% export_ifcb("NESLTER_broadscale", ...
%     CountMatFile="path/to/count.mat", ...
%     CarbonMatFile="path/to/carbon.mat", ...
%     OutputDir="path/to/output")
%
% Default MAT filenames:
%   count_group_class_withTS.mat
%   carbon_group_class_withTS.mat
%
% Outputs:
%   <OutputDir>/ifcb_class.csv
%   <OutputDir>/ifcb_metadata.csv
%   <OutputDir>/ifcb_count_raw.csv
%   <OutputDir>/ifcb_carbon_raw.csv
%
% Scientific data model:
%   rows    = IFCB samples or sample-derived records
%   columns = metadata fields or taxon/class biomass/count variables

arguments
    dataset (1, 1) string = "NESLTER_transect"
    options.SummaryDir (1, 1) string = ""
    options.CountMatFile (1, 1) string = ""
    options.CarbonMatFile (1, 1) string = ""
    options.OutputDir (1, 1) string = ""
end

%% --- RESOLVE INPUT AND OUTPUT PATHS ---

if strlength(options.SummaryDir) == 0
    summaryDir = fullfile("\\sosiknas1", "IFCB_products", dataset, "summary");
else
    summaryDir = options.SummaryDir;
end

if strlength(options.CountMatFile) == 0
    countMatFile = fullfile(summaryDir, "count_group_class_withTS.mat");
else
    countMatFile = options.CountMatFile;
end

if strlength(options.CarbonMatFile) == 0
    carbonMatFile = fullfile(summaryDir, "carbon_group_class_withTS.mat");
else
    carbonMatFile = options.CarbonMatFile;
end

if strlength(options.OutputDir) == 0
    outputDir = fullfile(pwd, "data", dataset);
else
    outputDir = options.OutputDir;
end

fileList = {countMatFile, carbonMatFile};

%% --- ENSURE OUTPUT DIRECTORY EXISTS ---

if ~exist(outputDir, 'dir')
    fprintf('Creating output directory: %s\n', outputDir);
    mkdir(outputDir);
end

%% --- CONFIGURE WORKFLOW LOGGING ---

logDir = fullfile(outputDir, 'logs');
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
fprintf("  dataset: %s\n", dataset);
fprintf("  summary_dir: %s\n", summaryDir);
fprintf("  output_dir: %s\n", outputDir);
fprintf("  input_files: %s\n", strjoin(string(fileList), ", "));
fprintf("  log_dir: %s\n", logDir);

savedMetadataAndClass = false;

%% --- PROCESS EACH MAT FILE ---

try
for f = 1:length(fileList)

    fprintf('\n[%d/%d] Loading file: %s\n', f, length(fileList), fileList{f});
    loadedData = load(fileList{f});

    %% Convert metadata timestamps if present.
    %
    % IFCB timestamps are stored as UTC strings. Removing the '+00:00'
    % suffix allows datetime() to parse the values while preserving UTC:
    %   sample_time -> t_sample, with timezone = UTC
    if isfield(loadedData, 'meta_data') && isfield(loadedData.meta_data, 'sample_time')
        rawTimes = loadedData.meta_data.sample_time;
        trimmedTimes = erase(rawTimes, '+00:00');
        loadedData.meta_data.sample_time = datetime(trimmedTimes, ...
            'InputFormat', 'yyyy-MM-dd HH:mm:ss', ...
            'TimeZone','UTC');
    end

    %% Extract class labels and metadata once.
    %
    % The *_label variables define the class groups used by the raw count
    % and carbon tables. Their exported table records:
    %   class = IFCB class/taxon name
    %   label = source label group in the MATLAB product
    if ~savedMetadataAndClass && isfield(loadedData, 'meta_data')
        fprintf('Extracting class labels...\n');

        loadedNames = fieldnames(loadedData);
        labelNames = loadedNames(endsWith(loadedNames, '_label'));
        classTable = table();

        for i = 1:numel(labelNames)
            thisName = labelNames{i};
            labelData = loadedData.(thisName);

            if ~iscell(labelData)
                fprintf('  Skipping non-cell label variable: %s\n', thisName);
                continue
            end

            groupName = erase(thisName, '_label');
            groupLabels = repmat({groupName}, size(labelData, 1), 1);

            tmp = table(labelData, groupLabels, ...
                'VariableNames', {'class', 'label'});

            classTable = [classTable; tmp]; %#ok<AGROW>
            fprintf('  Processed label group: %s (%d entries)\n', ...
                groupName, size(labelData, 1));
        end

        fprintf('Writing IFCB class table: %s\n', ...
            fullfile(outputDir, 'ifcb_class.csv'));
        writetable(classTable, fullfile(outputDir, 'ifcb_class.csv'));

        fprintf('Writing IFCB metadata table: %s\n', ...
            fullfile(outputDir, 'ifcb_metadata.csv'));
        writetable(loadedData.meta_data, fullfile(outputDir, 'ifcb_metadata.csv'));

        savedMetadataAndClass = true;
    end

    %% Save raw count table when present.
    %
    % Count data are cell abundance observations:
    %   N_{sample,class} = number of classified IFCB images/cells
    if isfield(loadedData, 'classcount_opt_adhoc_merge') && ...
            istable(loadedData.classcount_opt_adhoc_merge)
        outFile = fullfile(outputDir, 'ifcb_count_raw.csv');
        fprintf('Saving raw cell count table: %s\n', outFile);
        writetable(loadedData.classcount_opt_adhoc_merge, outFile);
    else
        fprintf('Count table not found in this file. Skipping.\n');
    end

    %% Save raw carbon table when present.
    %
    % Carbon data are biomass observations:
    %   C_{sample,class} = estimated carbon biomass for each IFCB class
    if isfield(loadedData, 'classC_opt_adhoc_merge') && ...
            istable(loadedData.classC_opt_adhoc_merge)
        outFile = fullfile(outputDir, 'ifcb_carbon_raw.csv');
        fprintf('Saving raw carbon table: %s\n', outFile);
        writetable(loadedData.classC_opt_adhoc_merge, outFile);
    else
        fprintf('Carbon table not found in this file. Skipping.\n');
    end
end

fprintf('\nDataset %s processed successfully.\n', dataset);
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
