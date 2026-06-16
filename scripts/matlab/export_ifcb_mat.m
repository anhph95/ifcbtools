% ========================================================================
% IFCB MATLAB PRODUCT EXPORT
%
% Export IFCB summary MAT products into CSV files used by the Python
% data processing step.
%
% Outputs:
%   - data/<dataset>/ifcb_class.csv
%   - data/<dataset>/ifcb_metadata.csv
%   - data/<dataset>/ifcb_count_raw.csv
%   - data/<dataset>/ifcb_carbon_raw.csv
%
% Scientific data model:
%   rows    = IFCB samples or sample-derived records
%   columns = metadata fields or taxon/class biomass/count variables
% ========================================================================

%% --- USER SETTINGS ---

dataset = 'NESLTER_broadscale';

summaryDir = fullfile( ...
    '\\sosiknas1', ...
    'IFCB_products', ...
    dataset, ...
    'summary');

fileList = {
    fullfile(summaryDir, 'count_group_class_withTS.mat')
    fullfile(summaryDir, 'carbon_group_class_withTS.mat')
};

%% --- REPOSITORY-RELATIVE OUTPUT SETTINGS ---

% This file lives in:
%   <repo>/scripts/matlab/export_ifcb_mat.m
%
% Move three levels upward to recover <repo>, then write local exports to:
%   <repo>/data/<dataset>/
repoDir = fileparts(fileparts(fileparts(mfilename('fullpath'))));
baseDir = fullfile(repoDir, 'data', dataset);

%% --- ENSURE OUTPUT DIRECTORY EXISTS ---

if ~exist(baseDir, 'dir')
    fprintf('Creating output directory: %s\n', baseDir);
    mkdir(baseDir);
end

savedMetadataAndClass = false;

%% --- PROCESS EACH MAT FILE ---

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
            fullfile(baseDir, 'ifcb_class.csv'));
        writetable(classTable, fullfile(baseDir, 'ifcb_class.csv'));

        fprintf('Writing IFCB metadata table: %s\n', ...
            fullfile(baseDir, 'ifcb_metadata.csv'));
        writetable(loadedData.meta_data, fullfile(baseDir, 'ifcb_metadata.csv'));

        savedMetadataAndClass = true;
    end

    %% Save raw count table when present.
    %
    % Count data are cell abundance observations:
    %   N_{sample,class} = number of classified IFCB images/cells
    if isfield(loadedData, 'classcount_opt_adhoc_merge') && ...
            istable(loadedData.classcount_opt_adhoc_merge)
        outFile = fullfile(baseDir, 'ifcb_count_raw.csv');
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
        outFile = fullfile(baseDir, 'ifcb_carbon_raw.csv');
        fprintf('Saving raw carbon table: %s\n', outFile);
        writetable(loadedData.classC_opt_adhoc_merge, outFile);
    else
        fprintf('Carbon table not found in this file. Skipping.\n');
    end
end

fprintf('\nAll IFCB files processed successfully.\n');
