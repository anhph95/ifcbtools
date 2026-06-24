%% Shared IFCB community-variability workflow settings for MATLAB scripts
%
% Run analysis scripts from analysis/community-variability.
% The shared ecological array convention is:
%   X(time, site, taxon)

analysisDir = pwd;
resultsDir = fullfile(analysisDir, 'results');
if ~exist(resultsDir, 'dir')
    mkdir(resultsDir);
end

dataDir = fullfile(analysisDir, 'data', 'NESLTER_transect');

seasons = ["JFM", "AMJ", "JAS", "OND"];
stationList = ["L1", "L2", "L3", "L4", "L5", "L6", "L7", "L8", "L9", "L10", "L11"];
mainCruise = [
    "EN608", "AR28B", "EN617", "AR32", "EN627", "AR34B", "EN644", ...
    "AR39B", "EN649", "AR44", "EN655", "EN657", "EN661", "AR52B", ...
    "EN668", "AR61B", "AT46", "AR66B", "EN687", "AR70B", "EN695", ...
    "HRS2303", "EN706", "AR77", "EN712", "EN715", "EN720", "AE2426", ...
    "EN727", "AR88", "AR92", "AR95", "AR99"
];

%% Shared dependency-free workflow logging
%
% The calling script defines workflowName before running this file. MATLAB's
% diary records command-window output and warnings under ./logs from the
% invocation directory.
if ~exist("workflowName", "var")
    workflowName = "ifcb_matlab";
end
logDir = fullfile(pwd, "logs");
if ~exist(logDir, "dir")
    mkdir(logDir);
end
logTimestamp = string(datetime("now", "Format", "yyyyMMdd_HHmmss"));
outLogPath = fullfile(logDir, workflowName + "_" + logTimestamp + ".out.log");
errLogPath = fullfile(logDir, workflowName + "_" + logTimestamp + ".err.log");
fclose(fopen(errLogPath, "a"));
diary(outLogPath);
fprintf("%s | INFO | %s | Logging to: %s\n", ...
    string(datetime("now", "Format", "yyyy-MM-dd HH:mm:ss")), workflowName, outLogPath);
fprintf("%s | INFO | %s | Errors to: %s\n", ...
    string(datetime("now", "Format", "yyyy-MM-dd HH:mm:ss")), workflowName, errLogPath);

% logRunConfiguration writes user-controlled workflow settings in a uniform
% block. Callers pass ordinary names and values after defining their options.
logRunConfiguration = @log_run_configuration_local;

function log_run_configuration_local(settings)
secretTerms = ["password", "token", "secret", "api_key", "apikey"];
fprintf("Run configuration:\n");
names = fieldnames(settings);
for idx = 1:numel(names)
    name = names{idx};
    value = settings.(name);
    if any(contains(lower(string(name)), secretTerms))
        text = "<redacted>";
    elseif iscell(value)
        text = strjoin(string(value), ", ");
    elseif numel(value) > 1
        text = strjoin(string(value), ", ");
    else
        text = string(value);
    end
    fprintf("  %s: %s\n", name, text);
end
end
