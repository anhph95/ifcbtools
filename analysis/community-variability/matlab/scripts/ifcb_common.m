%% Shared IFCB community-variability workflow settings for MATLAB scripts
%
% Run analysis scripts from analysis/community-variability.
% The shared ecological array convention is:
%   X(time, site, taxon)

analysisDir = pwd;
repoDir = fileparts(fileparts(analysisDir));
resultsDir = fullfile(analysisDir, 'results');
if ~exist(resultsDir, 'dir')
    mkdir(resultsDir);
end

dataDir = fullfile(repoDir, 'data', 'NESLTER_transect');
communityVariabilityMatlabDir = fullfile(repoDir, 'community.variability', 'matlab', 'community_variability');
addpath(communityVariabilityMatlabDir);

seasons = ["JFM", "AMJ", "JAS", "OND"];
stationList = ["L1", "L2", "L3", "L4", "L5", "L6", "L7", "L8", "L9", "L10", "L11"];
mainCruise = [
    "EN608", "AR28B", "EN617", "AR32", "EN627", "AR34B", "EN644", ...
    "AR39B", "EN649", "AR44", "EN655", "EN657", "EN661", "AR52B", ...
    "EN668", "AR61B", "AT46", "AR66B", "EN687", "AR70B", "EN695", ...
    "HRS2303", "EN706", "AR77", "EN712", "EN715", "EN720", "AE2426", ...
    "EN727", "AR88", "AR92", "AR95", "AR99"
];
