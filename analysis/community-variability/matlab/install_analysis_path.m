function install_analysis_path(savePath)
%INSTALL_ANALYSIS_PATH Add MATLAB analysis scripts and metric functions.
%
% Run this from anywhere inside a checkout containing:
%   analysis/community-variability/matlab
%   community.variability/matlab/community_variability

if nargin < 1
    savePath = false;
end

analysisMatlabDir = fileparts(mfilename("fullpath"));
repoDir = fileparts(fileparts(fileparts(analysisMatlabDir)));
scriptDir = fullfile(analysisMatlabDir, "scripts");
metricDir = fullfile(repoDir, "community.variability", "matlab", "community_variability");

addpath(scriptDir);
if exist(metricDir, "dir")
    addpath(metricDir);
else
    warning("Community variability MATLAB folder not found: %s", metricDir);
end

fprintf("Added MATLAB analysis scripts:\n  %s\n", scriptDir);
fprintf("Added MATLAB metric functions:\n  %s\n", metricDir);

if savePath
    status = savepath;
    if status == 0
        fprintf("Saved MATLAB path.\n");
    else
        warning("MATLAB could not save the path. Add these folders manually in startup.m if needed.");
    end
end
end
