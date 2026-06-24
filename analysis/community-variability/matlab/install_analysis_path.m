function install_analysis_path(savePath)
%INSTALL_ANALYSIS_PATH Add MATLAB analysis scripts and metric functions.
%
% Run this from anywhere after installing or adding the community variability
% metric functions separately.

if nargin < 1
    savePath = false;
end

analysisMatlabDir = fileparts(mfilename("fullpath"));
scriptDir = fullfile(analysisMatlabDir, "scripts");

addpath(scriptDir);

fprintf("Added MATLAB analysis scripts:\n  %s\n", scriptDir);
fprintf("Install or add community variability MATLAB metric functions separately.\n");

if savePath
    status = savepath;
    if status == 0
        fprintf("Saved MATLAB path.\n");
    else
        warning("MATLAB could not save the path. Add these folders manually in startup.m if needed.");
    end
end
end
