function install_community_variability(savePath)
%INSTALL_COMMUNITY_VARIABILITY Add MATLAB community variability functions.
%
% The installed path contains functions for metacommunity arrays:
%   X(time, site, taxon)
%
% savePath = true saves the path for future MATLAB sessions.

if nargin < 1
    savePath = false;
end

thisFile = mfilename("fullpath");
matlabDir = fileparts(thisFile);
functionDir = fullfile(matlabDir, "community_variability");

if ~exist(functionDir, "dir")
    error("Cannot find community_variability function folder: %s", functionDir);
end

addpath(functionDir);
fprintf("Added community variability MATLAB functions:\n  %s\n", functionDir);

if savePath
    status = savepath;
    if status == 0
        fprintf("Saved MATLAB path.\n");
    else
        warning("MATLAB could not save the path. Add this folder manually in startup.m if needed.");
    end
end
end
