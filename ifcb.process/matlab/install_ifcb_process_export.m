function install_ifcb_process_export(savePath)
%INSTALL_IFCB_PROCESS_EXPORT Add the IFCB MATLAB export script folder.
%
% This makes export_ifcb_mat.m available on the MATLAB path.
% savePath = true saves the path for future MATLAB sessions.

if nargin < 1
    savePath = false;
end

exportDir = fileparts(mfilename("fullpath"));
addpath(exportDir);
fprintf("Added IFCB process MATLAB export folder:\n  %s\n", exportDir);

if savePath
    status = savepath;
    if status == 0
        fprintf("Saved MATLAB path.\n");
    else
        warning("MATLAB could not save the path. Add this folder manually in startup.m if needed.");
    end
end
end
