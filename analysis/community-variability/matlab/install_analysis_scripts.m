% Download editable IFCB MATLAB analysis starter scripts into the cwd.
%
% Set overwrite = true before running this script to replace existing copies.

if ~exist("overwrite", "var")
    overwrite = false;
end

baseUrl = "https://raw.githubusercontent.com/anhph95/ifcbtools/main/analysis/community-variability/matlab/scripts";
scriptFiles = [
    "ifcb_common.m"
    "ifcb_single_season.m"
    "ifcb_power_analysis.m"
    "ifcb_sensitivity_analysis.m"
    "ifcb_seasonal_comparison.m"
];

destDir = pwd;
fprintf("Installing MATLAB analysis scripts into:\n  %s\n", destDir);

for idx = 1:numel(scriptFiles)
    scriptFile = scriptFiles(idx);
    destFile = fullfile(destDir, scriptFile);
    if exist(destFile, "file") && ~overwrite
        fprintf("Skipping existing file: %s\n", scriptFile);
        continue
    end

    url = baseUrl + "/" + scriptFile;
    websave(destFile, url);
    fprintf("Downloaded: %s\n", scriptFile);
end

fprintf("Done. Open the copied scripts in this working folder and build on them.\n");
