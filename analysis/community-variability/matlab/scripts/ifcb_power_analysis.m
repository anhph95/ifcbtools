%% IFCB seasonal bootstrap power analysis in MATLAB
%
% Bootstrap timesteps t in X(time, site, taxon), preserving the site x taxon
% community matrix within each sampled year.

workflowName = "ifcb_power_analysis_MATLAB";
run("ifcb_common.m")
fprintf("Starting MATLAB power-analysis workflow\n");

try

inputFile = fullfile(dataDir, "ifcb_carbon_clean_station_bottle_nutrient_fill.csv");
nBoot = 1000;
bootstrapSeed = 123;
logRunConfiguration(struct( ...
    "working_directory", pwd, ...
    "data_dir", dataDir, ...
    "results_dir", resultsDir, ...
    "input_file", inputFile, ...
    "seasons", seasons, ...
    "n_boot", nBoot, ...
    "bootstrap_seed", bootstrapSeed, ...
    "station_list", stationList, ...
    "main_cruise", mainCruise));

[df, taxaCols] = load_ifcb_carbon_power_local(inputFile, dataDir);
for seasonFilter = seasons
    fprintf("Processing season: %s\n", seasonFilter);
    ds = select_season_metacommunity_power_local(df, seasonFilter, stationList, mainCruise);
    fprintf("Complete years retained: %d\n", numel(unique(ds.year)));
    communityWide = ds(:, ["nearest_station", "year", taxaCols]);
    communityWide.Properties.VariableNames(1:2) = ["site", "timestep"];
    communityArray = make_community_array(communityWide, taxaCols);
    bootRes = bootstrap_by_dimension(communityArray, "timestep", nBoot, bootstrapSeed, true);
    writetable(bootRes.boot, fullfile(resultsDir, "boot_" + seasonFilter + ".csv"));
    writetable(bootRes.summary, fullfile(resultsDir, "summary_" + seasonFilter + ".csv"));
    writetable(spatial_bd_by_time(communityArray), fullfile(resultsDir, "spatial_variance_" + seasonFilter + ".csv"));
end
fprintf("MATLAB power-analysis workflow completed\n");
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

function [df, taxaCols] = load_ifcb_carbon_power_local(inputFile, dataDir)
carbonPath = inputFile;
taxonomyPath = fullfile(dataDir, "ifcb_taxonomy.csv");
df = readtable(carbonPath, "VariableNamingRule", "preserve");
taxonomy = readtable(taxonomyPath, "VariableNamingRule", "preserve");
taxonNameIdx = find(ismember(taxonomy.Properties.VariableNames, ["Label", "Annotations"]), 1);
taxonNames = string(taxonomy.(taxonomy.Properties.VariableNames{taxonNameIdx}));
taxaCols = taxonNames(ismember(taxonNames, string(df.Properties.VariableNames)))';
for j = 1:numel(taxaCols); df.(taxaCols(j)) = str2double(string(df.(taxaCols(j)))); end
df.sample_time = datetime(df.sample_time, "InputFormat", "yyyy-MM-dd HH:mm:ssXXX", "TimeZone", "UTC");
df = df(df.year ~= 2026, :);
[~, ~, groupIdx] = unique(string(df.cruise) + "|" + string(df.cast));
df.n_obs = accumarray(groupIdx, 1);
standIn = df(df.year == 2020 & string(df.season) == "JAS" & string(df.nearest_station) == "L9" & string(df.sample_type) == "underway", :);
if height(standIn) > 0
    [~, idx] = max(standIn.latitude); standIn = standIn(idx, :);
    standIn.nearest_station = "L8"; standIn.cast = "L8"; standIn.sample_type = "cast_from_udw";
    df = [df; standIn];
end
end

function ds = select_season_metacommunity_power_local(df, seasonFilter, stationList, mainCruise)
keep = ismember(string(df.nearest_station), stationList) & ismember(string(df.sample_type), ["cast", "cast_from_udw"]) & string(df.season) == seasonFilter;
ds = df(keep, :);
ds.station_order = double(categorical(string(ds.nearest_station), stationList, "Ordinal", true));
ds.main_cruise_priority = ismember(string(ds.cruise), mainCruise);
ds = sortrows(ds, ["year", "nearest_station", "cast", "depth"], ["ascend", "ascend", "ascend", "ascend"]);
[~, idx] = unique(string(ds.year) + "|" + string(ds.nearest_station) + "|" + string(ds.cast), "stable"); ds = ds(idx, :);
ds = sortrows(ds, ["year", "station_order", "n_obs", "main_cruise_priority"], ["ascend", "ascend", "descend", "descend"]);
[~, idx] = unique(string(ds.year) + "|" + string(ds.nearest_station), "stable"); ds = ds(idx, :);
years = unique(ds.year);
completeYears = years(arrayfun(@(y) numel(unique(string(ds.nearest_station(ds.year == y)))) == numel(stationList), years));
ds = ds(ismember(ds.year, completeYears), :);
ds.station_order = []; ds.main_cruise_priority = [];
end
