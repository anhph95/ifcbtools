%% IFCB seasonal leave-one-out sensitivity analysis in MATLAB
%
% Leave-one-year removes X(t, :, :). Leave-one-taxon removes X(:, :, j).
% Metrics are recalculated after each removal and written for comparison.

run("matlab/scripts/ifcb_common.m")

dataVersion = "fill";
[df, taxaCols] = load_ifcb_carbon_sensitivity_local(dataDir, dataVersion);
for seasonFilter = seasons
    fprintf("Processing season: %s\n", seasonFilter);
    ds = select_season_metacommunity_sensitivity_local(df, seasonFilter, stationList, mainCruise);
    communityWide = ds(:, ["nearest_station", "year", taxaCols]);
    communityWide.Properties.VariableNames(1:2) = ["site", "timestep"];
    communityArray = make_community_array(communityWide, taxaCols);
    yearOut = leave_one_out(communityArray, "timestep");
    yearOut.Properties.VariableNames(strcmp(yearOut.Properties.VariableNames, "timestep_removed")) = "year_removed";
    yearOut.season = repmat(seasonFilter, height(yearOut), 1);
    yearOut = movevars(yearOut, "season", "Before", 1);
    writetable(yearOut, fullfile(resultsDir, "leave_one_year_out_" + seasonFilter + ".csv"));
    taxonOut = leave_one_out(communityArray, "taxon");
    taxonOut.season = repmat(seasonFilter, height(taxonOut), 1);
    taxonOut = movevars(taxonOut, "season", "Before", 1);
    writetable(taxonOut, fullfile(resultsDir, "leave_one_taxon_out_" + seasonFilter + ".csv"));
end

function [df, taxaCols] = load_ifcb_carbon_sensitivity_local(dataDir, dataVersion)
carbonPath = fullfile(dataDir, "ifcb_carbon_" + dataVersion + ".csv");
taxonomyPath = fullfile(dataDir, "ifcb_taxonomy_" + dataVersion + ".csv");
if dataVersion == "clean"; taxonomyPath = fullfile(dataDir, "ifcb_taxonomy.csv"); end
df = readtable(carbonPath, "VariableNamingRule", "preserve");
taxonomy = readtable(taxonomyPath, "VariableNamingRule", "preserve");
taxonNameIdx = find(ismember(taxonomy.Properties.VariableNames, ["Label", "Annotations"]), 1);
taxonNames = string(taxonomy.(taxonomy.Properties.VariableNames{taxonNameIdx}));
taxaCols = taxonNames(ismember(taxonNames, string(df.Properties.VariableNames)))';
for j = 1:numel(taxaCols); df.(taxaCols(j)) = str2double(string(df.(taxaCols(j)))); end
df = df(df.year ~= 2026, :);
[~, ~, groupIdx] = unique(string(df.cruise) + "|" + string(df.cast)); df.n_obs = accumarray(groupIdx, 1);
standIn = df(df.year == 2020 & string(df.season) == "JAS" & string(df.nearest_station) == "L9" & string(df.sample_type) == "underway", :);
if height(standIn) > 0
    [~, idx] = max(standIn.latitude); standIn = standIn(idx, :);
    standIn.nearest_station = "L8"; standIn.cast = "L8"; standIn.sample_type = "cast_from_udw"; df = [df; standIn];
end
end

function ds = select_season_metacommunity_sensitivity_local(df, seasonFilter, stationList, mainCruise)
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
