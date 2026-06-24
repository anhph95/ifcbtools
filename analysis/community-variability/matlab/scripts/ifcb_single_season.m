%% IFCB single-season community-variability workflow in MATLAB
%
% This script prepares one balanced seasonal metacommunity:
%   X(time, site, taxon)
% and writes observed alpha-gamma-phi metric estimates plus composition data.

workflowName = "ifcb_single_season_MATLAB";
run("ifcb_common.m")
fprintf("Starting MATLAB single-season workflow\n");

try

inputFile = fullfile(dataDir, "ifcb_carbon_clean_station_bottle_nutrient_fill.csv");
seasonFilter = "JAS";
topTaxaPerStation = 3;
logRunConfiguration(struct( ...
    "working_directory", pwd, ...
    "data_dir", dataDir, ...
    "results_dir", resultsDir, ...
    "input_file", inputFile, ...
    "season", seasonFilter, ...
    "top_taxa_per_station", topTaxaPerStation, ...
    "station_list", stationList, ...
    "main_cruise", mainCruise));

[df, taxaCols] = load_ifcb_carbon_local(inputFile, dataDir);
ds = select_season_metacommunity_local(df, seasonFilter, stationList, mainCruise);
fprintf("Complete years retained for %s: %d\n", seasonFilter, numel(unique(ds.year)));

communityWide = ds(:, ["nearest_station", "year", taxaCols]);
communityWide.Properties.VariableNames(1:2) = ["site", "timestep"];
communityArray = make_community_array(communityWide, taxaCols);
metrics = calc_metacommunity_metrics(communityArray);
disp(metrics)
writetable(rows2vars(metrics, "VariableNamesSource", "varname", "DataVariable", "estimate"), ...
    fullfile(resultsDir, "estimate_" + seasonFilter + ".csv"));

dsWithMeta = add_metacommunity_rows_local(ds, taxaCols);
topSpecies = dominant_taxa_by_station_local(ds, taxaCols, topTaxaPerStation);
composition = prepare_composition_long_local(dsWithMeta, taxaCols, topSpecies);
writetable(composition, fullfile(resultsDir, "composition_" + seasonFilter + ".csv"));
fprintf("MATLAB single-season workflow completed\n");
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

function [df, taxaCols] = load_ifcb_carbon_local(inputFile, dataDir)
carbonPath = inputFile;
taxonomyPath = fullfile(dataDir, "ifcb_taxonomy.csv");
df = readtable(carbonPath, "VariableNamingRule", "preserve");
taxonomy = readtable(taxonomyPath, "VariableNamingRule", "preserve");
if any(strcmp("Label", taxonomy.Properties.VariableNames))
    taxonNames = string(taxonomy.Label);
else
    taxonNames = string(taxonomy.Annotations);
end
taxaCols = taxonNames(ismember(taxonNames, string(df.Properties.VariableNames)))';
for j = 1:numel(taxaCols)
    df.(taxaCols(j)) = str2double(string(df.(taxaCols(j))));
end
df.sample_time = datetime(df.sample_time, "InputFormat", "yyyy-MM-dd HH:mm:ssXXX", "TimeZone", "UTC");
df = df(df.year ~= 2026, :);
[~, ~, groupIdx] = unique(string(df.cruise) + "|" + string(df.cast));
df.n_obs = accumarray(groupIdx, 1);
standIn = df(df.year == 2020 & string(df.season) == "JAS" & string(df.nearest_station) == "L9" & string(df.sample_type) == "underway", :);
if height(standIn) > 0
    [~, idx] = max(standIn.latitude);
    standIn = standIn(idx, :);
    standIn.nearest_station = repmat("L8", height(standIn), 1);
    standIn.cast = repmat("L8", height(standIn), 1);
    standIn.sample_type = repmat("cast_from_udw", height(standIn), 1);
    df = [df; standIn];
end
end

function ds = select_season_metacommunity_local(df, seasonFilter, stationList, mainCruise)
keep = ismember(string(df.nearest_station), stationList) & ...
    ismember(string(df.sample_type), ["cast", "cast_from_udw"]) & ...
    string(df.season) == seasonFilter;
ds = df(keep, :);
ds.station_order = double(categorical(string(ds.nearest_station), stationList, "Ordinal", true));
ds.main_cruise_priority = ismember(string(ds.cruise), mainCruise);
ds = sortrows(ds, ["year", "nearest_station", "cast", "depth"], ["ascend", "ascend", "ascend", "ascend"]);
[~, idx] = unique(string(ds.year) + "|" + string(ds.nearest_station) + "|" + string(ds.cast), "stable");
ds = ds(idx, :);
ds = sortrows(ds, ["year", "station_order", "n_obs", "main_cruise_priority"], ["ascend", "ascend", "descend", "descend"]);
[~, idx] = unique(string(ds.year) + "|" + string(ds.nearest_station), "stable");
ds = ds(idx, :);
years = unique(ds.year);
completeYears = years(arrayfun(@(y) numel(unique(string(ds.nearest_station(ds.year == y)))) == numel(stationList), years));
ds = ds(ismember(ds.year, completeYears), :);
ds.station_order = [];
ds.main_cruise_priority = [];
end

function out = add_metacommunity_rows_local(ds, taxaCols)
years = unique(ds.year);
totals = table();
for y = years'
    rows = ds(ds.year == y, :);
    one = rows(1, :);
    one.nearest_station = "Metacommunity";
    one.sample_time = min(rows.sample_time);
    for j = 1:numel(taxaCols)
        one.(taxaCols(j)) = sum(rows.(taxaCols(j)), "omitnan");
    end
    totals = [totals; one]; %#ok<AGROW>
end
out = [ds; totals];
end

function topSpecies = dominant_taxa_by_station_local(ds, taxaCols, nPerStation)
selected = strings(1, 0);
stations = unique(string(ds.nearest_station));
for i = 1:numel(stations)
    rows = ds(string(ds.nearest_station) == stations(i), :);
    totals = zeros(1, numel(taxaCols));
    for j = 1:numel(taxaCols)
        totals(j) = sum(rows.(taxaCols(j)), "omitnan");
    end
    [~, idx] = sort(totals, "descend");
    selected = unique([selected, taxaCols(idx(1:min(nPerStation, numel(idx))))], "stable");
end
pooled = zeros(1, numel(selected));
for j = 1:numel(selected)
    pooled(j) = sum(ds.(selected(j)), "omitnan");
end
[~, idx] = sort(pooled, "descend");
topSpecies = selected(idx);
end

function out = prepare_composition_long_local(ds, taxaCols, topSpecies)
rows = table();
for j = 1:numel(taxaCols)
    species = repmat(taxaCols(j), height(ds), 1);
    speciesGrp = species;
    speciesGrp(~ismember(speciesGrp, topSpecies)) = "Others";
    tmp = table(string(ds.nearest_station), ds.year, species, speciesGrp, ds.(taxaCols(j)), ...
        "VariableNames", ["site", "year", "species", "species_grp", "count"]);
    rows = [rows; tmp]; %#ok<AGROW>
end
out = groupsummary(rows, ["site", "year", "species_grp"], "sum", "count");
out.Properties.VariableNames(end) = "count";
end
