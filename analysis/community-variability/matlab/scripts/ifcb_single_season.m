%% IFCB single-season community-variability workflow in MATLAB
%
% This script prepares one balanced seasonal metacommunity:
%   X(time, site, taxon)
% and writes observed alpha-gamma-phi metric estimates plus composition data.

run("matlab/scripts/ifcb_common.m")

dataVersion = "fill";
seasonFilter = "JAS";
topTaxaPerStation = 3;

[df, taxaCols] = load_ifcb_carbon_local(dataDir, dataVersion);
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

function [df, taxaCols] = load_ifcb_carbon_local(dataDir, dataVersion)
carbonPath = fullfile(dataDir, "ifcb_carbon_" + dataVersion + ".csv");
if dataVersion == "fill"
    taxonomyPath = fullfile(dataDir, "ifcb_taxonomy_fill.csv");
else
    taxonomyPath = fullfile(dataDir, "ifcb_taxonomy.csv");
end
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
