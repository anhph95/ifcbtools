function community = make_community_array(data_wide, taxa_cols, time_step_col_name, site_id_col_name)
%MAKE_COMMUNITY_ARRAY Build a labeled metacommunity array from a wide table.
%
% community.values(t, i, j) = biomass of taxon j at time t and site i.
% Duplicate time x site rows are summed because multiple samples represent
% additive biomass observations for the same community cell.

if nargin < 3 || isempty(time_step_col_name)
    time_step_col_name = "timestep";
end
if nargin < 4 || isempty(site_id_col_name)
    site_id_col_name = "site";
end

if ~istable(data_wide)
    error("data_wide must be a MATLAB table.");
end

taxa_cols = string(taxa_cols);
time_step_col_name = string(time_step_col_name);
site_id_col_name = string(site_id_col_name);

time_values = data_wide.(time_step_col_name);
site_values = data_wide.(site_id_col_name);
time_ids = unique(string(time_values), "sorted");
site_ids = unique(string(site_values), "sorted");

X = zeros(numel(time_ids), numel(site_ids), numel(taxa_cols));

% Convert taxon columns to numeric biomass and treat missing values as zero.
% This preserves the ecological meaning: no measured biomass contributes 0
% to X(t,i,j), while duplicate rows still add to the same cell.
biomass = zeros(height(data_wide), numel(taxa_cols));
for j = 1:numel(taxa_cols)
    col = data_wide.(taxa_cols(j));
    if isnumeric(col)
        values = double(col);
    else
        values = str2double(string(col));
    end
    values(~isfinite(values)) = 0;
    biomass(:, j) = values;
end

for row_idx = 1:height(data_wide)
    t = find(time_ids == string(time_values(row_idx)), 1);
    i = find(site_ids == string(site_values(row_idx)), 1);
    X(t, i, :) = squeeze(X(t, i, :)) + biomass(row_idx, :)';
end

community = struct();
community.values = X;
community.timestep = time_ids;
community.site = site_ids;
community.taxon = taxa_cols(:)';
end
