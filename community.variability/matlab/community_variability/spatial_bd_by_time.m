function out = spatial_bd_by_time(X)
%SPATIAL_BD_BY_TIME Spatial compositional variability at each timestep.
%
% For each time t, calculate local Hellinger composition z(t,i,j), then
% BD_t^h = sum_j Var_i(z(t,i,j)).
% The returned weights are regional biomass weights through time.

if isstruct(X)
    values = double(X.values);
    timestep = string(X.timestep(:));
else
    values = double(X);
    timestep = string((1:size(values, 1))');
end
if ndims(values) ~= 3
    error("Community data must have shape X(time, site, taxon).");
end

site_biomass = sum(values, 3, "omitnan");
site_relative = values ./ site_biomass;
site_relative(~isfinite(site_relative)) = 0;
site_hellinger = sqrt(site_relative);

taxon_spatial_var = squeeze(var(site_hellinger, 0, 2, "omitnan"));
bd = sum(taxon_spatial_var, 2, "omitnan");

total_metacommunity_biomass = squeeze(sum(site_biomass, 2, "omitnan"));
weights = total_metacommunity_biomass ./ sum(total_metacommunity_biomass, "omitnan");
out = table(timestep, bd, total_metacommunity_biomass, weights, bd .* weights, ...
    "VariableNames", ["timestep", "BD", "total_metacommunity_biomass", "weights", "BD_x_wt"]);
end
