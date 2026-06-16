function out = bd_spatial_weighted(X)
%BD_SPATIAL_WEIGHTED Biomass-weighted spatial beta diversity through time.
%
% For each timestep, BD_t^h = sum_j Var_i(z(t,i,j)), where
% z(t,i,j) = sqrt(X(t,i,j) / X(t,i,.)).
% BD_beta^h = sum_t W_t BD_t^h, with W_t proportional to total biomass.

if isstruct(X)
    values = double(X.values);
else
    values = double(X);
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
out = sum(bd(:) .* weights(:), "omitnan");
end
