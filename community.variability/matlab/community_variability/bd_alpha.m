function out = bd_alpha(X)
%BD_ALPHA Biomass-weighted local compositional variability.
%
% Local relative biomass is p(t,i,j) = X(t,i,j) / X(t,i,.).
% Hellinger composition is z(t,i,j) = sqrt(p(t,i,j)).
% Site-level BD_i^h = sum_j Var_t(z(t,i,j)).
% BD_alpha^h = sum_i w_i BD_i^h, where
% w_i = mean_t(X(t,i,.)) / sum_i mean_t(X(t,i,.)).

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

taxon_temporal_var = squeeze(var(site_hellinger, 0, 1, "omitnan"));
site_bd = sum(taxon_temporal_var, 2, "omitnan");

mean_site_biomass = squeeze(mean(site_biomass, 1, "omitnan"));
site_weights = mean_site_biomass ./ sum(mean_site_biomass, "omitnan");
out = sum(site_bd(:) .* site_weights(:), "omitnan");
end
