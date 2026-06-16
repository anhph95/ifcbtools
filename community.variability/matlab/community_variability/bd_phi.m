function out = bd_phi(X)
%BD_PHI Spatial compositional synchrony.
%
% BD_phi^h = BD_gamma^h / BD_alpha^h.
% It compares regional compositional turnover with biomass-weighted local
% compositional turnover.

if isstruct(X)
    values = double(X.values);
else
    values = double(X);
end
if ndims(values) ~= 3
    error("Community data must have shape X(time, site, taxon).");
end

regional_taxon_biomass = squeeze(sum(values, 2, "omitnan"));
regional_total = sum(regional_taxon_biomass, 2, "omitnan");
regional_relative = regional_taxon_biomass ./ regional_total;
regional_relative(~isfinite(regional_relative)) = 0;
regional_hellinger = sqrt(regional_relative);
bd_g = sum(var(regional_hellinger, 0, 1, "omitnan"), "omitnan");

site_biomass = sum(values, 3, "omitnan");
site_relative = values ./ site_biomass;
site_relative(~isfinite(site_relative)) = 0;
site_hellinger = sqrt(site_relative);
taxon_temporal_var = squeeze(var(site_hellinger, 0, 1, "omitnan"));
site_bd = sum(taxon_temporal_var, 2, "omitnan");
mean_site_biomass = squeeze(mean(site_biomass, 1, "omitnan"));
site_weights = mean_site_biomass ./ sum(mean_site_biomass, "omitnan");
bd_a = sum(site_bd(:) .* site_weights(:), "omitnan");

out = bd_g / bd_a;
end
