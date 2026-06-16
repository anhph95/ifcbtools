function out = bd_gamma(X)
%BD_GAMMA Regional compositional variability using Hellinger composition.
%
% Regional taxon biomass is X_(t.j) = sum_i X(t,i,j).
% Regional relative biomass is p_(t.j) = X_(t.j) / X_(t..).
% Hellinger composition is z_(t.j) = sqrt(p_(t.j)).
% BD_gamma^h = sum_j Var_t(z_(t.j)).

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

taxon_temporal_var = var(regional_hellinger, 0, 1, "omitnan");
out = sum(taxon_temporal_var(:), "omitnan");
end
