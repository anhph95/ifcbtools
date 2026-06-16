function metrics = calc_metacommunity_metrics(X)
%CALC_METACOMMUNITY_METRICS Calculate core aggregate and composition metrics.
%
% Output follows the R/Python long-table style:
%   varname:  CV_gamma, CV_alpha, CV_phi, BD_gamma, BD_alpha, BD_phi, BD_beta
%   estimate: numeric metric value

if isstruct(X)
    values = double(X.values);
else
    values = double(X);
end
if ndims(values) ~= 3
    error("Community data must have shape X(time, site, taxon).");
end

% Aggregate biomass partition.
% CV_gamma^2 tracks regional temporal variability.
% CV_alpha^2 sums local temporal variability before synchrony is applied.
site_biomass = sum(values, 3, "omitnan");
regional_biomass = squeeze(sum(site_biomass, 2, "omitnan"));
mu = mean(regional_biomass, "omitnan");
cv_g = (std(regional_biomass, 0, 1, "omitnan") / mu)^2;
site_sd = std(site_biomass, 0, 1, "omitnan");
cv_a = (sum(site_sd(:), "omitnan") / mu)^2;
cv_p = cv_g / cv_a;

% Regional compositional variability, BD_gamma^h.
regional_taxon_biomass = squeeze(sum(values, 2, "omitnan"));
regional_total = sum(regional_taxon_biomass, 2, "omitnan");
regional_relative = regional_taxon_biomass ./ regional_total;
regional_relative(~isfinite(regional_relative)) = 0;
regional_hellinger = sqrt(regional_relative);
bd_g = sum(var(regional_hellinger, 0, 1, "omitnan"), "omitnan");

% Local compositional variability, BD_alpha^h.
site_relative = values ./ site_biomass;
site_relative(~isfinite(site_relative)) = 0;
site_hellinger = sqrt(site_relative);
taxon_temporal_var = squeeze(var(site_hellinger, 0, 1, "omitnan"));
site_bd = sum(taxon_temporal_var, 2, "omitnan");
mean_site_biomass = squeeze(mean(site_biomass, 1, "omitnan"));
site_weights = mean_site_biomass ./ sum(mean_site_biomass, "omitnan");
bd_a = sum(site_bd(:) .* site_weights(:), "omitnan");
bd_p = bd_g / bd_a;

% Spatial beta diversity, BD_beta^h.
taxon_spatial_var = squeeze(var(site_hellinger, 0, 2, "omitnan"));
bd_by_time = sum(taxon_spatial_var, 2, "omitnan");
time_weights = regional_biomass ./ sum(regional_biomass, "omitnan");
bd_b = sum(bd_by_time(:) .* time_weights(:), "omitnan");

varname = ["CV_gamma"; "CV_alpha"; "CV_phi"; "BD_gamma"; "BD_alpha"; "BD_phi"; "BD_beta"];
estimate = [cv_g; cv_a; cv_p; bd_g; bd_a; bd_p; bd_b];
metrics = table(varname, estimate);
end
