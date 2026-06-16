function out = cv_phi(X)
%CV_PHI Spatial aggregate synchrony.
%
% CV_phi = CV_gamma^2 / CV_alpha^2.
% Values near 1 mean local aggregate biomass fluctuations are synchronous;
% smaller values mean local sites compensate for one another through time.

if isstruct(X)
    values = double(X.values);
else
    values = double(X);
end
if ndims(values) ~= 3
    error("Community data must have shape X(time, site, taxon).");
end

regional_biomass = sum(sum(values, 3, "omitnan"), 2, "omitnan");
regional_biomass = squeeze(regional_biomass);
cv_g = (std(regional_biomass, 0, 1, "omitnan") / mean(regional_biomass, "omitnan"))^2;

site_biomass = sum(values, 3, "omitnan");
site_sd = std(site_biomass, 0, 1, "omitnan");
cv_a = (sum(site_sd(:), "omitnan") / mean(regional_biomass, "omitnan"))^2;

out = cv_g / cv_a;
end
