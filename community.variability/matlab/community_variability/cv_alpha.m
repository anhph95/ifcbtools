function out = cv_alpha(X)
%CV_ALPHA Biomass-weighted local aggregate variability.
%
% Site biomass is X_(ti.) = sum_j X(t,i,j).
% CV_alpha^2 = [sum_i sd_t(X_(ti.)) / mean_t(X_(t..))]^2.
% This measures the total local temporal fluctuation before accounting for
% whether sites fluctuate synchronously.

if isstruct(X)
    values = double(X.values);
else
    values = double(X);
end
if ndims(values) ~= 3
    error("Community data must have shape X(time, site, taxon).");
end

site_biomass = sum(values, 3, "omitnan");
regional_biomass = sum(site_biomass, 2, "omitnan");
regional_biomass = squeeze(regional_biomass);
mu = mean(regional_biomass, "omitnan");
site_sd = std(site_biomass, 0, 1, "omitnan");
out = (sum(site_sd(:), "omitnan") / mu)^2;
end
