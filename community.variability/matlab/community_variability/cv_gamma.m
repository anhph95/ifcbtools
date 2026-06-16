function out = cv_gamma(X)
%CV_GAMMA Regional aggregate variability.
%
% Let X(t,i,j) be biomass at time t, site i, taxon j.
% Regional biomass is X_(t..) = sum_i sum_j X(t,i,j).
% CV_gamma^2 = [sd_t(X_(t..)) / mean_t(X_(t..))]^2.

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
mu = mean(regional_biomass, "omitnan");
sigma = std(regional_biomass, 0, 1, "omitnan");
out = (sigma / mu)^2;
end
