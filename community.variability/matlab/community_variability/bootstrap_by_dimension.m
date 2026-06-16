function result = bootstrap_by_dimension(X, margin, n_boot, seed, baseline_in_boot)
%BOOTSTRAP_BY_DIMENSION Bootstrap by resampling one array dimension.
%
% margin can be "timestep", "site", "taxon", or integer 1, 2, 3.
% Each bootstrap replicate samples slices with replacement and recalculates
% the core metacommunity metrics.

if nargin < 3 || isempty(n_boot)
    n_boot = 1000;
end
if nargin < 4 || isempty(seed)
    seed = 123;
end
if nargin < 5 || isempty(baseline_in_boot)
    baseline_in_boot = true;
end

if isstruct(X)
    values = double(X.values);
    labels = X;
else
    values = double(X);
    labels = struct();
    labels.timestep = string(1:size(values, 1));
    labels.site = string(1:size(values, 2));
    labels.taxon = string(1:size(values, 3));
end
if ndims(values) ~= 3
    error("Community data must have shape X(time, site, taxon).");
end

names = ["timestep", "site", "taxon"];
if isstring(margin) || ischar(margin)
    margin_idx = find(names == string(margin), 1);
    if isempty(margin_idx)
        error("margin must be timestep, site, taxon, 1, 2, or 3.");
    end
else
    margin_idx = double(margin);
    if ~ismember(margin_idx, [1, 2, 3])
        error("Integer margin must be 1, 2, or 3.");
    end
end

rng(seed);
n_margin = size(values, margin_idx);
baseline = calc_metacommunity_metrics(struct("values", values, ...
    "timestep", labels.timestep, "site", labels.site, "taxon", labels.taxon));

boot_replicates = table();
for boot_id = 1:n_boot
    sampled_idx = randi(n_margin, [1, n_margin]);
    if margin_idx == 1
        boot_values = values(sampled_idx, :, :);
        boot_labels = struct("values", boot_values, "timestep", labels.timestep(sampled_idx), ...
            "site", labels.site, "taxon", labels.taxon);
    elseif margin_idx == 2
        boot_values = values(:, sampled_idx, :);
        boot_labels = struct("values", boot_values, "timestep", labels.timestep, ...
            "site", labels.site(sampled_idx), "taxon", labels.taxon);
    else
        boot_values = values(:, :, sampled_idx);
        boot_labels = struct("values", boot_values, "timestep", labels.timestep, ...
            "site", labels.site, "taxon", labels.taxon(sampled_idx));
    end

    metric_long = calc_metacommunity_metrics(boot_labels);
    metric_wide = array2table(metric_long.estimate', "VariableNames", cellstr(metric_long.varname'));
    metric_wide.boot_id = boot_id;
    metric_wide.sample_type = "Bootstrap";
    metric_wide = movevars(metric_wide, ["boot_id", "sample_type"], "Before", 1);
    boot_replicates = [boot_replicates; metric_wide]; %#ok<AGROW>
end

baseline_wide = array2table(baseline.estimate', "VariableNames", cellstr(baseline.varname'));
baseline_wide.boot_id = 0;
baseline_wide.sample_type = "Baseline";
baseline_wide = movevars(baseline_wide, ["boot_id", "sample_type"], "Before", 1);
if baseline_in_boot
    boot = [baseline_wide; boot_replicates];
else
    boot = boot_replicates;
end

summary = baseline;
summary.lwr = nan(height(summary), 1);
summary.upr = nan(height(summary), 1);
for row_idx = 1:height(summary)
    metric_name = summary.varname(row_idx);
    values_for_metric = sort(boot_replicates.(metric_name));
    values_for_metric = values_for_metric(isfinite(values_for_metric));
    if isempty(values_for_metric)
        summary.lwr(row_idx) = NaN;
        summary.upr(row_idx) = NaN;
    else
        % Dependency-free empirical quantiles with linear interpolation.
        % q = 0.025 and 0.975 define the bootstrap confidence interval.
        n = numel(values_for_metric);
        lo_pos = 1 + (n - 1) * 0.025;
        hi_pos = 1 + (n - 1) * 0.975;
        lo_floor = floor(lo_pos);
        hi_floor = floor(hi_pos);
        lo_ceil = ceil(lo_pos);
        hi_ceil = ceil(hi_pos);
        summary.lwr(row_idx) = values_for_metric(lo_floor) + ...
            (lo_pos - lo_floor) * (values_for_metric(lo_ceil) - values_for_metric(lo_floor));
        summary.upr(row_idx) = values_for_metric(hi_floor) + ...
            (hi_pos - hi_floor) * (values_for_metric(hi_ceil) - values_for_metric(hi_floor));
    end
end

result = struct("baseline", baseline, "boot", boot, ...
    "boot_replicates", boot_replicates, "summary", summary);
end
