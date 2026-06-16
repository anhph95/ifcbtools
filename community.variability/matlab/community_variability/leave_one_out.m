function out = leave_one_out(X, margin)
%LEAVE_ONE_OUT Recalculate metrics after removing one slice at a time.
%
% margin can be "timestep", "site", "taxon", or integer 1, 2, 3.
% This estimates sensitivity: delta = metric_without_slice - baseline.

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

margin_name = names(margin_idx);
removed_col = margin_name + "_removed";
slice_labels = labels.(margin_name);

baseline = calc_metacommunity_metrics(struct("values", values, ...
    "timestep", labels.timestep, "site", labels.site, "taxon", labels.taxon));
baseline.(removed_col) = repmat("Baseline", height(baseline), 1);
baseline = movevars(baseline, removed_col, "Before", 1);
out = baseline;

for label_idx = 1:numel(slice_labels)
    keep_idx = setdiff(1:numel(slice_labels), label_idx, "stable");
    if margin_idx == 1
        omitted = struct("values", values(keep_idx, :, :), "timestep", labels.timestep(keep_idx), ...
            "site", labels.site, "taxon", labels.taxon);
    elseif margin_idx == 2
        omitted = struct("values", values(:, keep_idx, :), "timestep", labels.timestep, ...
            "site", labels.site(keep_idx), "taxon", labels.taxon);
    else
        omitted = struct("values", values(:, :, keep_idx), "timestep", labels.timestep, ...
            "site", labels.site, "taxon", labels.taxon(keep_idx));
    end

    metrics = calc_metacommunity_metrics(omitted);
    metrics.(removed_col) = repmat(string(slice_labels(label_idx)), height(metrics), 1);
    metrics = movevars(metrics, removed_col, "Before", 1);
    out = [out; metrics]; %#ok<AGROW>
end
end
