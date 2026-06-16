function out = add_baseline_delta(sensitivity_results, removed_col, baseline_label, group_cols)
%ADD_BASELINE_DELTA Add metric change relative to the baseline row.
%
% delta = estimate_without_component - estimate_baseline.
% abs_delta is useful for ranking influential timesteps, sites, or taxa.

if nargin < 3 || isempty(baseline_label)
    baseline_label = "Baseline";
end
if nargin < 4
    group_cols = strings(1, 0);
end

removed_col = string(removed_col);
baseline_label = string(baseline_label);
group_cols = string(group_cols);

out = sensitivity_results;
out.baseline = nan(height(out), 1);
out.delta = nan(height(out), 1);
out.abs_delta = nan(height(out), 1);

for row_idx = 1:height(out)
    same_metric = string(out.varname) == string(out.varname(row_idx));
    same_group = true(height(out), 1);
    for group_idx = 1:numel(group_cols)
        group_values = out.(group_cols(group_idx));
        if isnumeric(group_values) || islogical(group_values)
            same_group = same_group & group_values == group_values(row_idx);
        else
            same_group = same_group & string(group_values) == string(group_values(row_idx));
        end
    end
    baseline_row = same_metric & same_group & string(out.(removed_col)) == baseline_label;
    if any(baseline_row)
        out.baseline(row_idx) = out.estimate(find(baseline_row, 1));
        out.delta(row_idx) = out.estimate(row_idx) - out.baseline(row_idx);
        out.abs_delta(row_idx) = abs(out.delta(row_idx));
    end
end
end
