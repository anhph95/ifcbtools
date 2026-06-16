%% IFCB seasonal comparison summaries in MATLAB
%
% Reads saved bootstrap, sensitivity, and spatial compositional variability
% outputs, then writes cross-season summary CSVs used by plotting workflows.

run("matlab/scripts/ifcb_common.m")

communityMetricOrder = ["CV_gamma", "CV_alpha", "CV_phi", "BD_gamma", "BD_alpha", "BD_phi", "BD_beta"];
plotVarsSensitivity = ["BD_gamma", "CV_gamma", "BD_phi", "CV_phi"];
topNTaxa = 10;

summaryAll = table();
bootAll = table();
for seasonFilter = seasons
    s = readtable(fullfile(resultsDir, "summary_" + seasonFilter + ".csv"));
    s.season = repmat(seasonFilter, height(s), 1);
    summaryAll = [summaryAll; s]; %#ok<AGROW>
    b = readtable(fullfile(resultsDir, "boot_" + seasonFilter + ".csv"));
    b.season = repmat(seasonFilter, height(b), 1);
    bootAll = [bootAll; b]; %#ok<AGROW>
end
bootAll = bootAll(string(bootAll.sample_type) == "Bootstrap", :);

summaryWide = unstack(summaryAll(:, ["season", "varname", "estimate"]), "estimate", "varname");
writetable(summaryWide, fullfile(resultsDir, "seasonal_metric_summary_wide.csv"));

looYearAll = table();
looTaxonAll = table();
for seasonFilter = seasons
    looYearAll = [looYearAll; readtable(fullfile(resultsDir, "leave_one_year_out_" + seasonFilter + ".csv"))]; %#ok<AGROW>
    looTaxonAll = [looTaxonAll; readtable(fullfile(resultsDir, "leave_one_taxon_out_" + seasonFilter + ".csv"))]; %#ok<AGROW>
end
looYearAll = add_baseline_delta_local(looYearAll, "year_removed", "season");
looTaxonAll = add_baseline_delta_local(looTaxonAll, "taxon_removed", "season");
writetable(looYearAll, fullfile(resultsDir, "leave_one_year_all_seasons_with_delta.csv"));
writetable(looTaxonAll, fullfile(resultsDir, "leave_one_taxon_all_seasons_with_delta.csv"));
writetable(classify_taxon_effects_local(looTaxonAll, plotVarsSensitivity, topNTaxa), ...
    fullfile(resultsDir, "leave_one_taxon_top_effects_all_seasons.csv"));

ellipsePairs = ["BD_alpha", "CV_alpha"; "BD_phi", "CV_phi"; "BD_gamma", "CV_gamma"; "BD_beta", "BD_gamma"];
ellipses = table();
for k = 1:size(ellipsePairs, 1)
    ellipses = [ellipses; make_ellipse_local(bootAll, seasons, ellipsePairs(k, 1), ellipsePairs(k, 2))]; %#ok<AGROW>
end
writetable(ellipses, fullfile(resultsDir, "bootstrap_metric_ellipses.csv"));

spavarAll = table();
for idx = 1:numel(seasons)
    seasonFilter = seasons(idx);
    s = readtable(fullfile(resultsDir, "spatial_variance_" + seasonFilter + ".csv"));
    s.season = repmat(seasonFilter, height(s), 1);
    s.year = str2double(string(s.timestep));
    s.time_index = s.year + (idx - 1) / numel(seasons);
    spavarAll = [spavarAll; s]; %#ok<AGROW>
end
writetable(spavarAll, fullfile(resultsDir, "spatial_variance_all_seasons.csv"));

function out = add_baseline_delta_local(dat, removedCol, groupCol)
out = dat;
out.baseline = nan(height(out), 1);
out.delta = nan(height(out), 1);
out.abs_delta = nan(height(out), 1);
for i = 1:height(out)
    sameMetric = string(out.varname) == string(out.varname(i));
    sameGroup = string(out.(groupCol)) == string(out.(groupCol)(i));
    baselineRow = sameMetric & sameGroup & string(out.(removedCol)) == "Baseline";
    if any(baselineRow)
        baselineValue = out.estimate(find(baselineRow, 1));
        out.baseline(i) = baselineValue;
        out.delta(i) = out.estimate(i) - baselineValue;
        out.abs_delta(i) = abs(out.delta(i));
    end
end
end

function out = classify_taxon_effects_local(dat, plotVars, topN)
rows = table();
for season = unique(string(dat.season))'
    for varname = plotVars
        d = dat(string(dat.season) == season & string(dat.varname) == varname & string(dat.taxon_removed) ~= "Baseline", :);
        if height(d) == 0; continue; end
        q = quantile(d.delta, [0.05, 0.95]);
        d.effect_class = repmat("Not significant", height(d), 1);
        d.effect_class(d.delta < q(1)) = "Negative";
        d.effect_class(d.delta > q(2)) = "Positive";
        d = sortrows(d, "abs_delta", "descend");
        rows = [rows; d(1:min(topN, height(d)), :)]; %#ok<AGROW>
    end
end
out = rows;
end

function out = make_ellipse_local(dat, seasons, xvar, yvar)
theta = linspace(0, 2*pi, 200)';
r2 = 5.991464547107979;
out = table();
for season = seasons
    d = dat(string(dat.season) == season, :);
    xy = [d.(xvar), d.(yvar)];
    xy = xy(all(isfinite(xy), 2), :);
    if size(xy, 1) < 3; continue; end
    S = cov(xy);
    if any(~isfinite(S), "all") || det(S) <= 0; continue; end
    mu = mean(xy, 1);
    [V, D] = eig(S);
    A = V * diag(sqrt(max(diag(D), 0) * r2));
    pts = mu + [cos(theta), sin(theta)] * A';
    tmp = table(repmat(season, size(pts, 1), 1), pts(:, 1), pts(:, 2), repmat(xvar, size(pts, 1), 1), repmat(yvar, size(pts, 1), 1), ...
        "VariableNames", ["season", "x", "y", "xvar", "yvar"]);
    out = [out; tmp]; %#ok<AGROW>
end
end
