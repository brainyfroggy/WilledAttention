%% Reviewer 4: EEG time-resolved decoding cluster inference
%
% Adds family-wise error control for time-resolved decoding curves. The
% original curve scripts save subject x time decoding accuracies, and older
% plotting code used pointwise tests against 50%. This script tests the same
% subject-level curves against chance using cluster-based permutation over
% time, reports effect sizes/CIs at every time point, and exports figures.

clear; clc; close all;

cfg = struct();
cfg.rootDir = fileparts(fileparts(fileparts(mfilename('fullpath'))));
cfg.outDir = fullfile(cfg.rootDir, 'Result Figures', 'Reviewer4_DecodingInference', 'EEG');
cfg.chance = 50;
cfg.alpha = 0.05;
cfg.clusterAlpha = 0.05;
cfg.nPerm = 2000;          % increase to 5000 or 10000 for final manuscript
cfg.randomSeed = 20260705;
cfg.preferFieldTrip = true;
cfg.smoothWindow = 1;      % set to 20 only to reproduce older smoothed plots

if strcmp(getenv('REVIEWER4_FAST'), '1')
    cfg.nPerm = 100;
end

if ~exist(cfg.outDir, 'dir')
    mkdir(cfg.outDir);
end

rng(cfg.randomSeed, 'twister');

fprintf('Reviewer 4 EEG cluster inference\n');
fprintf('Root: %s\n', cfg.rootDir);
fprintf('Output: %s\n', cfg.outDir);
fprintf('FieldTrip ft_timelockstatistics available: %d\n', exist('ft_timelockstatistics', 'file') == 2);
cfg.fieldTripClusterUsable = cfg.preferFieldTrip && checkFieldTripClusterUsable(cfg);
fprintf('FieldTrip cluster smoke test usable: %d\n', cfg.fieldTripClusterUsable);

analyses = defineAnalyses(cfg);
allResults = {};
timepointTables = {};
clusterTables = {};

for a = 1:numel(analyses)
    fprintf('\n%s\n', analyses(a).name);
    data = load(analyses(a).sourceFile);
    groups = makeGroups(data, analyses(a));
    for g = 1:numel(groups)
        result = runOneClusterTest(groups(g), analyses(a), cfg);
        allResults{end + 1, 1} = result; %#ok<SAGROW>
        timepointTables{end + 1} = result.timepointTable; %#ok<SAGROW>
        clusterTables{end + 1} = result.clusterTable; %#ok<SAGROW>
        plotDecodingCluster(result, cfg);
    end
end

timepointStats = vertcat(timepointTables{:});
clusterStats = vertcat(clusterTables{:});
responseText = makeResponseText(clusterStats, cfg);

writetable(timepointStats, fullfile(cfg.outDir, 'reviewer4_eeg_timepoint_stats.csv'));
writetable(clusterStats, fullfile(cfg.outDir, 'reviewer4_eeg_cluster_stats.csv'));

fid = fopen(fullfile(cfg.outDir, 'reviewer4_eeg_response_text.txt'), 'w');
fprintf(fid, '%s\n', responseText);
fclose(fid);

save(fullfile(cfg.outDir, 'reviewer4_eeg_cluster_inference.mat'), ...
    'cfg', 'analyses', 'allResults', 'timepointStats', 'clusterStats', 'responseText');

fprintf('\nDone. Results saved to:\n%s\n', cfg.outDir);
disp(clusterStats);

%% Local functions

function analyses = defineAnalyses(cfg)
analyses = struct([]);
analyses(1).name = 'PreCue_Choice';
analyses(1).sourceFile = fullfile(cfg.rootDir, 'EEG', 'precuedecoding.mat');
analyses(1).timeMs = -2500:4:0;
analyses(1).ufVar = 'AccuracyChoUF';
analyses(1).ucdVar = 'AccuracyChoUCD';
analyses(1).condition = 'Choice';
analyses(1).epoch = 'Pre-cue';

analyses(2).name = 'PreCue_Instructed';
analyses(2).sourceFile = fullfile(cfg.rootDir, 'EEG', 'precuedecoding.mat');
analyses(2).timeMs = -2500:4:0;
analyses(2).ufVar = 'AccuracyInsUF';
analyses(2).ucdVar = 'AccuracyInsUCD';
analyses(2).condition = 'Instructed';
analyses(2).epoch = 'Pre-cue';

analyses(3).name = 'PostCue_ERPRemoved_Choice';
analyses(3).sourceFile = fullfile(cfg.rootDir, 'EEG', 'erpremoved_postcuedecoding.mat');
analyses(3).timeMs = 0:4:1000;
analyses(3).ufVar = 'AccuracyChoUF';
analyses(3).ucdVar = 'AccuracyChoUCD';
analyses(3).condition = 'Choice';
analyses(3).epoch = 'Post-cue ERP-removed';

analyses(4).name = 'PostCue_ERPRemoved_Instructed';
analyses(4).sourceFile = fullfile(cfg.rootDir, 'EEG', 'erpremoved_postcuedecoding.mat');
analyses(4).timeMs = 0:4:1000;
analyses(4).ufVar = 'AccuracyInsUF';
analyses(4).ucdVar = 'AccuracyInsUCD';
analyses(4).condition = 'Instructed';
analyses(4).epoch = 'Post-cue ERP-removed';
end

function groups = makeGroups(data, analysis)
if ~isfield(data, analysis.ufVar) || ~isfield(data, analysis.ucdVar)
    error('Missing variables %s and/or %s in %s.', analysis.ufVar, analysis.ucdVar, analysis.sourceFile);
end
uf = double(data.(analysis.ufVar));
ucd = double(data.(analysis.ucdVar));
if size(uf, 2) ~= numel(analysis.timeMs) || size(ucd, 2) ~= numel(analysis.timeMs)
    error('Time axis length mismatch for %s.', analysis.name);
end

groups = struct([]);
groups(1).dataset = 'UF';
groups(1).accuracy = uf;
groups(1).subjectIds = arrayfun(@(x) sprintf('UF_S%02d', x), 1:size(uf, 1), 'UniformOutput', false)';
groups(2).dataset = 'UCD';
groups(2).accuracy = ucd;
groups(2).subjectIds = arrayfun(@(x) sprintf('UCD_S%02d', x), 1:size(ucd, 1), 'UniformOutput', false)';
groups(3).dataset = 'Combined';
groups(3).accuracy = [uf; ucd];
groups(3).subjectIds = [groups(1).subjectIds; groups(2).subjectIds];
end

function result = runOneClusterTest(group, analysis, cfg)
acc = group.accuracy;
if cfg.smoothWindow > 1
    for s = 1:size(acc, 1)
        acc(s, :) = smooth(acc(s, :), cfg.smoothWindow);
    end
end

fallback = oneDimClusterPermutation(acc, cfg.chance, cfg);
method = 'sign-flip max-cluster fallback';
fieldTripStat = [];

if cfg.fieldTripClusterUsable
    try
        fieldTripStat = runFieldTripCluster(acc, analysis.timeMs, cfg);
        [fallback.significantMask, fallback.clusters] = clustersFromFieldTrip(fieldTripStat, analysis.timeMs, cfg, fallback);
        method = 'FieldTrip ft_timelockstatistics';
    catch ME
        warning('FieldTrip cluster test failed for %s %s: %s. Using fallback.', ...
            analysis.name, group.dataset, ME.message);
    end
end

timepointTable = makeTimepointTable(acc, analysis, group, cfg);
clusterTable = makeClusterTable(fallback.clusters, acc, analysis, group, cfg, method);

result = struct();
result.analysis = analysis.name;
result.epoch = analysis.epoch;
result.condition = analysis.condition;
result.dataset = group.dataset;
result.subjectIds = group.subjectIds;
result.timeMs = analysis.timeMs;
result.accuracy = acc;
result.tObserved = fallback.tObserved;
result.significantMask = fallback.significantMask;
result.clusterThresholdT = fallback.clusterThresholdT;
result.maxClusterNull = fallback.maxClusterNull;
result.clusters = fallback.clusters;
result.method = method;
result.fieldTripStat = fieldTripStat;
result.timepointTable = timepointTable;
result.clusterTable = clusterTable;
end

function usable = checkFieldTripClusterUsable(cfg)
usable = false;
if ~cfg.preferFieldTrip || exist('ft_timelockstatistics', 'file') ~= 2
    return;
end
try
    tmpCfg = cfg;
    tmpCfg.nPerm = min(20, cfg.nPerm);
    smokeAcc = cfg.chance + [1:12; 2:13; 1:12; 2:13; 1:12; 2:13];
    runFieldTripCluster(smokeAcc, 0:11, tmpCfg);
    usable = true;
catch ME
    warning('FieldTrip cluster smoke test failed: %s. Built-in sign-flip cluster fallback will be used.', ME.message);
end
end

function stat = runFieldTripCluster(acc, timeMs, cfg)
if exist('ft_defaults', 'file') == 2
    ft_defaults;
end

nSub = size(acc, 1);
observed = [];
observed.label = {'decoding'};
observed.time = timeMs ./ 1000;
observed.dimord = 'subj_chan_time';
observed.individual = reshape(acc, nSub, 1, []);

chance = observed;
chance.individual = reshape(cfg.chance .* ones(size(acc)), nSub, 1, []);

design = [1:nSub 1:nSub; ones(1, nSub) 2 .* ones(1, nSub)];

ftCfg = [];
ftCfg.method = 'montecarlo';
ftCfg.statistic = 'depsamplesT';
ftCfg.correctm = 'cluster';
ftCfg.clusteralpha = cfg.clusterAlpha;
ftCfg.clusterstatistic = 'maxsum';
ftCfg.minnbchan = 0;
ftCfg.neighbours = [];
ftCfg.tail = 1;
ftCfg.clustertail = 1;
ftCfg.alpha = cfg.alpha;
ftCfg.numrandomization = cfg.nPerm;
ftCfg.design = design;
ftCfg.uvar = 1;
ftCfg.ivar = 2;
ftCfg.parameter = 'individual';
ftCfg.channel = {'decoding'};
ftCfg.latency = 'all';

stat = ft_timelockstatistics(ftCfg, observed, chance);
end

function [sigMask, clusters] = clustersFromFieldTrip(stat, timeMs, cfg, fallback)
sigMask = false(1, numel(timeMs));
clusters = fallback.clusters;

if isfield(stat, 'mask') && ~isempty(stat.mask)
    sigMask = squeeze(stat.mask);
    sigMask = sigMask(:)';
end

if ~isfield(stat, 'posclusters') || isempty(stat.posclusters) || ~isfield(stat, 'posclusterslabelmat')
    return;
end

labelMat = squeeze(stat.posclusterslabelmat);
labelMat = labelMat(:)';
clusters = struct('startIdx', {}, 'endIdx', {}, 'startMs', {}, 'endMs', {}, ...
    'mass', {}, 'p', {}, 'significant', {});
for c = 1:numel(stat.posclusters)
    idx = find(labelMat == c);
    if isempty(idx)
        continue;
    end
    clusters(end + 1).startIdx = idx(1); %#ok<AGROW>
    clusters(end).endIdx = idx(end);
    clusters(end).startMs = timeMs(idx(1));
    clusters(end).endMs = timeMs(idx(end));
    if isfield(stat.posclusters(c), 'clusterstat')
        clusters(end).mass = stat.posclusters(c).clusterstat;
    else
        clusters(end).mass = sum(fallback.tObserved(idx));
    end
    clusters(end).p = stat.posclusters(c).prob;
    clusters(end).significant = stat.posclusters(c).prob < cfg.alpha;
end
end

function result = oneDimClusterPermutation(acc, chance, cfg)
delta = acc - chance;
validTime = all(isfinite(delta), 1);
delta(:, ~validTime) = nan;
[tObserved, df] = oneSampleT(delta);
if exist('tinv', 'file') == 2
    clusterThresholdT = tinv(1 - cfg.clusterAlpha, df);
else
    clusterThresholdT = 1.7;
end

clusters = findPositiveClusters(tObserved, clusterThresholdT);
nSub = size(delta, 1);
maxClusterNull = zeros(cfg.nPerm, 1);
for p = 1:cfg.nPerm
    signs = (rand(nSub, 1) > 0.5) .* 2 - 1;
    permDelta = delta .* signs;
    permT = oneSampleT(permDelta);
    permClusters = findPositiveClusters(permT, clusterThresholdT);
    if isempty(permClusters)
        maxClusterNull(p) = 0;
    else
        maxClusterNull(p) = max([permClusters.mass]);
    end
end

sigMask = false(1, size(acc, 2));
for c = 1:numel(clusters)
    clusters(c).p = (1 + sum(maxClusterNull >= clusters(c).mass)) ./ (cfg.nPerm + 1);
    clusters(c).significant = clusters(c).p < cfg.alpha;
    if clusters(c).significant
        sigMask(clusters(c).startIdx:clusters(c).endIdx) = true;
    end
end

result = struct();
result.tObserved = tObserved;
result.clusterThresholdT = clusterThresholdT;
result.maxClusterNull = maxClusterNull;
result.clusters = clusters;
result.significantMask = sigMask;
end

function [tVals, df] = oneSampleT(delta)
n = sum(isfinite(delta), 1);
mu = mean(delta, 1, 'omitnan');
sd = std(delta, 0, 1, 'omitnan');
tVals = mu ./ (sd ./ sqrt(n));
tVals(sd == 0 | n < 2) = nan;
df = max(n) - 1;
end

function clusters = findPositiveClusters(tVals, thresholdT)
above = tVals > thresholdT;
clusters = struct('startIdx', {}, 'endIdx', {}, 'startMs', {}, 'endMs', {}, ...
    'mass', {}, 'p', {}, 'significant', {});
if ~any(above)
    return;
end
d = diff([false above false]);
starts = find(d == 1);
ends = find(d == -1) - 1;
for c = 1:numel(starts)
    idx = starts(c):ends(c);
    clusters(c).startIdx = starts(c); %#ok<AGROW>
    clusters(c).endIdx = ends(c);
    clusters(c).startMs = nan;
    clusters(c).endMs = nan;
    clusters(c).mass = sum(tVals(idx), 'omitnan');
    clusters(c).p = nan;
    clusters(c).significant = false;
end
end

function tableOut = makeTimepointTable(acc, analysis, group, cfg)
rows = {};
for t = 1:numel(analysis.timeMs)
    st = oneSampleStats(acc(:, t), cfg.chance, cfg);
    rows(end + 1, :) = {analysis.name, analysis.epoch, analysis.condition, group.dataset, ...
        analysis.timeMs(t), st.n, st.mean, st.ciLow, st.ciHigh, ...
        st.t, st.df, st.pRight, st.pTwoSided, st.dz, st.hedgesGz}; %#ok<AGROW>
end
tableOut = cell2table(rows, 'VariableNames', ...
    {'Analysis', 'Epoch', 'Condition', 'Dataset', 'TimeMs', 'N', ...
    'MeanAccuracy', 'CI95_Low', 'CI95_High', 'T', 'DF', ...
    'P_RightTail_Uncorrected', 'P_TwoSided_Uncorrected', 'CohenDz', 'HedgesGz'});
end

function tableOut = makeClusterTable(clusters, acc, analysis, group, cfg, method)
rows = {};
for c = 1:numel(clusters)
    idx = clusters(c).startIdx:clusters(c).endIdx;
    st = oneSampleStats(mean(acc(:, idx), 2, 'omitnan'), cfg.chance, cfg);
    rows(end + 1, :) = {analysis.name, analysis.epoch, analysis.condition, group.dataset, ...
        c, analysis.timeMs(clusters(c).startIdx), analysis.timeMs(clusters(c).endIdx), ...
        numel(idx), clusters(c).mass, clusters(c).p, clusters(c).significant, ...
        st.n, st.mean, st.ciLow, st.ciHigh, st.dz, st.hedgesGz, method}; %#ok<AGROW>
end
if isempty(rows)
    rows = {analysis.name, analysis.epoch, analysis.condition, group.dataset, ...
        nan, nan, nan, 0, nan, nan, false, size(acc, 1), nan, nan, nan, nan, nan, method};
end
tableOut = cell2table(rows, 'VariableNames', ...
    {'Analysis', 'Epoch', 'Condition', 'Dataset', 'ClusterID', 'StartMs', 'EndMs', ...
    'NTimepoints', 'ClusterMass', 'ClusterP_FWER', 'Significant', 'N', ...
    'ClusterMeanAccuracy', 'ClusterCI95_Low', 'ClusterCI95_High', ...
    'ClusterCohenDz', 'ClusterHedgesGz', 'Method'});
end

function st = oneSampleStats(x, chance, cfg)
x = x(:);
x = x(isfinite(x));
st = struct('n', numel(x), 'mean', nan, 'ciLow', nan, 'ciHigh', nan, ...
    't', nan, 'df', numel(x) - 1, 'pRight', nan, 'pTwoSided', nan, ...
    'dz', nan, 'hedgesGz', nan);
if numel(x) < 2
    return;
end
d = x - chance;
sd = std(d, 0, 'omitnan');
if sd == 0 || ~isfinite(sd)
    return;
end
st.mean = mean(x, 'omitnan');
sem = std(x, 0, 'omitnan') ./ sqrt(numel(x));
meanDiff = mean(d, 'omitnan');
st.t = meanDiff ./ (sd ./ sqrt(numel(x)));
st.df = numel(x) - 1;
if exist('tcdf', 'file') == 2
    st.pRight = 1 - tcdf(st.t, st.df);
    st.pTwoSided = 2 .* tcdf(-abs(st.t), st.df);
end
if exist('tinv', 'file') == 2
    tcrit = tinv(1 - cfg.alpha / 2, st.df);
    st.ciLow = st.mean - tcrit .* sem;
    st.ciHigh = st.mean + tcrit .* sem;
end
st.dz = meanDiff ./ sd;
J = 1 - 3 ./ (4 .* st.df - 1);
st.hedgesGz = st.dz .* J;
end

function plotDecodingCluster(result, cfg)
fig = figure('Color', 'w', 'Position', [100 100 930 430]);
hold on;
meanAcc = mean(result.accuracy, 1, 'omitnan');
n = sum(isfinite(result.accuracy), 1);
sem = std(result.accuracy, 0, 1, 'omitnan') ./ sqrt(n);
fill([result.timeMs fliplr(result.timeMs)], ...
    [meanAcc + sem fliplr(meanAcc - sem)], [0.18 0.42 0.68], ...
    'FaceAlpha', 0.18, 'EdgeColor', 'none');
plot(result.timeMs, meanAcc, 'Color', [0.12 0.32 0.55], 'LineWidth', 2.2);
yline(cfg.chance, '--', 'Chance', 'Color', [0.25 0.25 0.25]);

if any(result.significantMask)
    yLimits = ylim;
    sigRuns = maskToRuns(result.significantMask);
    for r = 1:size(sigRuns, 1)
        x1 = result.timeMs(sigRuns(r, 1));
        x2 = result.timeMs(sigRuns(r, 2));
        plot([x1 x2], [yLimits(1) + 0.06 .* range(yLimits) yLimits(1) + 0.06 .* range(yLimits)], ...
            'k-', 'LineWidth', 4);
    end
end

xlabel('Time (ms)');
ylabel('Decoding accuracy (%)');
title(sprintf('%s %s %s (%s)', result.dataset, result.epoch, result.condition, result.method), ...
    'Interpreter', 'none');
box off;
set(gca, 'FontName', 'Arial', 'FontSize', 12, 'LineWidth', 1.1);

base = fullfile(cfg.outDir, sprintf('reviewer4_eeg_%s_%s_%s', ...
    result.analysis, result.dataset, result.condition));
saveFigure(fig, base);
close(fig);
end

function runs = maskToRuns(mask)
d = diff([false mask(:)' false]);
runs = [find(d == 1)' find(d == -1)' - 1];
end

function txt = makeResponseText(clusterStats, cfg)
sig = clusterStats(clusterStats.Significant == true, :);
txt = sprintf(['Reviewer 4 EEG inference draft:\n\n' ...
    'For each EEG time-resolved decoding analysis, subject-level accuracy curves were tested against %.0f%% chance using a one-dimensional nonparametric cluster framework over time. ' ...
    'At each time point we computed a one-sample t statistic on accuracy minus chance. Candidate positive clusters were defined by a cluster-forming threshold of one-tailed p < %.3f, and cluster mass was the sum of t values across contiguous time points. ' ...
    'Family-wise error was controlled by sign-flipping each participant''s accuracy-minus-chance curve and retaining the maximum positive cluster mass across time for each of %d permutations. ' ...
    'FieldTrip ft_timelockstatistics was detected and used when available; if FieldTrip failed for a given contrast, the script used the equivalent built-in one-dimensional sign-flip max-cluster implementation. ' ...
    'The exported tables also report pointwise descriptive statistics, 95%% confidence intervals, Cohen''s dz, and Hedges gz. ' ...
    'Across the generated analyses, %d cluster(s) survived FWER correction at alpha=%.2f.'], ...
    cfg.chance, cfg.clusterAlpha, cfg.nPerm, height(sig), cfg.alpha);
end

function saveFigure(fig, basePath)
if exist('exportgraphics', 'file') == 2
    exportgraphics(fig, [basePath '.png'], 'Resolution', 300);
    exportgraphics(fig, [basePath '.tiff'], 'Resolution', 300);
else
    print(fig, [basePath '.png'], '-dpng', '-r300');
    print(fig, [basePath '.tiff'], '-dtiff', '-r300');
end
savefig(fig, [basePath '.fig']);
end
