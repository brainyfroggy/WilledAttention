%% Reviewer 2: post-cue alpha power for Choice vs Instructed cues
%
% Rationale:
% The reviewer asks whether the willed/choice condition recruits additional
% control regions without a behavioral cost because arousal/effort increases.
% A direct EEG check is to compare cue-locked alpha modulation after choice
% cues and instructed cues. Lower alpha after the cue is commonly interpreted
% as stronger cortical engagement/arousal. This script computes induced
% post-cue alpha power, baseline-normalized to the pre-cue interval.

clear; clc; close all;

cfg = struct();
cfg.rootDir = fileparts(fileparts(fileparts(mfilename('fullpath'))));
cfg.fs = 250;
cfg.nfft = 1024;
cfg.alphaHz = [8 12];
cfg.baselineIdx = 1:125;       % -500 to 0 ms
cfg.cueSample = 126;           % 0 ms
cfg.winSamples = 125;          % 500 ms windows
cfg.stepSamples = 5;           % 20 ms steps
cfg.postStarts = 126:5:376;    % windows from 0-500 to 1000-1500 ms
cfg.summaryWindowMs = [300 1000];
cfg.removeERP = true;          % induced alpha, matching prior ERP-removal logic
cfg.fzIdx = 17;                % used in prior project scripts
cfg.ozIdx = 20;                % used in prior project scripts
cfg.outDir = fullfile(cfg.rootDir, 'Result Figures', 'Reviewer2_PostCueAlpha');
cfg.ufLocFile = fullfile(cfg.rootDir, 'brainproduct31.locs');

if ~exist(cfg.outDir, 'dir')
    mkdir(cfg.outDir);
end

fprintf('Post-cue alpha analysis for Reviewer 2\n');
fprintf('Root: %s\n', cfg.rootDir);
fprintf('ERP removal before FFT: %d\n', cfg.removeERP);

datasets = struct([]);
datasets(1).name = 'UF';
datasets(1).choiceDir = fullfile(cfg.rootDir, 'EEG', 'UF', 'EEG processed data', 'Filtered', 'Epoched ICA', 'Choice');
datasets(1).instructedDir = fullfile(cfg.rootDir, 'EEG', 'UF', 'EEG processed data', 'Filtered', 'Epoched ICA', 'Instructed');
datasets(2).name = 'UCD';
datasets(2).choiceDir = fullfile(cfg.rootDir, 'EEG', 'UCD', 'EEG processed data', 'Filtered', 'Epoched ICA', 'Choice');
datasets(2).instructedDir = fullfile(cfg.rootDir, 'EEG', 'UCD', 'EEG processed data', 'Filtered', 'Epoched ICA', 'Instructed');

datasetResults = analyzeDataset(datasets(1), cfg);
for d = 2:numel(datasets)
    datasetResults(d) = analyzeDataset(datasets(d), cfg);
end

combined = combineDatasets(datasetResults, cfg);
allGroups = [datasetResults, combined];

statsWindow = makeSummaryStats(allGroups, cfg);
statsCurve = makeTimepointStats(allGroups, cfg);

save(fullfile(cfg.outDir, 'reviewer2_postcue_alpha_choice_vs_instructed.mat'), ...
    'cfg', 'datasets', 'datasetResults', 'combined', 'statsWindow', 'statsCurve');
writetable(statsWindow, fullfile(cfg.outDir, 'reviewer2_alpha_summary_window_stats.csv'));
writetable(statsCurve, fullfile(cfg.outDir, 'reviewer2_alpha_timepoint_stats.csv'));

for d = 1:numel(datasetResults)
    plotAlphaTimecourse(datasetResults(d), cfg, fullfile(cfg.outDir, [datasetResults(d).name '_alpha_timecourse']));
    plotAlphaDifference(datasetResults(d), statsCurve, cfg, fullfile(cfg.outDir, [datasetResults(d).name '_alpha_difference']));
    plotSummaryWindow(datasetResults(d), statsWindow, cfg, fullfile(cfg.outDir, [datasetResults(d).name '_alpha_300to1000ms']));
end

plotAlphaTimecourse(combined, cfg, fullfile(cfg.outDir, 'combined_alpha_timecourse'));
plotAlphaDifference(combined, statsCurve, cfg, fullfile(cfg.outDir, 'combined_alpha_difference'));
plotSummaryWindow(combined, statsWindow, cfg, fullfile(cfg.outDir, 'combined_alpha_300to1000ms'));
plotUFTopography(datasetResults, cfg, fullfile(cfg.outDir, 'UF_alpha_topography_choice_minus_instructed'));

fprintf('\nDone. Results saved to:\n%s\n', cfg.outDir);
disp(statsWindow);

%% Local functions

function dataset = analyzeDataset(dataset, cfg)
choiceFiles = dir(fullfile(dataset.choiceDir, '*.mat'));
instructedFiles = dir(fullfile(dataset.instructedDir, '*.mat'));

choiceNames = {choiceFiles.name};
instructedNames = {instructedFiles.name};
[subjectNames, choiceIdx, instructedIdx] = intersect(choiceNames, instructedNames, 'stable');

if isempty(subjectNames)
    error('No matched Choice/Instructed files found for %s.', dataset.name);
end

timeMs = windowCenterMs(cfg);
nSub = numel(subjectNames);
nTime = numel(timeMs);

metricNames = {'global', 'fz', 'oz'};
for m = 1:numel(metricNames)
    dataset.choice.(metricNames{m}) = nan(nSub, nTime);
    dataset.instructed.(metricNames{m}) = nan(nSub, nTime);
    dataset.choiceBase.(metricNames{m}) = nan(nSub, 1);
    dataset.instructedBase.(metricNames{m}) = nan(nSub, 1);
end

dataset.choiceTopo = {};
dataset.instructedTopo = {};
dataset.subjectNames = subjectNames(:);
dataset.timeMs = timeMs;

fprintf('\n%s: %d matched subjects\n', dataset.name, nSub);
for s = 1:nSub
    fprintf('  %s (%d/%d)\n', subjectNames{s}, s, nSub);

    choiceData = load(fullfile(dataset.choiceDir, choiceFiles(choiceIdx(s)).name));
    instructedData = load(fullfile(dataset.instructedDir, instructedFiles(instructedIdx(s)).name));

    choiceEpochs = catTrials(choiceData, 'choice_left', 'choice_right');
    instructedEpochs = catTrials(instructedData, 'instructed_left', 'instructed_right');

    [choiceDb, ~, choiceBase] = computeAlphaCondition(choiceEpochs, cfg);
    [instructedDb, ~, instructedBase] = computeAlphaCondition(instructedEpochs, cfg);

    dataset.choiceTopo{s, 1} = choiceDb;
    dataset.instructedTopo{s, 1} = instructedDb;

    dataset.choice.global(s, :) = mean(choiceDb, 1, 'omitnan');
    dataset.instructed.global(s, :) = mean(instructedDb, 1, 'omitnan');
    dataset.choice.fz(s, :) = channelTrace(choiceDb, cfg.fzIdx);
    dataset.instructed.fz(s, :) = channelTrace(instructedDb, cfg.fzIdx);
    dataset.choice.oz(s, :) = channelTrace(choiceDb, cfg.ozIdx);
    dataset.instructed.oz(s, :) = channelTrace(instructedDb, cfg.ozIdx);

    dataset.choiceBase.global(s, :) = mean(choiceBase, 1, 'omitnan');
    dataset.instructedBase.global(s, :) = mean(instructedBase, 1, 'omitnan');
    dataset.choiceBase.fz(s, :) = channelTrace(choiceBase, cfg.fzIdx);
    dataset.instructedBase.fz(s, :) = channelTrace(instructedBase, cfg.fzIdx);
    dataset.choiceBase.oz(s, :) = channelTrace(choiceBase, cfg.ozIdx);
    dataset.instructedBase.oz(s, :) = channelTrace(instructedBase, cfg.ozIdx);
end
end

function epochs = catTrials(dataStruct, leftName, rightName)
if ~isfield(dataStruct, leftName) || ~isfield(dataStruct, rightName)
    error('Missing expected variables %s and/or %s.', leftName, rightName);
end
epochs = cat(3, double(dataStruct.(leftName)), double(dataStruct.(rightName)));
end

function [alphaDb, alphaRaw, alphaBase] = computeAlphaCondition(data, cfg)
if cfg.removeERP
    data = data - mean(data, 3, 'omitnan');
end

alphaBaseByTrial = alphaInWindow(data, cfg.baselineIdx, cfg);
alphaBase = mean(alphaBaseByTrial, 2, 'omitnan');

nChan = size(data, 1);
nWin = numel(cfg.postStarts);
alphaDb = nan(nChan, nWin);
alphaRaw = nan(nChan, nWin);

for w = 1:nWin
    idx = cfg.postStarts(w):(cfg.postStarts(w) + cfg.winSamples - 1);
    alphaThisByTrial = alphaInWindow(data, idx, cfg);
    alphaRaw(:, w) = mean(alphaThisByTrial, 2, 'omitnan');
    alphaDb(:, w) = mean(10 .* log10(alphaThisByTrial ./ max(alphaBaseByTrial, eps)), 2, 'omitnan');
end
end

function alphaByTrial = alphaInWindow(data, idx, cfg)
window = hamming(numel(idx), 'periodic')';
x = data(:, idx, :) .* reshape(window, 1, [], 1);
y = fft(x, cfg.nfft, 2) ./ numel(idx);
powerSpectrum = (2 .* abs(y(:, 1:(cfg.nfft / 2 + 1), :))).^2;
freq = cfg.fs / 2 .* linspace(0, 1, cfg.nfft / 2 + 1);
alphaIdx = freq >= cfg.alphaHz(1) & freq <= cfg.alphaHz(2);
alphaByTrial = squeeze(mean(powerSpectrum(:, alphaIdx, :), 2, 'omitnan'));
if isvector(alphaByTrial)
    alphaByTrial = reshape(alphaByTrial, size(data, 1), []);
end
end

function trace = channelTrace(data, channelIdx)
if size(data, 1) >= channelIdx
    trace = data(channelIdx, :);
else
    trace = nan(1, size(data, 2));
end
end

function timeMs = windowCenterMs(cfg)
windowCenters = cfg.postStarts + (cfg.winSamples - 1) / 2;
timeMs = (windowCenters - cfg.cueSample) ./ cfg.fs .* 1000;
end

function combined = combineDatasets(datasets, cfg)
combined = struct();
combined.name = 'Combined';
combined.choiceDir = '';
combined.instructedDir = '';
combined.timeMs = datasets(1).timeMs;
combined.subjectNames = {};

metricNames = {'global', 'fz', 'oz'};
for m = 1:numel(metricNames)
    metric = metricNames{m};
    combined.choice.(metric) = [];
    combined.instructed.(metric) = [];
    combined.choiceBase.(metric) = [];
    combined.instructedBase.(metric) = [];
end

for d = 1:numel(datasets)
    for s = 1:numel(datasets(d).subjectNames)
        combined.subjectNames{end + 1, 1} = sprintf('%s_%s', datasets(d).name, datasets(d).subjectNames{s}); %#ok<AGROW>
    end
    for m = 1:numel(metricNames)
        metric = metricNames{m};
        combined.choice.(metric) = [combined.choice.(metric); datasets(d).choice.(metric)]; %#ok<AGROW>
        combined.instructed.(metric) = [combined.instructed.(metric); datasets(d).instructed.(metric)]; %#ok<AGROW>
        combined.choiceBase.(metric) = [combined.choiceBase.(metric); datasets(d).choiceBase.(metric)]; %#ok<AGROW>
        combined.instructedBase.(metric) = [combined.instructedBase.(metric); datasets(d).instructedBase.(metric)]; %#ok<AGROW>
    end
end

combined.choiceTopo = {};
combined.instructedTopo = {};
end

function statsTable = makeSummaryStats(groups, cfg)
metricNames = {'global', 'fz', 'oz'};
metricLabels = {'All channels', 'Fz', 'Oz'};
rows = {};

for g = 1:numel(groups)
    inWin = groups(g).timeMs >= cfg.summaryWindowMs(1) & groups(g).timeMs <= cfg.summaryWindowMs(2);
    for m = 1:numel(metricNames)
        metric = metricNames{m};
        choice = mean(groups(g).choice.(metric)(:, inWin), 2, 'omitnan');
        instructed = mean(groups(g).instructed.(metric)(:, inWin), 2, 'omitnan');
        st = pairedStats(choice, instructed);
        rows(end + 1, :) = {groups(g).name, metricLabels{m}, numel(st.diff), ...
            mean(choice, 'omitnan'), mean(instructed, 'omitnan'), mean(st.diff, 'omitnan'), ...
            st.t, st.df, st.p, st.dz}; %#ok<AGROW>
    end
end

statsTable = cell2table(rows, 'VariableNames', ...
    {'Dataset', 'Metric', 'N', 'ChoiceMeanDb', 'InstructedMeanDb', ...
    'ChoiceMinusInstructedDb', 'T', 'DF', 'P', 'CohenDz'});
statsTable.P_FDR = fdrBH(statsTable.P);
end

function statsTable = makeTimepointStats(groups, cfg)
metricNames = {'global', 'fz', 'oz'};
metricLabels = {'All channels', 'Fz', 'Oz'};
rows = {};

for g = 1:numel(groups)
    for m = 1:numel(metricNames)
        metric = metricNames{m};
        p = nan(numel(groups(g).timeMs), 1);
        for t = 1:numel(groups(g).timeMs)
            choice = groups(g).choice.(metric)(:, t);
            instructed = groups(g).instructed.(metric)(:, t);
            st = pairedStats(choice, instructed);
            p(t) = st.p;
            rows(end + 1, :) = {groups(g).name, metricLabels{m}, groups(g).timeMs(t), ...
                mean(choice, 'omitnan'), mean(instructed, 'omitnan'), ...
                mean(st.diff, 'omitnan'), st.t, st.df, st.p, st.dz, nan}; %#ok<AGROW>
        end
        q = fdrBH(p);
        startIdx = size(rows, 1) - numel(groups(g).timeMs) + 1;
        for t = 1:numel(groups(g).timeMs)
            rows{startIdx + t - 1, 11} = q(t);
        end
    end
end

statsTable = cell2table(rows, 'VariableNames', ...
    {'Dataset', 'Metric', 'TimeMs', 'ChoiceMeanDb', 'InstructedMeanDb', ...
    'ChoiceMinusInstructedDb', 'T', 'DF', 'P', 'CohenDz', 'P_FDR'});
end

function st = pairedStats(choice, instructed)
valid = isfinite(choice) & isfinite(instructed);
choice = choice(valid);
instructed = instructed(valid);
diffVals = choice - instructed;
st.diff = diffVals;
st.df = numel(diffVals) - 1;
st.dz = mean(diffVals, 'omitnan') ./ std(diffVals, 0, 'omitnan');
if numel(diffVals) < 2 || std(diffVals, 0, 'omitnan') == 0
    st.t = nan;
    st.p = nan;
else
    [~, st.p, ~, stats] = ttest(choice, instructed);
    st.t = stats.tstat;
    st.df = stats.df;
end
end

function q = fdrBH(p)
p = p(:);
q = nan(size(p));
valid = isfinite(p);
pv = p(valid);
if isempty(pv)
    return;
end
[ps, order] = sort(pv);
m = numel(ps);
qs = ps .* m ./ (1:m)';
qs = flipud(cummin(flipud(qs)));
qs(qs > 1) = 1;
tmp = nan(size(pv));
tmp(order) = qs;
q(valid) = tmp;
end

function plotAlphaTimecourse(group, cfg, outBase)
metricNames = {'global', 'fz', 'oz'};
metricLabels = {'All channels', 'Fz', 'Oz'};
choiceColor = [0.74 0.20 0.16];
instructedColor = [0.08 0.35 0.68];

figure('Color', 'w', 'Position', [100 100 850 900]);
tiledlayout(3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
for m = 1:numel(metricNames)
    nexttile;
    metric = metricNames{m};
    hold on;
    hChoice = plotMeanSem(group.timeMs, group.choice.(metric), choiceColor);
    hInstructed = plotMeanSem(group.timeMs, group.instructed.(metric), instructedColor);
    yline(0, ':', 'Color', [0.35 0.35 0.35]);
    xline(cfg.summaryWindowMs(1), ':', 'Color', [0.50 0.50 0.50]);
    xline(cfg.summaryWindowMs(2), ':', 'Color', [0.50 0.50 0.50]);
    title(metricLabels{m});
    ylabel('Alpha change (dB)');
    set(gca, 'Box', 'off', 'FontName', 'Arial', 'FontSize', 12, 'LineWidth', 1.2);
    if m == 1
        legend([hChoice hInstructed], {'Choice', 'Instructed'}, 'Location', 'best', 'Box', 'off');
    end
    if m == numel(metricNames)
        xlabel('Time after cue (ms)');
    end
end
exportgraphics(gcf, [outBase '.png'], 'Resolution', 300);
exportgraphics(gcf, [outBase '.tiff'], 'Resolution', 300);
end

function plotAlphaDifference(group, statsCurve, cfg, outBase)
metricNames = {'global', 'fz', 'oz'};
metricLabels = {'All channels', 'Fz', 'Oz'};
diffColor = [0.25 0.25 0.25];

figure('Color', 'w', 'Position', [120 120 850 900]);
tiledlayout(3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
for m = 1:numel(metricNames)
    nexttile;
    metric = metricNames{m};
    diffData = group.choice.(metric) - group.instructed.(metric);
    hold on;
    plotMeanSem(group.timeMs, diffData, diffColor);
    yline(0, '-', 'Color', [0.25 0.25 0.25]);
    xline(cfg.summaryWindowMs(1), ':', 'Color', [0.50 0.50 0.50]);
    xline(cfg.summaryWindowMs(2), ':', 'Color', [0.50 0.50 0.50]);

    rowIdx = strcmp(statsCurve.Dataset, group.name) & strcmp(statsCurve.Metric, metricLabels{m}) & statsCurve.P_FDR < 0.05;
    sigTimes = statsCurve.TimeMs(rowIdx);
    if ~isempty(sigTimes)
        yLimits = ylim;
        scatter(sigTimes, repmat(yLimits(1) + 0.06 * range(yLimits), size(sigTimes)), 18, ...
            'filled', 'MarkerFaceColor', [0.10 0.10 0.10], 'MarkerEdgeColor', 'none');
    end

    title([metricLabels{m} ': Choice - Instructed']);
    ylabel('\Delta alpha (dB)');
    set(gca, 'Box', 'off', 'FontName', 'Arial', 'FontSize', 12, 'LineWidth', 1.2);
    if m == numel(metricNames)
        xlabel('Time after cue (ms)');
    end
end
exportgraphics(gcf, [outBase '.png'], 'Resolution', 300);
exportgraphics(gcf, [outBase '.tiff'], 'Resolution', 300);
end

function plotSummaryWindow(group, statsWindow, cfg, outBase)
metricNames = {'global', 'fz', 'oz'};
metricLabels = {'All channels', 'Fz', 'Oz'};
inWin = group.timeMs >= cfg.summaryWindowMs(1) & group.timeMs <= cfg.summaryWindowMs(2);
choiceColor = [0.74 0.20 0.16];
instructedColor = [0.08 0.35 0.68];

figure('Color', 'w', 'Position', [150 150 900 340]);
tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
for m = 1:numel(metricNames)
    nexttile;
    metric = metricNames{m};
    choice = mean(group.choice.(metric)(:, inWin), 2, 'omitnan');
    instructed = mean(group.instructed.(metric)(:, inWin), 2, 'omitnan');
    hold on;
    for s = 1:numel(choice)
        plot([1 2], [choice(s) instructed(s)], '-', 'Color', [0.78 0.78 0.78], 'LineWidth', 0.75);
    end
    scatter(ones(size(choice)) - 0.06, choice, 18, 'filled', 'MarkerFaceColor', choiceColor, 'MarkerFaceAlpha', 0.65);
    scatter(2 .* ones(size(instructed)) + 0.06, instructed, 18, 'filled', 'MarkerFaceColor', instructedColor, 'MarkerFaceAlpha', 0.65);
    errorbar([1 2], [mean(choice, 'omitnan') mean(instructed, 'omitnan')], ...
        [std(choice, 0, 'omitnan') ./ sqrt(sum(isfinite(choice))) std(instructed, 0, 'omitnan') ./ sqrt(sum(isfinite(instructed)))], ...
        'k', 'LineStyle', 'none', 'LineWidth', 1.5, 'CapSize', 0);
    plot([1 2], [mean(choice, 'omitnan') mean(instructed, 'omitnan')], 'k-', 'LineWidth', 2.2);
    xlim([0.5 2.5]);
    xticks([1 2]);
    xticklabels({'Choice', 'Instructed'});
    ylabel('Alpha change (dB)');
    title(metricLabels{m});
    row = statsWindow(strcmp(statsWindow.Dataset, group.name) & strcmp(statsWindow.Metric, metricLabels{m}), :);
    if ~isempty(row)
        yLimits = ylim;
        text(1.5, yLimits(2) - 0.08 * range(yLimits), sprintf('p = %.3g, dz = %.2f', row.P, row.CohenDz), ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'FontSize', 10);
    end
    set(gca, 'Box', 'off', 'FontName', 'Arial', 'FontSize', 12, 'LineWidth', 1.2);
end
exportgraphics(gcf, [outBase '.png'], 'Resolution', 300);
exportgraphics(gcf, [outBase '.tiff'], 'Resolution', 300);
end

function lineHandle = plotMeanSem(x, data, color)
meanData = mean(data, 1, 'omitnan');
n = sum(isfinite(data), 1);
semData = std(data, 0, 1, 'omitnan') ./ sqrt(n);
fill([x fliplr(x)], [meanData + semData fliplr(meanData - semData)], color, ...
    'FaceAlpha', 0.18, 'EdgeColor', 'none');
lineHandle = plot(x, meanData, 'Color', color, 'LineWidth', 2.5);
end

function plotUFTopography(datasets, cfg, outBase)
ufIdx = find(strcmp({datasets.name}, 'UF'), 1);
if isempty(ufIdx) || exist('topoplot', 'file') ~= 2 || ~exist(cfg.ufLocFile, 'file')
    fprintf('Skipping UF topography: topoplot or %s not available on MATLAB path.\n', cfg.ufLocFile);
    return;
end

uf = datasets(ufIdx);
inWin = uf.timeMs >= cfg.summaryWindowMs(1) & uf.timeMs <= cfg.summaryWindowMs(2);
nSub = numel(uf.subjectNames);
nChan = size(uf.choiceTopo{1}, 1);
diffTopo = nan(nSub, nChan);

for s = 1:nSub
    diffTopo(s, :) = mean(uf.choiceTopo{s}(:, inWin) - uf.instructedTopo{s}(:, inWin), 2, 'omitnan');
end

figure('Color', 'w', 'Position', [200 200 450 380]);
topoplot(mean(diffTopo, 1, 'omitnan'), cfg.ufLocFile, ...
    'style', 'map', 'electrodes', 'on', 'plotrad', 0.84, 'headrad', 0.75);
title('UF Choice - Instructed alpha, 300-1000 ms');
colorbar;
colormap('jet');
exportgraphics(gcf, [outBase '.png'], 'Resolution', 300);
exportgraphics(gcf, [outBase '.tiff'], 'Resolution', 300);
end
