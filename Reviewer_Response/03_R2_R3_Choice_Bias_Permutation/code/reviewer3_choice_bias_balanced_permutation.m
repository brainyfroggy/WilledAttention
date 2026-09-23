%% Reviewer 3: choice bias, balanced decoding, and permutation baseline
%
% Reviewer concern:
% If subjects choose Left and Right unequally often in willed/choice trials,
% an unbalanced left-vs-right classifier can exceed 50% by learning the
% majority class. This analysis quantifies left/right choice proportions,
% reruns choice decoding after balancing trial counts within subject, and
% estimates empirical chance by permuting labels within subject.
%
% Decoding features:
% Pre-cue voltage alpha power, matching the existing project scripts:
% EEG/<dataset>/EEG processed data/Filtered/Voltage PreCue Alpha/Choice.
% These files contain Pxx_choice_left_final and Pxx_choice_right_final,
% computed from -500 to 0 ms and summed across the 8-12 Hz FFT bins.

clear; clc; close all;

cfg = struct();
cfg.rootDir = fileparts(fileparts(fileparts(mfilename('fullpath'))));
cfg.outDir = fullfile(cfg.rootDir, 'Result Figures', 'Reviewer3_ChoiceBias');
cfg.randomSeed = 20260705;

% Practical defaults. For the final manuscript run, set cfg.nPerm to 1000
% or 10000. If runtime is high, reduce cfg.nObservedRepeats first.
cfg.nPerm = 10000;
cfg.nObservedRepeats = 20;
cfg.nPermRepeats = 1;
cfg.kFold = 10;
cfg.runUnbalancedObserved = true;

% Optional environment override for a quick smoke test from the command line:
% set REVIEWER3_FAST=1 before launching MATLAB.
if strcmp(getenv('REVIEWER3_FAST'), '1')
    cfg.nPerm = 10;
    cfg.nObservedRepeats = 5;
    cfg.nPermRepeats = 1;
end

if ~exist(cfg.outDir, 'dir')
    mkdir(cfg.outDir);
end

rng(cfg.randomSeed, 'twister');

datasets = struct([]);
datasets(1).name = 'UF';
datasets(1).choiceEpochDir = fullfile(cfg.rootDir, 'EEG', 'UF', 'EEG processed data', 'Filtered', 'Epoched ICA', 'Choice');
datasets(1).choiceAlphaDir = fullfile(cfg.rootDir, 'EEG', 'UF', 'EEG processed data', 'Filtered', 'Voltage PreCue Alpha', 'Choice');
datasets(2).name = 'UCD';
datasets(2).choiceEpochDir = fullfile(cfg.rootDir, 'EEG', 'UCD', 'EEG processed data', 'Filtered', 'Epoched ICA', 'Choice');
datasets(2).choiceAlphaDir = fullfile(cfg.rootDir, 'EEG', 'UCD', 'EEG processed data', 'Filtered', 'Voltage PreCue Alpha', 'Choice');

fprintf('Reviewer 3 choice-bias analysis\n');
fprintf('Root: %s\n', cfg.rootDir);
fprintf('Output: %s\n', cfg.outDir);
fprintf('Permutations: %d\n', cfg.nPerm);

datasetResults = analyzeDataset(datasets(1), cfg);
for d = 2:numel(datasets)
    datasetResults(d) = analyzeDataset(datasets(d), cfg); %#ok<SAGROW>
end

[biasTable, decodingTable, permutationTable] = makeSubjectTables(datasetResults);
groupStatsTable = makeGroupStats(datasetResults, cfg);
groupPermutationTable = makeGroupPermutationTable(datasetResults, cfg);
responseText = makeResponseText(groupStatsTable, decodingTable, cfg);

writetable(biasTable, fullfile(cfg.outDir, 'reviewer3_choice_bias_by_subject.csv'));
writetable(decodingTable, fullfile(cfg.outDir, 'reviewer3_balanced_decoding_by_subject.csv'));
writetable(permutationTable, fullfile(cfg.outDir, 'reviewer3_balanced_permutation_by_subject.csv'));
writetable(groupStatsTable, fullfile(cfg.outDir, 'reviewer3_group_stats.csv'));
writetable(groupPermutationTable, fullfile(cfg.outDir, 'reviewer3_group_permutation_distribution.csv'));

fid = fopen(fullfile(cfg.outDir, 'reviewer3_response_text.txt'), 'w');
fprintf(fid, '%s\n', responseText);
fclose(fid);

save(fullfile(cfg.outDir, 'reviewer3_choice_bias_balanced_permutation.mat'), ...
    'cfg', 'datasets', 'datasetResults', 'biasTable', 'decodingTable', ...
    'permutationTable', 'groupStatsTable', 'groupPermutationTable', 'responseText');

plotChoiceBias(biasTable, cfg);
plotPermutationBaselines(groupPermutationTable, groupStatsTable, cfg);
plotSubjectDeltas(decodingTable, cfg);

fprintf('\nDone. Results saved to:\n%s\n\n', cfg.outDir);
disp(groupStatsTable);
fprintf('\nDraft response text saved as reviewer3_response_text.txt\n');

%% Local functions

function result = analyzeDataset(dataset, cfg)
epochFiles = dir(fullfile(dataset.choiceEpochDir, '*.mat'));
if isempty(epochFiles)
    error('No choice epoch files found for %s: %s', dataset.name, dataset.choiceEpochDir);
end

result = struct();
result.name = dataset.name;
result.subject = {};
result.fileName = {};
result.nLeft = [];
result.nRight = [];
result.leftProp = [];
result.leftPct = [];
result.majorityPct = [];
result.binomialP = [];
result.unbalancedAcc = [];
result.balancedAcc = [];
result.balancedPerm = nan(numel(epochFiles), cfg.nPerm);
result.subjectPermP = [];
result.decodeIncluded = false(numel(epochFiles), 1);

fprintf('\n%s: %d choice epoch files\n', dataset.name, numel(epochFiles));

for s = 1:numel(epochFiles)
    fileName = epochFiles(s).name;
    [~, subjectName] = fileparts(fileName);
    fprintf('  %s (%d/%d)\n', subjectName, s, numel(epochFiles));

    epochData = load(fullfile(dataset.choiceEpochDir, fileName), 'choice_left', 'choice_right');
    if ~isfield(epochData, 'choice_left') || ~isfield(epochData, 'choice_right')
        warning('Skipping %s %s: missing choice_left and/or choice_right.', dataset.name, fileName);
        continue;
    end

    nLeft = size(epochData.choice_left, 3);
    nRight = size(epochData.choice_right, 3);
    nTotal = nLeft + nRight;
    leftProp = nLeft ./ nTotal;

    result.subject{s, 1} = subjectName; %#ok<AGROW>
    result.fileName{s, 1} = fileName; %#ok<AGROW>
    result.nLeft(s, 1) = nLeft;
    result.nRight(s, 1) = nRight;
    result.leftProp(s, 1) = leftProp;
    result.leftPct(s, 1) = 100 .* leftProp;
    result.majorityPct(s, 1) = 100 .* max(nLeft, nRight) ./ nTotal;
    result.binomialP(s, 1) = twoSidedBinomialP(nLeft, nTotal, 0.5);
    result.unbalancedAcc(s, 1) = nan;
    result.balancedAcc(s, 1) = nan;
    result.subjectPermP(s, 1) = nan;

    alphaPath = fullfile(dataset.choiceAlphaDir, fileName);
    if ~exist(alphaPath, 'file')
        warning('No pre-cue alpha feature file for %s %s. Bias retained; decoding set to NaN.', dataset.name, fileName);
        continue;
    end

    alphaData = load(alphaPath, 'Pxx_choice_left_final', 'Pxx_choice_right_final');
    if ~isfield(alphaData, 'Pxx_choice_left_final') || ~isfield(alphaData, 'Pxx_choice_right_final')
        warning('Skipping decoding for %s %s: missing Pxx_choice_*_final.', dataset.name, fileName);
        continue;
    end

    leftFeatures = trialFeatures(alphaData.Pxx_choice_left_final);
    rightFeatures = trialFeatures(alphaData.Pxx_choice_right_final);

    if cfg.runUnbalancedObserved
        result.unbalancedAcc(s, 1) = repeatedDecoding(leftFeatures, rightFeatures, cfg, ...
            cfg.nObservedRepeats, false, false);
    end

    result.balancedAcc(s, 1) = repeatedDecoding(leftFeatures, rightFeatures, cfg, ...
        cfg.nObservedRepeats, true, false);

    permAcc = nan(1, cfg.nPerm);
    for p = 1:cfg.nPerm
        permAcc(p) = repeatedDecoding(leftFeatures, rightFeatures, cfg, ...
            cfg.nPermRepeats, true, true);
    end

    result.balancedPerm(s, :) = permAcc;
    result.subjectPermP(s, 1) = (1 + sum(permAcc >= result.balancedAcc(s), 'omitnan')) ./ ...
        (1 + sum(~isnan(permAcc)));
    result.decodeIncluded(s, 1) = true;
end

keep = ~cellfun(@isempty, result.subject);
result.subject = result.subject(keep);
result.fileName = result.fileName(keep);
result.nLeft = result.nLeft(keep);
result.nRight = result.nRight(keep);
result.leftProp = result.leftProp(keep);
result.leftPct = result.leftPct(keep);
result.majorityPct = result.majorityPct(keep);
result.binomialP = result.binomialP(keep);
result.unbalancedAcc = result.unbalancedAcc(keep);
result.balancedAcc = result.balancedAcc(keep);
result.balancedPerm = result.balancedPerm(keep, :);
result.subjectPermP = result.subjectPermP(keep);
result.decodeIncluded = result.decodeIncluded(keep);
end

function X = trialFeatures(data)
nTrial = size(data, 3);
X = permute(double(data), [3 1 2]);
X = reshape(X, nTrial, []);
badColumn = all(isnan(X), 1) | all(isinf(X), 1);
X(:, badColumn) = [];
X(~isfinite(X)) = nan;
colMean = mean(X, 1, 'omitnan');
for c = 1:size(X, 2)
    bad = isnan(X(:, c));
    X(bad, c) = colMean(c);
end
X(~isfinite(X)) = 0;
end

function acc = repeatedDecoding(leftFeatures, rightFeatures, cfg, nRepeats, doBalance, doPermute)
accEach = nan(nRepeats, 1);
for r = 1:nRepeats
    [X, y] = makeDesign(leftFeatures, rightFeatures, doBalance);
    if isempty(X) || numel(unique(y)) < 2
        continue;
    end
    if doPermute
        y = y(randperm(numel(y)));
    end
    accEach(r) = crossValidatedAccuracy(X, y, cfg.kFold);
end
acc = mean(accEach, 'omitnan');
end

function [X, y] = makeDesign(leftFeatures, rightFeatures, doBalance)
nLeft = size(leftFeatures, 1);
nRight = size(rightFeatures, 1);

if doBalance
    nMin = min(nLeft, nRight);
    if nMin < 2
        X = [];
        y = [];
        return;
    end
    leftIdx = randperm(nLeft, nMin);
    rightIdx = randperm(nRight, nMin);
    X = [leftFeatures(leftIdx, :); rightFeatures(rightIdx, :)];
    y = [-ones(nMin, 1); ones(nMin, 1)];
else
    X = [leftFeatures; rightFeatures];
    y = [-ones(nLeft, 1); ones(nRight, 1)];
end

X = normalizeRows(X);
end

function X = normalizeRows(X)
rowMean = mean(X, 2, 'omitnan');
rowSd = std(X, 0, 2, 'omitnan');
rowSd(rowSd == 0 | isnan(rowSd)) = 1;
X = (X - rowMean) ./ rowSd;
X(~isfinite(X)) = 0;
end

function acc = crossValidatedAccuracy(X, y, kFold)
minClassN = min(arrayfun(@(lab) sum(y == lab), unique(y)));
kFold = min(kFold, minClassN);
if kFold < 2
    acc = nan;
    return;
end

foldId = stratifiedFoldId(y, kFold);
predAll = nan(size(y));

for k = 1:kFold
    test = foldId == k;
    train = ~test;
    predAll(test) = trainPredictLinearSvm(X(train, :), y(train), X(test, :));
end

acc = 100 .* mean(predAll == y, 'omitnan');
end

function foldId = stratifiedFoldId(y, kFold)
foldId = zeros(size(y));
labels = unique(y(:))';
for lab = labels
    idx = find(y == lab);
    idx = idx(randperm(numel(idx)));
    foldId(idx) = mod(0:(numel(idx) - 1), kFold) + 1;
end
end

function pred = trainPredictLinearSvm(xTrain, yTrain, xTest)
if exist('fitcsvm', 'file') == 2
    model = fitcsvm(xTrain, yTrain, ...
        'KernelFunction', 'linear', ...
        'BoxConstraint', 1, ...
        'Standardize', false);
    pred = predict(model, xTest);
elseif exist('svmtrain', 'file') == 2 && exist('svmpredict', 'file') == 2
    model = svmtrain(yTrain, double(xTrain), '-q -t 0 -c 1'); %#ok<SVMTRAIN>
    pred = svmpredict(zeros(size(xTest, 1), 1), double(xTest), model, '-q'); %#ok<SVMPREDICT>
else
    error('No supported SVM function found. Install Statistics and Machine Learning Toolbox or libsvm.');
end
pred = pred(:);
end

function p = twoSidedBinomialP(k, n, p0)
if n == 0
    p = nan;
    return;
end
allK = 0:n;
logProb = gammaln(n + 1) - gammaln(allK + 1) - gammaln(n - allK + 1) + ...
    allK .* log(p0) + (n - allK) .* log(1 - p0);
prob = exp(logProb);
observedProb = prob(k + 1);
p = min(1, sum(prob(prob <= observedProb + eps)));
end

function [biasTable, decodingTable, permutationTable] = makeSubjectTables(datasetResults)
biasRows = {};
decodeRows = {};
permRows = {};

for d = 1:numel(datasetResults)
    res = datasetResults(d);
    for s = 1:numel(res.subject)
        biasRows(end + 1, :) = {res.name, res.subject{s}, res.fileName{s}, ...
            res.nLeft(s), res.nRight(s), res.nLeft(s) + res.nRight(s), ...
            res.leftProp(s), res.leftPct(s), res.majorityPct(s), res.binomialP(s)}; %#ok<AGROW>

        permMean = mean(res.balancedPerm(s, :), 'omitnan');
        permSd = std(res.balancedPerm(s, :), 0, 'omitnan');
        decodeRows(end + 1, :) = {res.name, res.subject{s}, res.fileName{s}, ...
            res.decodeIncluded(s), res.unbalancedAcc(s), res.balancedAcc(s), ...
            permMean, permSd, res.subjectPermP(s), res.balancedAcc(s) - permMean, ...
            res.majorityPct(s)}; %#ok<AGROW>

        for p = 1:size(res.balancedPerm, 2)
            permRows(end + 1, :) = {res.name, res.subject{s}, p, res.balancedPerm(s, p)}; %#ok<AGROW>
        end
    end
end

biasTable = cell2table(biasRows, 'VariableNames', ...
    {'Dataset', 'Subject', 'FileName', 'NLeft', 'NRight', 'NTotal', ...
    'LeftProportion', 'LeftPercent', 'MajorityClassPercent', 'BinomialP'});

decodingTable = cell2table(decodeRows, 'VariableNames', ...
    {'Dataset', 'Subject', 'FileName', 'DecodeIncluded', 'UnbalancedObservedAccuracy', ...
    'BalancedObservedAccuracy', 'BalancedPermutationMean', 'BalancedPermutationSD', ...
    'SubjectPermutationP_ObservedGreater', 'BalancedMinusPermutationMean', ...
    'MajorityClassPercent'});

permutationTable = cell2table(permRows, 'VariableNames', ...
    {'Dataset', 'Subject', 'Permutation', 'BalancedPermutationAccuracy'});
end

function statsTable = makeGroupStats(datasetResults, cfg)
rows = {};
groupNames = [{datasetResults.name}, {'Combined'}];

for g = 1:numel(groupNames)
    if g <= numel(datasetResults)
        name = datasetResults(g).name;
        leftProp = datasetResults(g).leftProp;
        balancedAcc = datasetResults(g).balancedAcc;
        unbalancedAcc = datasetResults(g).unbalancedAcc;
        majorityPct = datasetResults(g).majorityPct;
        perm = datasetResults(g).balancedPerm;
    else
        name = 'Combined';
        leftProp = vertcat(datasetResults.leftProp);
        balancedAcc = vertcat(datasetResults.balancedAcc);
        unbalancedAcc = vertcat(datasetResults.unbalancedAcc);
        majorityPct = vertcat(datasetResults.majorityPct);
        perm = vertcat(datasetResults.balancedPerm);
    end

    stBias = oneSampleStats(leftProp, 0.5);
    rows(end + 1, :) = {name, 'Left choice proportion vs 0.5', numel(leftProp), ...
        mean(leftProp, 'omitnan'), 0.5, stBias.t, stBias.df, stBias.p, stBias.d, nan}; %#ok<AGROW>

    stBal = oneSampleStats(balancedAcc, 50);
    obsMean = mean(balancedAcc, 'omitnan');
    permGroupMean = mean(perm, 1, 'omitnan');
    groupPermP = (1 + sum(permGroupMean >= obsMean, 'omitnan')) ./ ...
        (1 + sum(~isnan(permGroupMean)));
    rows(end + 1, :) = {name, 'Balanced decoding accuracy vs 50', sum(~isnan(balancedAcc)), ...
        obsMean, 50, stBal.t, stBal.df, stBal.p, stBal.d, groupPermP}; %#ok<AGROW>

    stUnbal = oneSampleStats(unbalancedAcc, 50);
    rows(end + 1, :) = {name, 'Unbalanced observed accuracy vs 50', sum(~isnan(unbalancedAcc)), ...
        mean(unbalancedAcc, 'omitnan'), 50, stUnbal.t, stUnbal.df, stUnbal.p, stUnbal.d, nan}; %#ok<AGROW>

    stMajority = oneSampleStats(majorityPct, 50);
    rows(end + 1, :) = {name, 'Majority class baseline vs 50', numel(majorityPct), ...
        mean(majorityPct, 'omitnan'), 50, stMajority.t, stMajority.df, stMajority.p, stMajority.d, nan}; %#ok<AGROW>
end

statsTable = cell2table(rows, 'VariableNames', ...
    {'Dataset', 'Test', 'N', 'Mean', 'Reference', 'T', 'DF', 'P', 'CohenD', ...
    'GroupPermutationP_ObservedGreater'});
end

function st = oneSampleStats(x, mu)
x = x(:);
x = x(~isnan(x));
st = struct('t', nan, 'df', numel(x) - 1, 'p', nan, 'd', nan);
if numel(x) < 2
    return;
end
sd = std(x, 0);
if sd == 0
    return;
end
st.t = (mean(x) - mu) ./ (sd ./ sqrt(numel(x)));
st.df = numel(x) - 1;
if exist('tcdf', 'file') == 2
    st.p = 2 .* tcdf(-abs(st.t), st.df);
end
st.d = (mean(x) - mu) ./ sd;
end

function groupPermutationTable = makeGroupPermutationTable(datasetResults, cfg)
rows = {};
for d = 1:numel(datasetResults)
    permMean = mean(datasetResults(d).balancedPerm, 1, 'omitnan');
    for p = 1:cfg.nPerm
        rows(end + 1, :) = {datasetResults(d).name, p, permMean(p)}; %#ok<AGROW>
    end
end

combinedPerm = mean(vertcat(datasetResults.balancedPerm), 1, 'omitnan');
for p = 1:cfg.nPerm
    rows(end + 1, :) = {'Combined', p, combinedPerm(p)}; %#ok<AGROW>
end

groupPermutationTable = cell2table(rows, 'VariableNames', ...
    {'Dataset', 'Permutation', 'MeanBalancedPermutationAccuracy'});
end

function txt = makeResponseText(groupStatsTable, decodingTable, cfg)
combinedBias = groupStatsTable(strcmp(groupStatsTable.Dataset, 'Combined') & ...
    strcmp(groupStatsTable.Test, 'Left choice proportion vs 0.5'), :);
combinedBal = groupStatsTable(strcmp(groupStatsTable.Dataset, 'Combined') & ...
    strcmp(groupStatsTable.Test, 'Balanced decoding accuracy vs 50'), :);
combinedMajority = groupStatsTable(strcmp(groupStatsTable.Dataset, 'Combined') & ...
    strcmp(groupStatsTable.Test, 'Majority class baseline vs 50'), :);

strongBiasN = sum(abs(decodingTable.MajorityClassPercent - 50) >= 10);
decodeN = sum(decodingTable.DecodeIncluded);
if combinedBal.GroupPermutationP_ObservedGreater < 0.05
    permSentence = sprintf(['The group-level permutation p-value for observed balanced accuracy ' ...
        'exceeding the permuted baseline was %.4g.'], combinedBal.GroupPermutationP_ObservedGreater);
else
    permSentence = sprintf(['With cfg.nPerm=%d, the group-level permutation p-value for observed balanced accuracy ' ...
        'exceeding the permuted baseline was %.4g; this should be rerun with 1000 or 10000 permutations ' ...
        'for the final manuscript value.'], cfg.nPerm, combinedBal.GroupPermutationP_ObservedGreater);
end

txt = sprintf(['Reviewer 3 response draft:\n\n' ...
    'We agree that individual left/right choice imbalance could inflate an unbalanced classifier. ' ...
    'We therefore quantified each participant''s choice proportions and added two control analyses. ' ...
    'Across the combined EEG sample, the mean left-choice proportion was %.3f (one-sample t-test vs 0.5: t(%d)=%.2f, p=%.4g). ' ...
    'The mean subject-level majority-class baseline was %.2f%%, and %d participants had a majority-class percentage at least 10 percentage points away from 50%%.\n\n' ...
    'To remove the trivial majority-class strategy, we reran left-vs-right decoding after downsampling the larger choice class within each subject to match the smaller class. ' ...
    'The analysis used the same pre-cue voltage alpha power features as the original decoding (-500 to 0 ms, 8-12 Hz, linear SVM, stratified %d-fold cross-validation). ' ...
    'The balanced observed decoding mean was %.2f%% across %d decoded participants. ' ...
    'We also generated an empirical chance distribution by shuffling left/right labels within subject and rerunning the same balanced decoding pipeline %d times. %s ' ...
    'Because the balanced analysis uses equal numbers of left and right trials for every subject, the classifier cannot obtain above-chance performance by exploiting a subject-specific majority side. ' ...
    'The permutation test directly addresses the reviewer''s concern that nominal 50%% chance may be inappropriate under label imbalance.'], ...
    combinedBias.Mean, combinedBias.DF, combinedBias.T, combinedBias.P, ...
    combinedMajority.Mean, strongBiasN, cfg.kFold, combinedBal.Mean, decodeN, ...
    cfg.nPerm, permSentence);
end

function plotChoiceBias(biasTable, cfg)
fig = figure('Color', 'w', 'Position', [100 100 1100 450]);
tiledlayout(1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

nexttile;
hold on;
datasets = unique(biasTable.Dataset, 'stable');
colors = lines(numel(datasets));
for d = 1:numel(datasets)
    idx = strcmp(biasTable.Dataset, datasets{d});
    x = d + 0.12 .* randn(sum(idx), 1);
    scatter(x, biasTable.LeftPercent(idx), 46, colors(d, :), 'filled', ...
        'MarkerFaceAlpha', 0.75, 'MarkerEdgeColor', 'k');
end
yline(50, '--k', '50%');
set(gca, 'XTick', 1:numel(datasets), 'XTickLabel', datasets);
xlim([0.5 numel(datasets) + 0.5]);
ylabel('Left choices (%)');
title('Subject-level choice bias');
box off;

nexttile;
histogram(biasTable.LeftPercent, 12, 'FaceColor', [0.25 0.45 0.75], 'EdgeColor', 'w');
xline(50, '--k', '50%');
xlabel('Left choices (%)');
ylabel('Number of subjects');
title('Combined distribution');
box off;

saveFigure(fig, fullfile(cfg.outDir, 'reviewer3_choice_bias_subjects'));
close(fig);
end

function plotPermutationBaselines(groupPermutationTable, groupStatsTable, cfg)
datasets = unique(groupPermutationTable.Dataset, 'stable');
fig = figure('Color', 'w', 'Position', [100 100 1200 420]);
tiledlayout(1, numel(datasets), 'Padding', 'compact', 'TileSpacing', 'compact');
for d = 1:numel(datasets)
    nexttile;
    idx = strcmp(groupPermutationTable.Dataset, datasets{d});
    histogram(groupPermutationTable.MeanBalancedPermutationAccuracy(idx), 18, ...
        'FaceColor', [0.65 0.65 0.65], 'EdgeColor', 'w');
    obs = groupStatsTable.Mean(strcmp(groupStatsTable.Dataset, datasets{d}) & ...
        strcmp(groupStatsTable.Test, 'Balanced decoding accuracy vs 50'));
    xline(obs, '-r', sprintf('Observed %.2f%%', obs), 'LineWidth', 1.5);
    xline(50, '--k', '50%');
    xlabel('Group mean balanced accuracy (%)');
    ylabel('Permutations');
    title(datasets{d});
    box off;
end
saveFigure(fig, fullfile(cfg.outDir, 'reviewer3_group_permutation_baseline'));
close(fig);
end

function plotSubjectDeltas(decodingTable, cfg)
plotTable = decodingTable(decodingTable.DecodeIncluded, :);
[~, order] = sort(plotTable.BalancedMinusPermutationMean, 'descend');
plotTable = plotTable(order, :);

fig = figure('Color', 'w', 'Position', [100 100 1200 470]);
bar(plotTable.BalancedMinusPermutationMean, 'FaceColor', [0.25 0.55 0.75], 'EdgeColor', 'none');
yline(0, '--k');
xlabel('Subject');
ylabel('Observed balanced accuracy - permutation mean (%)');
title('Subject-level balanced decoding relative to empirical chance');
set(gca, 'XTick', 1:height(plotTable), 'XTickLabel', plotTable.Subject, 'XTickLabelRotation', 60);
box off;
saveFigure(fig, fullfile(cfg.outDir, 'reviewer3_subject_decoding_minus_permutation'));
close(fig);
end

function saveFigure(fig, basePath)
if exist('exportgraphics', 'file') == 2
    exportgraphics(fig, [basePath '.png'], 'Resolution', 300);
else
    print(fig, [basePath '.png'], '-dpng', '-r300');
end
savefig(fig, [basePath '.fig']);
end
