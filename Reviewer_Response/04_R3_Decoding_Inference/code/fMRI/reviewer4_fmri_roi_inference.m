%% Reviewer 4: fMRI ROI decoding inference
%
% Adds explicit group-level inference for ROI decoding accuracy without
% modifying the original decoding scripts. The original fMRI scripts save
% subject x CV-repeat x ROI accuracies; this script averages repeats within
% subject, tests subject-level accuracies against chance, reports CIs and
% within-subject effect sizes, and corrects across ROIs.

clear; clc; close all;

cfg = struct();
cfg.rootDir = fileparts(fileparts(fileparts(mfilename('fullpath'))));
cfg.roiSet = '7 Frontal';
cfg.conditions = {'Choice', 'Instructed'};
cfg.datasets = {'UF', 'UCD'};
cfg.chance = 50;
cfg.alpha = 0.05;
cfg.nPerm = 2000;          % increase to 5000 or 10000 for final manuscript
cfg.nBoot = 5000;
cfg.randomSeed = 20260705;
cfg.includeCombinedRoiInCorrection = true;
cfg.generateSubjectPermutationNull = false; % prototype support; expensive
cfg.outDir = fullfile(cfg.rootDir, 'Result Figures', 'Reviewer4_DecodingInference', 'fMRI');

if strcmp(getenv('REVIEWER4_FAST'), '1')
    cfg.nPerm = 100;
    cfg.nBoot = 500;
end

if ~exist(cfg.outDir, 'dir')
    mkdir(cfg.outDir);
end

rng(cfg.randomSeed, 'twister');

fprintf('Reviewer 4 fMRI ROI decoding inference\n');
fprintf('Root: %s\n', cfg.rootDir);
fprintf('ROI set: %s\n', cfg.roiSet);
fprintf('Output: %s\n', cfg.outDir);

roiLabels = getRoiLabels(cfg.rootDir, cfg.datasets{1}, cfg.roiSet);
datasetResults = loadAllResults(cfg, roiLabels);
summaryTable = makeSummaryTable(datasetResults, cfg);
subjectTable = makeSubjectTable(datasetResults, cfg);
responseText = makeResponseText(summaryTable, cfg);

writetable(summaryTable, fullfile(cfg.outDir, 'reviewer4_fmri_roi_group_stats.csv'));
writetable(subjectTable, fullfile(cfg.outDir, 'reviewer4_fmri_roi_subject_accuracies.csv'));

fid = fopen(fullfile(cfg.outDir, 'reviewer4_fmri_response_text.txt'), 'w');
fprintf(fid, '%s\n', responseText);
fclose(fid);

save(fullfile(cfg.outDir, 'reviewer4_fmri_roi_inference.mat'), ...
    'cfg', 'roiLabels', 'datasetResults', 'summaryTable', 'subjectTable', 'responseText');

plotRoiInference(datasetResults, summaryTable, cfg);

if cfg.generateSubjectPermutationNull
    fprintf('\nGenerating subject-level label-shuffle nulls. This can be slow.\n');
    nullResults = generateFmriSubjectPermutationNull(cfg, roiLabels); %#ok<UNRCH>
    save(fullfile(cfg.outDir, 'reviewer4_fmri_subject_permutation_null.mat'), ...
        'cfg', 'roiLabels', 'nullResults', '-v7.3');
end

fprintf('\nDone. Results saved to:\n%s\n', cfg.outDir);
disp(summaryTable);

%% Local functions

function labels = getRoiLabels(rootDir, dataset, roiSet)
roiDir = fullfile(rootDir, 'fMRI', 'ROI', dataset, roiSet);
roiFiles = dir(fullfile(roiDir, '*.mat'));
if isempty(roiFiles)
    error('No ROI files found: %s', roiDir);
end
labels = erase({roiFiles.name}, '.mat');
labels = matlab.lang.makeUniqueStrings(labels);
end

function results = loadAllResults(cfg, roiLabels)
template = struct('dataset', '', 'condition', '', 'roiLabels', {{}}, ...
    'accuracy', [], 'subjectIds', {{}}, 'sourceFile', '');
results = repmat(template, 0, 1);
idx = 0;

for d = 1:numel(cfg.datasets)
    for c = 1:numel(cfg.conditions)
        idx = idx + 1;
        results(idx) = loadOneResult(cfg.datasets{d}, cfg.conditions{c}, cfg, roiLabels); %#ok<AGROW>
    end
end

for c = 1:numel(cfg.conditions)
    condition = cfg.conditions{c};
    idx = idx + 1;
    results(idx).dataset = 'Combined';
    results(idx).condition = condition;
    results(idx).roiLabels = roiLabels;
    results(idx).subjectIds = {};
    results(idx).accuracy = [];
    results(idx).sourceFile = '';
    for d = 1:numel(cfg.datasets)
        one = results(strcmp({results.dataset}, cfg.datasets{d}) & strcmp({results.condition}, condition));
        results(idx).accuracy = [results(idx).accuracy; one.accuracy]; %#ok<AGROW>
        for s = 1:numel(one.subjectIds)
            results(idx).subjectIds{end + 1, 1} = sprintf('%s_%s', cfg.datasets{d}, one.subjectIds{s}); %#ok<AGROW>
        end
    end
end
end

function result = loadOneResult(dataset, condition, cfg, roiLabels)
matPath = fullfile(cfg.rootDir, 'fMRI', 'fMRI_BetaDecoding', dataset, cfg.roiSet, ...
    sprintf('BetaDecoding_zscored_%s.mat', condition));
if ~exist(matPath, 'file')
    error('Missing fMRI decoding file: %s', matPath);
end

S = load(matPath, 'AccuracyInst');
if ~isfield(S, 'AccuracyInst')
    error('Missing AccuracyInst in %s', matPath);
end

acc = double(S.AccuracyInst);
if ndims(acc) == 4
    acc = squeeze(acc(:, :, :, 1));
end
if ndims(acc) ~= 3
    error('Expected AccuracyInst as subject x repeat x ROI in %s.', matPath);
end
acc = squeeze(mean(acc, 2, 'omitnan'));
if isvector(acc)
    acc = acc(:);
end

if size(acc, 2) ~= numel(roiLabels)
    warning('ROI label count (%d) differs from accuracy columns (%d) for %s.', ...
        numel(roiLabels), size(acc, 2), matPath);
    roiLabels = makeFallbackRoiLabels(size(acc, 2), roiLabels);
end

result = struct();
result.dataset = dataset;
result.condition = condition;
result.roiLabels = roiLabels;
result.accuracy = acc;
result.subjectIds = arrayfun(@(x) sprintf('S%02d', x), 1:size(acc, 1), 'UniformOutput', false)';
result.sourceFile = matPath;
end

function labels = makeFallbackRoiLabels(n, labels)
labels = labels(:)';
if numel(labels) >= n
    labels = labels(1:n);
else
    for i = (numel(labels) + 1):n
        labels{i} = sprintf('ROI_%02d', i); %#ok<AGROW>
    end
end
end

function summaryTable = makeSummaryTable(results, cfg)
rows = {};

for r = 1:numel(results)
    acc = results(r).accuracy;
    roiLabels = results(r).roiLabels;
    pRight = nan(1, size(acc, 2));
    tVals = nan(1, size(acc, 2));

    roiMask = true(1, size(acc, 2));
    if ~cfg.includeCombinedRoiInCorrection
        roiMask = ~contains(lower(roiLabels), {'whole', 'combined'});
    end

    for roi = 1:size(acc, 2)
        st = oneSampleStats(acc(:, roi), cfg.chance, cfg);
        pRight(roi) = st.pRight;
        tVals(roi) = st.t;
        rows(end + 1, :) = {results(r).dataset, results(r).condition, roiLabels{roi}, ...
            st.n, st.mean, st.sd, st.sem, st.ciLow, st.ciHigh, st.bootCiLow, st.bootCiHigh, ...
            st.t, st.df, st.pRight, st.pTwoSided, st.dz, st.hedgesGz, nan, nan}; %#ok<AGROW>
    end

    q = nan(size(pRight));
    q(roiMask) = fdrBH(pRight(roiMask));
    pMax = nan(size(pRight));
    pMax(roiMask) = maxStatPermutation(acc(:, roiMask), cfg.chance, tVals(roiMask), cfg);

    startIdx = size(rows, 1) - size(acc, 2) + 1;
    for roi = 1:size(acc, 2)
        rows{startIdx + roi - 1, 18} = q(roi);
        rows{startIdx + roi - 1, 19} = pMax(roi);
    end
end

summaryTable = cell2table(rows, 'VariableNames', ...
    {'Dataset', 'Condition', 'ROI', 'N', 'MeanAccuracy', 'SD', 'SEM', ...
    'CI95_Low', 'CI95_High', 'BootstrapCI95_Low', 'BootstrapCI95_High', ...
    'T', 'DF', 'P_RightTail', 'P_TwoSided', 'CohenDz', 'HedgesGz', ...
    'P_FDR_RightTail', 'P_MaxStat_RightTail'});
end

function subjectTable = makeSubjectTable(results, cfg)
rows = {};
for r = 1:numel(results)
    for s = 1:size(results(r).accuracy, 1)
        for roi = 1:size(results(r).accuracy, 2)
            rows(end + 1, :) = {results(r).dataset, results(r).condition, ...
                results(r).subjectIds{s}, results(r).roiLabels{roi}, ...
                results(r).accuracy(s, roi), results(r).accuracy(s, roi) - cfg.chance}; %#ok<AGROW>
        end
    end
end
subjectTable = cell2table(rows, 'VariableNames', ...
    {'Dataset', 'Condition', 'Subject', 'ROI', 'Accuracy', 'AccuracyMinusChance'});
end

function st = oneSampleStats(x, chance, cfg)
x = x(:);
x = x(isfinite(x));
st = struct('n', numel(x), 'mean', nan, 'sd', nan, 'sem', nan, ...
    'ciLow', nan, 'ciHigh', nan, 'bootCiLow', nan, 'bootCiHigh', nan, ...
    't', nan, 'df', numel(x) - 1, 'pRight', nan, 'pTwoSided', nan, ...
    'dz', nan, 'hedgesGz', nan);
if numel(x) < 2
    return;
end

d = x - chance;
st.mean = mean(x, 'omitnan');
st.sd = std(x, 0, 'omitnan');
st.sem = st.sd ./ sqrt(st.n);
sdDiff = std(d, 0, 'omitnan');
if sdDiff == 0 || ~isfinite(sdDiff)
    return;
end

meanDiff = mean(d, 'omitnan');
st.t = meanDiff ./ (sdDiff ./ sqrt(st.n));
st.df = st.n - 1;
if exist('tcdf', 'file') == 2
    st.pRight = 1 - tcdf(st.t, st.df);
    st.pTwoSided = 2 .* tcdf(-abs(st.t), st.df);
end
if exist('tinv', 'file') == 2
    tcrit = tinv(1 - cfg.alpha / 2, st.df);
    st.ciLow = st.mean - tcrit .* st.sem;
    st.ciHigh = st.mean + tcrit .* st.sem;
end
st.dz = meanDiff ./ sdDiff;
J = 1 - 3 ./ (4 .* st.df - 1);
st.hedgesGz = st.dz .* J;
if cfg.nBoot > 0
    bootMean = nan(cfg.nBoot, 1);
    for b = 1:cfg.nBoot
        bootMean(b) = mean(x(randi(st.n, st.n, 1)), 'omitnan');
    end
    st.bootCiLow = prctile(bootMean, 2.5);
    st.bootCiHigh = prctile(bootMean, 97.5);
end
end

function pCorr = maxStatPermutation(acc, chance, observedT, cfg)
delta = acc - chance;
validRoi = all(isfinite(delta), 1);
delta = delta(:, validRoi);
observedTValid = observedT(validRoi);
pCorr = nan(size(observedT));
if isempty(delta) || size(delta, 1) < 2
    return;
end

nSub = size(delta, 1);
maxT = nan(cfg.nPerm, 1);
for p = 1:cfg.nPerm
    signs = (rand(nSub, 1) > 0.5) .* 2 - 1;
    permDelta = delta .* signs;
    permT = mean(permDelta, 1, 'omitnan') ./ (std(permDelta, 0, 1, 'omitnan') ./ sqrt(nSub));
    maxT(p) = max(permT);
end
tmp = (1 + sum(maxT >= observedTValid, 1)) ./ (cfg.nPerm + 1);
pCorr(validRoi) = tmp;
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

function plotRoiInference(results, summaryTable, cfg)
for r = 1:numel(results)
    fig = figure('Color', 'w', 'Position', [80 80 1150 470]);
    acc = results(r).accuracy;
    roiLabels = results(r).roiLabels;
    meanAcc = mean(acc, 1, 'omitnan');
    ciLow = nan(size(meanAcc));
    ciHigh = nan(size(meanAcc));
    for roi = 1:numel(roiLabels)
        row = summaryTable(strcmp(summaryTable.Dataset, results(r).dataset) & ...
            strcmp(summaryTable.Condition, results(r).condition) & ...
            strcmp(summaryTable.ROI, roiLabels{roi}), :);
        ciLow(roi) = row.CI95_Low;
        ciHigh(roi) = row.CI95_High;
    end

    hold on;
    bar(1:numel(roiLabels), meanAcc, 0.72, 'FaceColor', [0.72 0.78 0.86], 'EdgeColor', 'none');
    errorbar(1:numel(roiLabels), meanAcc, meanAcc - ciLow, ciHigh - meanAcc, ...
        'k', 'LineStyle', 'none', 'LineWidth', 1.3, 'CapSize', 0);
    for roi = 1:numel(roiLabels)
        x = roi + linspace(-0.16, 0.16, size(acc, 1))';
        scatter(x, acc(:, roi), 24, 'filled', 'MarkerFaceColor', [0.20 0.36 0.55], ...
            'MarkerFaceAlpha', 0.55, 'MarkerEdgeColor', 'none');
    end
    yline(cfg.chance, '--', 'Chance', 'Color', [0.25 0.25 0.25], 'LineWidth', 1.2);
    ylabel('Decoding accuracy (%)');
    title(sprintf('%s %s fMRI ROI decoding', results(r).dataset, results(r).condition));
    xticks(1:numel(roiLabels));
    xticklabels(strrep(roiLabels, '_', '\_'));
    xtickangle(35);
    box off;
    set(gca, 'FontName', 'Arial', 'FontSize', 11, 'LineWidth', 1.1);
    ylim([min([cfg.chance - 10, acc(:)']) max([cfg.chance + 25, acc(:)'])]);

    sigRow = summaryTable(strcmp(summaryTable.Dataset, results(r).dataset) & ...
        strcmp(summaryTable.Condition, results(r).condition) & ...
        summaryTable.P_MaxStat_RightTail < cfg.alpha, :);
    for i = 1:height(sigRow)
        roi = find(strcmp(roiLabels, sigRow.ROI{i}));
        text(roi, max(ylim) - 0.06 .* range(ylim), '*', 'HorizontalAlignment', 'center', ...
            'FontSize', 18, 'FontWeight', 'bold');
    end

    base = fullfile(cfg.outDir, sprintf('reviewer4_fmri_%s_%s_roi_accuracy', ...
        results(r).dataset, results(r).condition));
    saveFigure(fig, base);
    close(fig);
end
end

function txt = makeResponseText(summaryTable, cfg)
combinedChoice = summaryTable(strcmp(summaryTable.Dataset, 'Combined') & ...
    strcmp(summaryTable.Condition, 'Choice'), :);
nSigFdr = sum(combinedChoice.P_FDR_RightTail < cfg.alpha, 'omitnan');
nSigMax = sum(combinedChoice.P_MaxStat_RightTail < cfg.alpha, 'omitnan');
txt = sprintf(['Reviewer 4 fMRI inference draft:\n\n' ...
    'For each ROI and participant, decoding accuracy was first averaged across the 100 cross-validation repetitions produced by the original beta-series SVM pipeline. ' ...
    'Group inference was then performed on these subject-level accuracies with one-sample tests against %.0f%% chance. ' ...
    'We report mean accuracy, 95%% t-based confidence intervals, bootstrap confidence intervals for the mean, Cohen''s dz, and small-sample-corrected Hedges gz. ' ...
    'Because multiple ROIs were tested, p-values were corrected within each dataset x condition family using Benjamini-Hochberg FDR and a sign-flip max-statistic permutation correction across ROIs (%d permutations). ' ...
    'In the combined-sample Choice analysis, %d ROIs survived FDR correction and %d ROIs survived max-statistic correction at alpha=%.2f. ' ...
    'No subject-level fMRI label-shuffle null files were found with the original outputs; therefore the primary implemented analysis uses the nominal binary-classification chance level of 50%%. ' ...
    'The script also contains an optional label-shuffling prototype to generate empirical subject-level null accuracies from the beta-series data and ROI masks if the final analysis requires empirical chance estimates.'], ...
    cfg.chance, cfg.nPerm, nSigFdr, nSigMax, cfg.alpha);
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

function nullResults = generateFmriSubjectPermutationNull(cfg, roiLabels)
% Prototype: rerun the same linear-SVM cross-validation after shuffling
% left/right labels within subject. This function is intentionally opt-in
% because it can be slow and requires libsvm or fitcsvm on the MATLAB path.
nullResults = struct([]);
for d = 1:numel(cfg.datasets)
    dataset = cfg.datasets{d};
    roiVoxels = loadUsefulVoxels(cfg, dataset);
    for c = 1:numel(cfg.conditions)
        condition = cfg.conditions{c};
        dataDir = fullfile(cfg.rootDir, 'fMRI', 'fMRI BetaSeries', [dataset ' Combined']);
        files = dir(fullfile(dataDir, '*.mat'));
        accNull = nan(numel(files), cfg.nPerm, numel(roiLabels));
        for s = 1:numel(files)
            data = load(fullfile(dataDir, files(s).name));
            [leftData, rightData] = conditionData(data, dataset, condition);
            for roi = 1:numel(roiLabels)
                X = [leftData(:, roiVoxels{roi}); rightData(:, roiVoxels{roi})];
                y = [-ones(size(leftData, 1), 1); ones(size(rightData, 1), 1)];
                X = zscore(X, 0, 2);
                X(~isfinite(X)) = 0;
                for p = 1:cfg.nPerm
                    accNull(s, p, roi) = crossValidatedAccuracy(X, y(randperm(numel(y))), 10);
                end
            end
        end
        idx = numel(nullResults) + 1;
        nullResults(idx).dataset = dataset;
        nullResults(idx).condition = condition;
        nullResults(idx).roiLabels = roiLabels;
        nullResults(idx).subjectFiles = {files.name}';
        nullResults(idx).accuracyNull = accNull;
    end
end
end

function roiVoxels = loadUsefulVoxels(cfg, dataset)
usefulFile = fullfile(cfg.rootDir, 'fMRI', sprintf('final_useful_voxels_%s.mat', dataset));
S = load(usefulFile, 'final_useful_voxels');
roiDir = fullfile(cfg.rootDir, 'fMRI', 'ROI', dataset, cfg.roiSet);
roiFiles = dir(fullfile(roiDir, '*.mat'));
roiVoxels = cell(numel(roiFiles), 1);
for r = 1:numel(roiFiles)
    R = load(fullfile(roiDir, roiFiles(r).name), 'roi');
    vox = find(R.roi > 0)';
    roiVoxels{r} = vox(ismember(vox, S.final_useful_voxels));
end
end

function [leftData, rightData] = conditionData(data, dataset, condition)
switch dataset
    case 'UF'
        switch condition
            case 'Choice'
                leftData = data.choice_left;
                rightData = data.choice_right;
            otherwise
                leftData = data.instructed_left;
                rightData = data.instructed_right;
        end
    case 'UCD'
        switch condition
            case 'Choice'
                leftData = data.BetaSeriesChooseLeft;
                rightData = data.BetaSeriesChooseRight;
            otherwise
                leftData = data.BetaSeriesAttendLeft;
                rightData = data.BetaSeriesAttendRight;
        end
end
end

function acc = crossValidatedAccuracy(X, y, kFold)
foldId = stratifiedFoldId(y, min(kFold, min(histcounts(categorical(y)))));
pred = nan(size(y));
for k = 1:max(foldId)
    test = foldId == k;
    train = ~test;
    pred(test) = trainPredictLinearSvm(X(train, :), y(train), X(test, :));
end
acc = 100 .* mean(pred == y, 'omitnan');
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
    model = fitcsvm(xTrain, yTrain, 'KernelFunction', 'linear', ...
        'BoxConstraint', 1, 'Standardize', false);
    pred = predict(model, xTest);
elseif exist('svmtrain', 'file') == 2 && exist('svmpredict', 'file') == 2
    model = svmtrain(yTrain, double(xTrain), '-q -t 0 -c 1'); %#ok<SVMTRAIN>
    pred = svmpredict(zeros(size(xTest, 1), 1), double(xTest), model, '-q'); %#ok<SVMPREDICT>
else
    error('No supported SVM function found. Install Statistics and Machine Learning Toolbox or libsvm.');
end
pred = pred(:);
end
