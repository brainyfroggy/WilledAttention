%% Plot saved non-CSD (voltage) decoding results only
% This script does not rerun decoding. It loads the saved
% Decoding_Voltage__*.mat files and regenerates the reviewer-facing plot.

clear; clc; close all;

SCRIPT_DIR = fileparts(mfilename('fullpath'));
if isempty(SCRIPT_DIR); SCRIPT_DIR = pwd; end
OUTDIR = fullfile(SCRIPT_DIR,'results_posterior');

caseLabels = {'UF_Choice','UF_Instructed','UCD_Choice','UCD_Instructed'};
plotLabels = {'UF Choice','UF Instructed','UCD Choice','UCD Instructed'};
cols = [0.75 0.75 0.75; 0.20 0.45 0.70; 0.85 0.45 0.25]; % all / posterior / frontal
YL = [45 60];

M = nan(numel(caseLabels),3);
E = nan(numel(caseLabels),3);
P = nan(numel(caseLabels),3);

rows = {};
for j = 1:numel(caseLabels)
    fpath = fullfile(OUTDIR, sprintf('Decoding_Voltage__%s.mat', caseLabels{j}));
    assert(exist(fpath,'file')==2, 'Saved voltage result not found: %s', fpath);
    S = load(fpath);

    sets = {S.subjAccAll, S.subjAccPosterior, S.subjAccFront};
    for s = 1:3
        a = sets{s};
        M(j,s) = mean(a);
        E(j,s) = std(a) / sqrt(numel(a));
        [~,P(j,s)] = ttest(a,50);
    end

    [~,pPostAll] = ttest(S.subjAccPosterior, S.subjAccAll);
    [~,pPostFront] = ttest(S.subjAccPosterior, S.subjAccFront);
    [~,pFrontAll] = ttest(S.subjAccFront, S.subjAccAll);
    rows(end+1,:) = {S.site, S.cond, S.nChanAll, numel(S.idxPost), numel(S.idxFront), ...
        M(j,1), P(j,1), M(j,2), P(j,2), M(j,3), P(j,3), ...
        M(j,2)-M(j,1), pPostAll, M(j,2)-M(j,3), pPostFront, M(j,3)-M(j,1), pFrontAll}; %#ok<SAGROW>
end

figure('Color','w','Position',[80 80 900 520]);
hold on;
X = 1:numel(caseLabels);
b = bar(X, M, 'grouped');
for s = 1:3
    b(s).FaceColor = cols(s,:);
end

for s = 1:3
    xc = b(s).XEndPoints;
    errorbar(xc, M(:,s), E(:,s), 'k', 'LineStyle', 'none', 'CapSize', 6);
    for j = 1:numel(caseLabels)
        starY = min(M(j,s)+E(j,s)+0.35, YL(2)-0.15);
        text(xc(j), starY, local_star(P(j,s)), ...
            'HorizontalAlignment','center', 'FontSize', 10, 'FontWeight','bold');
    end
end

set(gca, 'XTick', X, 'XTickLabel', plotLabels);
ylabel('Decoding accuracy (%)');
ylim(YL);
legend(b, {'All channels','Posterior ROI','Frontal ROI'}, 'Location','northeast');
title('Pre-cue alpha decoding');
box off;

saveas(gcf, fullfile(OUTDIR,'voltage_all_vs_posterior_vs_frontal_y45_60.png'));
saveas(gcf, fullfile(OUTDIR,'voltage_all_vs_posterior_vs_frontal_y45_60.fig'));

T = cell2table(rows, 'VariableNames', {'site','cond','nChanAll','nChanPost','nChanFront', ...
    'accAll_pct','pAll_vs50','accPost_pct','pPost_vs50','accFront_pct','pFront_vs50', ...
    'deltaPostMinusAll_pct','pPost_vsAll','deltaPostMinusFront_pct','pPost_vsFront', ...
    'deltaFrontMinusAll_pct','pFront_vsAll'});
writetable(T, fullfile(OUTDIR,'accuracy_summary_voltage_only.csv'));

fprintf('Saved voltage-only plot and CSV to:\n  %s\n', OUTDIR);

function s = local_star(p)
    if     p < 0.001, s = '***';
    elseif p < 0.01,  s = '**';
    elseif p < 0.05,  s = '*';
    else,             s = 'n.s.';
    end
end
