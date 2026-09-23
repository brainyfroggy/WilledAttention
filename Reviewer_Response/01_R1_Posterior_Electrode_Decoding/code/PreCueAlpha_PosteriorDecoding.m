%% =====================================================================
%  Pre-cue alpha-power decoding: all channels vs POSTERIOR and FRONTAL ROIs
%  -- response to Reviewer #1, Public Review, Comment 2.
%
%  Reviewer concern: "the EEG decoding approach utilized the entire
%  topography of electrodes rather than a biologically motivated posterior
%  region of interest. Given that alpha-mediated spatial attention is
%  traditionally localized to parieto-occipital sensors, using the full
%  electrode set risks the inclusion of non-neural artifacts such as
%  micro-saccades or muscle activity."
%
%  For each dataset/condition this script decodes attend-left vs attend-right
%  from pre-cue (-500..0 ms) alpha power using (a) ALL channels, (b) the
%  POSTERIOR parieto-occipital ROI, and (c) a FRONTAL control ROI, on the
%  SAME cross-validation folds, so the three are directly comparable. It
%  runs BOTH data types:
%     CSD     = current source density / surface Laplacian
%     Voltage = scalp voltage (this is what the paper's Figure 4A reported)
%  On the dense 58-ch UCD montage CSD sharpens topographies and raises
%  decoding ~6%% vs voltage; on the 31-ch UF montage the two are ~equal.
%
%  Pipeline matches the paper: z-score across channels, linear SVM (c=1),
%  100 x 10-fold CV.  Author: (Ding Lab), referee report 2-26-26.
%  Requires: Statistics & ML Toolbox (fitcsvm, crossvalind, ttest); LibSVM
%            optional (USE_LIBSVM=true).
%  =====================================================================

clear; clc; close all;

%% ---------------------------------------------------------------------
%  0.  USER PATHS
%  ---------------------------------------------------------------------
SCRIPT_DIR = fileparts(mfilename('fullpath'));
if isempty(SCRIPT_DIR); SCRIPT_DIR = pwd; end

% Changhao's shared data/code root.
DATA_ROOT = 'N:\Experimental_Data\Changhao Xiong\2022-2025\WilledAttention';

% Keep montage files with this reviewer-response analysis so the exact
% channel names/indices used below are documented and portable.
LAYOUT_DIR = fullfile(SCRIPT_DIR,'channel_layouts');
UF_LOCS   = fullfile(LAYOUT_DIR,'brainproduct31.locs'); % 31 ch, BrainAmp
UCD_CED   = fullfile(LAYOUT_DIR,'NeuroScan_58.ced');    % 58 ch, Neuroscan
UCD_DROP  = {};   % NeuroScan_58.ced already has exactly the 58 EEG channels

OUTDIR = fullfile(SCRIPT_DIR,'results_posterior');
if ~exist(OUTDIR,'dir'); mkdir(OUTDIR); end

assert(exist(DATA_ROOT,'dir')==7, 'Data root not found: %s', DATA_ROOT);
assert(exist(UF_LOCS,'file')==2, 'UF channel file not found: %s', UF_LOCS);
assert(exist(UCD_CED,'file')==2, 'UCD channel file not found: %s', UCD_CED);

%% ---------------------------------------------------------------------
%  1.  POSTERIOR (PARIETO-OCCIPITAL) ROI
%  ---------------------------------------------------------------------
%  Matching strips apostrophes + upper-cases, so old Neuroscan names in
%  NeuroScan_58.ced (P3', O1', T5', PzP, CB1, ...) match their plain forms.
posteriorLabels = upper({ ...
    'P7','P5','P3','P1','Pz','P2','P4','P6','P8','P9','P10', ...
    'PO7','PO5','PO3','PO1','POz','PO2','PO4','PO6','PO8','PO9','PO10', ...
    'O1','Oz','O2','Iz', ...
    'CP5','CP3','CP1','CPz','CP2','CP4','CP6','TP7','TP8','TP9','TP10', ...
    ... % --- Neuroscan-specific posterior labels (NeuroScan_58.ced) ---
    'T5','T6', 'CB1','CB2', 'PzA','PzP', ...
    'P1P','P2P','P3P','P4P', 'C1P','C2P','C3P','C4P', 'TCP1','TCP2'});

%  FRONTAL ROI (control): if the choice signal is genuinely posterior alpha,
%  frontal channels should decode much worse.  If frontal were instead doing
%  the work, that would point to ocular/muscle artifact -- the very thing the
%  reviewer worried about.  Frontal = Fp/AF/F/FC/FT line (+ Neuroscan F*A,C*A).
frontalLabels = upper({ ...
    'Fp1','Fp2','Fpz','AF7','AF3','AFz','AF4','AF8', ...
    'F7','F5','F3','F1','Fz','F2','F4','F6','F8', ...
    'FC5','FC3','FC1','FCz','FC2','FC4','FC6','FT7','FT8','FT9','FT10', ...
    ... % --- Neuroscan-specific frontal labels ---
    'F3A','F4A', 'C1A','C2A','C3A','C4A','CzA'});

%% ---------------------------------------------------------------------
%  2.  DATA TYPES x CASES
%  ---------------------------------------------------------------------
% {displayName , folder name under ...\Filtered\}
dataTypes = {
  'CSD',     'CSD PreCue Alpha'
  'Voltage', 'Voltage PreCue Alpha'
};
% {site , condition , leftVar , rightVar}
cases = {
  'UF', 'Choice',     'Pxx_choice_left_final',     'Pxx_choice_right_final'
  'UF', 'Instructed', 'Pxx_instructed_left_final', 'Pxx_instructed_right_final'
  'UCD','Choice',     'Pxx_choice_left_final',     'Pxx_choice_right_final'
  'UCD','Instructed', 'Pxx_instructed_left_final', 'Pxx_instructed_right_final'
};
caseLabels = {'UF_Choice','UF_Instructed','UCD_Choice','UCD_Instructed'};

nReps  = 100;   % CV repetitions (same as paper)
nFold  = 10;    % folds        (same as paper)

% SVM backend:  false -> MATLAB fitcsvm (no compiling);
%               true  -> LibSVM (exact paper pipeline; compile + addpath first).
USE_LIBSVM = false;

%% ---------------------------------------------------------------------
%  2b.  DOCUMENT CHANNEL INDICES AND GROUPS
%  ---------------------------------------------------------------------
siteNames = {'UF','UCD'};
channelTables = cell(numel(siteNames),1);
for s = 1:numel(siteNames)
    site = siteNames{s};
    labels = local_get_labels(site, UF_LOCS, UCD_CED, UCD_DROP);
    channelTables{s} = local_channel_table(site, labels, posteriorLabels, frontalLabels);
    writetable(channelTables{s}, fullfile(OUTDIR, sprintf('channel_groups_%s.csv', lower(site))));
end
channelTable = vertcat(channelTables{:});
writetable(channelTable, fullfile(OUTDIR,'channel_groups.csv'));
local_write_channel_markdown(channelTable, fullfile(OUTDIR,'channel_groups.md'));
fprintf('Channel index/group documentation saved to:\n  %s\n', fullfile(OUTDIR,'channel_groups.md'));

%% ---------------------------------------------------------------------
%  3.  MAIN LOOP  (data type x case)
%  ---------------------------------------------------------------------
summary = struct();
for d = 1 : size(dataTypes,1)
  dtName = dataTypes{d,1};  dtFolder = dataTypes{d,2};
  for c = 1 : size(cases,1)
    site = cases{c,1}; cond = cases{c,2}; leftVar = cases{c,3}; rightVar = cases{c,4};
    folder = fullfile(DATA_ROOT,'EEG',site,'EEG processed data','Filtered',dtFolder,cond);

    % --- channel labels in DATA order -> posterior & frontal indices ----
    labels   = local_get_labels(site, UF_LOCS, UCD_CED, UCD_DROP);
    [postIdx, frontIdx] = local_channel_indices(labels, posteriorLabels, frontalLabels);
    fprintf('\n[%s | %s %s] %d posterior, %d frontal of %d channels\n', dtName, site, cond, ...
        numel(postIdx), numel(frontIdx), numel(labels));
    fprintf('   posterior: %s\n   frontal:   %s\n', strjoin(labels(postIdx),','), strjoin(labels(frontIdx),','));

    Data_name = dir(fullfile(folder,'*.mat'));
    if isempty(Data_name)
        warning('No files in %s -- skipping.', folder); continue;
    end
    nSub = numel(Data_name);
    AccuracyAll = nan(nSub,nReps);  AccuracyPosterior = nan(nSub,nReps);  AccuracyFrontal = nan(nSub,nReps);
    nLeftRight  = nan(nSub,2);      chanUsed = labels(postIdx);  chanFront = labels(frontIdx);

    for m = 1 : nSub
        data = load(fullfile(folder, Data_name(m).name));
        L = data.(leftVar);  R = data.(rightVar);  nChan = size(L,1);
        assert(nChan == numel(labels), ...
            'Channel count mismatch for %s: data has %d rows but %s montage has %d labels.', ...
            Data_name(m).name, nChan, site, numel(labels));
        LeftData  = reshape(L, nChan, [])';   % nLeft  x nChan
        RightData = reshape(R, nChan, [])';   % nRight x nChan
        nLeftRight(m,:) = [size(LeftData,1) size(RightData,1)];

        Combined = cat(1, LeftData, RightData);
        Labels   = [-ones(size(LeftData,1),1); ones(size(RightData,1),1)];
        X_all   = zscore(Combined, 0, 2);
        X_post  = zscore(Combined(:, postIdx),  0, 2);
        X_front = zscore(Combined(:, frontIdx), 0, 2);

        for k = 1 : nReps
            Indices = crossvalind('Kfold', Labels, nFold);
            accA = nan(nFold,1);  accP = nan(nFold,1);  accF = nan(nFold,1);
            for i = 1 : nFold
                test = (Indices==i);  train = ~test;
                accA(i) = local_svm(X_all(train,:),   Labels(train), X_all(test,:),   Labels(test), USE_LIBSVM);
                accP(i) = local_svm(X_post(train,:),  Labels(train), X_post(test,:),  Labels(test), USE_LIBSVM);
                accF(i) = local_svm(X_front(train,:), Labels(train), X_front(test,:), Labels(test), USE_LIBSVM);
            end
            AccuracyAll(m,k)       = mean(accA);
            AccuracyPosterior(m,k) = mean(accP);
            AccuracyFrontal(m,k)   = mean(accF);
        end
        fprintf('  %-26s  L=%d R=%d   all=%.1f%%  post=%.1f%%  front=%.1f%%\n', Data_name(m).name, ...
            nLeftRight(m,1), nLeftRight(m,2), mean(AccuracyAll(m,:)), mean(AccuracyPosterior(m,:)), mean(AccuracyFrontal(m,:)));
    end

    tag = sprintf('%s__%s_%s', dtName, site, cond);   % e.g. CSD__UF_Choice
    clear S;
    S.dataType=dtName; S.site=site; S.cond=cond; S.caseLabel=sprintf('%s_%s',site,cond);
    S.chanAll=labels; S.chanPost=chanUsed; S.chanFront=chanFront;
    S.idxPost=postIdx; S.idxFront=frontIdx; S.nChanAll=numel(labels);
    S.chanUsed=chanUsed; % backward-compatible name for posterior ROI
    S.AccuracyAll=AccuracyAll; S.AccuracyPosterior=AccuracyPosterior; S.AccuracyFrontal=AccuracyFrontal;
    S.subjAccAll=mean(AccuracyAll,2);
    S.subjAccPosterior=mean(AccuracyPosterior,2);
    S.subjAcc=mean(AccuracyPosterior,2); % backward-compatible posterior name
    S.subjAccFront=mean(AccuracyFrontal,2);
    S.nLeftRight=nLeftRight;
    save(fullfile(OUTDIR,['Decoding_' tag '.mat']),'-struct','S');
    summary.(tag) = S;
  end
end

%% ---------------------------------------------------------------------
%  4.  GROUP STATS  (one-sample t-tests vs 50%, paired ROI contrasts)
%  ---------------------------------------------------------------------
fn = fieldnames(summary);
fprintf('\n=====================================  GROUP-LEVEL SUMMARY  =====================================\n');
fprintf('%-22s | %-16s | %-16s | %-16s | %-18s | %-18s | %s\n', ...
    'case (datatype)','ALL','POSTERIOR','FRONTAL','post-vs-all','post-vs-front','front-vs-all');
fprintf('%s\n', repmat('-',1,142));
for i = 1:numel(fn)
    aAll=summary.(fn{i}).subjAccAll; aPost=summary.(fn{i}).subjAcc; aFront=summary.(fn{i}).subjAccFront;
    [~,pAll ]=ttest(aAll ,50);
    [~,pPost]=ttest(aPost,50);
    [~,pFr  ]=ttest(aFront,50);
    [~,pPA]=ttest(aPost,aAll);    % paired posterior vs all
    [~,pPF]=ttest(aPost,aFront);   % paired posterior vs frontal
    [~,pFA]=ttest(aFront,aAll);    % paired frontal vs all
    nm=sprintf('%s_%s(%s)',summary.(fn{i}).site,summary.(fn{i}).cond,summary.(fn{i}).dataType);
    fprintf('%-22s | %5.2f%% p=%.2g | %5.2f%% p=%.2g | %5.2f%% p=%.2g | d=%+.2f p=%.3g | d=%+.2f p=%.3g | d=%+.2f p=%.3g\n', ...
        nm, mean(aAll),pAll, mean(aPost),pPost, mean(aFront),pFr, ...
        mean(aPost)-mean(aAll),pPA, mean(aPost)-mean(aFront),pPF, mean(aFront)-mean(aAll),pFA);
end
fprintf('%s\n', repmat('-',1,142));
fprintf('p under ALL/POSTERIOR/FRONTAL = one-sample vs 50%% chance.  d,p columns = paired t-tests.\n');

%% ---------------------------------------------------------------------
%  5.  BAR PLOT: one subplot per data type (all vs posterior, with stars)
%  ---------------------------------------------------------------------
YL = [45 68];
cols = [0.75 0.75 0.75; 0.20 0.45 0.70; 0.85 0.45 0.25];  % all / posterior / frontal
figure('Color','w','Position',[60 60 1280 480]);
for d = 1 : size(dataTypes,1)
    dtName = dataTypes{d,1};
    nC = numel(caseLabels);
    M=nan(nC,3); E=M; P=M; xlab=cell(1,nC);   % cols: 1=all 2=post 3=front
    for j = 1:nC
        tag = sprintf('%s__%s', dtName, caseLabels{j});
        xlab{j} = strrep(caseLabels{j},'_','\_');
        if ~isfield(summary,tag), continue; end
        sets = {summary.(tag).subjAccAll, summary.(tag).subjAcc, summary.(tag).subjAccFront};
        for s = 1:3
            a = sets{s};
            M(j,s)=mean(a); E(j,s)=std(a)/sqrt(numel(a)); [~,P(j,s)]=ttest(a,50);
        end
    end
    subplot(1,2,d); hold on;
    X=1:nC; b=bar(X,M,'grouped');
    for s=1:3, b(s).FaceColor=cols(s,:); end
    for s=1:3
        xc=b(s).XEndPoints;
        errorbar(xc,M(:,s),E(:,s),'k','linestyle','none','CapSize',6);
        for j=1:nC
            if ~isnan(M(j,s))
                text(xc(j),M(j,s)+E(j,s)+0.7, local_star(P(j,s)), ...
                    'HorizontalAlignment','center','FontSize',11,'FontWeight','bold');
            end
        end
    end
    set(gca,'XTick',X,'XTickLabel',xlab); ylabel('Decoding accuracy (%)'); ylim(YL);
    legend(b,{'All channels','Posterior ROI','Frontal ROI'},'Location','northwest');
    title(sprintf('%s pre-cue alpha decoding', dtName));
end
sgtitle('All vs posterior (parieto-occipital) vs frontal ROI   (* p<0.05  ** p<0.01  *** p<0.001;  one-sample vs 50%)');
saveas(gcf, fullfile(OUTDIR,'all_vs_posterior_vs_frontal.png'));

% tidy CSV (all cases x both data types, with p-values)
rows = {};
for i = 1:numel(fn)
    S=summary.(fn{i});
    [~,pa]=ttest(S.subjAccAll,50); [~,pp]=ttest(S.subjAcc,50); [~,pf]=ttest(S.subjAccFront,50);
    [~,ppa]=ttest(S.subjAcc,S.subjAccAll);
    [~,ppf]=ttest(S.subjAcc,S.subjAccFront);
    [~,pfa]=ttest(S.subjAccFront,S.subjAccAll);
    rows(end+1,:) = {S.dataType,S.site,S.cond,S.nChanAll,numel(S.chanUsed),numel(S.chanFront), ...
        mean(S.subjAccAll),pa,mean(S.subjAcc),pp,mean(S.subjAccFront),pf, ...
        mean(S.subjAcc)-mean(S.subjAccAll),ppa, ...
        mean(S.subjAcc)-mean(S.subjAccFront),ppf, ...
        mean(S.subjAccFront)-mean(S.subjAccAll),pfa}; %#ok<SAGROW>
end
T = cell2table(rows,'VariableNames',{'dataType','site','cond','nChanAll','nChanPost','nChanFront', ...
    'accAll_pct','pAll_vs50','accPost_pct','pPost_vs50','accFront_pct','pFront_vs50', ...
    'deltaPostMinusAll_pct','pPost_vsAll','deltaPostMinusFront_pct','pPost_vsFront', ...
    'deltaFrontMinusAll_pct','pFront_vsAll'});
writetable(T, fullfile(OUTDIR,'accuracy_summary.csv'));
fprintf('\nDone. Figure, accuracy_summary.csv, and channel_groups.* saved to:\n  %s\n', OUTDIR);

%% =====================================================================
%  LOCAL FUNCTIONS
%  =====================================================================
function s = local_star(p)
    if     p < 0.001, s = '***';
    elseif p < 0.01,  s = '**';
    elseif p < 0.05,  s = '*';
    else,             s = 'n.s.';
    end
end

function acc = local_svm(Xtr, ytr, Xte, yte, useLibSVM)
    if useLibSVM
        model = svmtrain(ytr, double(Xtr), '-t 0 -c 1 -q');   %#ok<SVMTRAIN> LibSVM
        [~,a,~] = svmpredict(yte, double(Xte), model, '-q');
        acc = a(1);
    else
        model = fitcsvm(double(Xtr), ytr, 'KernelFunction','linear', ...
            'BoxConstraint',1, 'Standardize',false);
        acc = 100 * mean(predict(model, double(Xte)) == yte);
    end
end

function labels = local_get_labels(site, ufLocs, ucdCed, ucdDrop)
% Channel labels in the SAME order as data rows.
% UF : 31 labels from .locs.   UCD: 58 labels from NeuroScan_58.ced.
    switch upper(site)
        case 'UF'
            labels = local_read_labels(ufLocs);
        case 'UCD'
            labels = local_read_labels(ucdCed);
            labels = labels(~ismember(upper(labels), upper(ucdDrop)));
        otherwise
            error('Unknown site %s', site);
    end
    labels = labels(:)';
end

function [postIdx, frontIdx] = local_channel_indices(labels, posteriorLabels, frontalLabels)
% Convert montage names into data-order indices for each reviewer-response ROI.
    labNorm  = local_norm_labels(labels);
    postIdx  = find(ismember(labNorm, posteriorLabels));
    frontIdx = find(ismember(labNorm, frontalLabels));
    assert(~isempty(postIdx), 'Posterior ROI matched zero channels.');
    assert(~isempty(frontIdx), 'Frontal ROI matched zero channels.');
    assert(isempty(intersect(postIdx,frontIdx)), 'Posterior and frontal ROIs overlap.');
end

function T = local_channel_table(site, labels, posteriorLabels, frontalLabels)
% Table documenting channel index, original label, normalized label, and group.
    labNorm = local_norm_labels(labels);
    postMask = ismember(labNorm, posteriorLabels);
    frontMask = ismember(labNorm, frontalLabels);
    group = repmat({'other'}, numel(labels), 1);
    group(postMask) = {'posterior'};
    group(frontMask) = {'frontal'};
    group(postMask & frontMask) = {'posterior+frontal'};
    T = table(repmat({upper(site)}, numel(labels), 1), (1:numel(labels))', labels(:), labNorm(:), group, ...
        'VariableNames', {'site','channel_index','channel_name','normalized_name','group'});
end

function labs = local_norm_labels(labels)
% Normalize only for matching; original names are preserved in output docs.
    labs = strtrim(regexprep(upper(labels), '''', ''));
end

function local_write_channel_markdown(T, fpath)
% Human-readable channel-index document to accompany channel_groups.csv.
    fid = fopen(fpath,'w');
    assert(fid > 0, 'Could not write channel markdown: %s', fpath);
    cleanup = onCleanup(@() fclose(fid));

    fprintf(fid, '# Channel groups for posterior-vs-frontal decoding\n\n');
    fprintf(fid, 'Generated by `PreCueAlpha_PosteriorDecoding.m` from the montage files in `channel_layouts/`.\n\n');
    fprintf(fid, 'Groups are assigned by channel label after upper-casing and stripping apostrophes. Channels not assigned to the posterior or frontal ROI are marked `other` and are included only in the all-channel decoder.\n\n');

    sites = unique(T.site, 'stable');
    for s = 1:numel(sites)
        site = sites{s};
        rows = find(strcmp(T.site, site));
        groups = T.group(rows);
        fprintf(fid, '## %s\n\n', site);
        fprintf(fid, '- Total channels: %d\n', numel(rows));
        fprintf(fid, '- Posterior channels: %d\n', sum(strcmp(groups,'posterior')));
        fprintf(fid, '- Frontal channels: %d\n', sum(strcmp(groups,'frontal')));
        fprintf(fid, '- Other channels: %d\n\n', sum(strcmp(groups,'other')));
        fprintf(fid, '| index | channel | group |\n');
        fprintf(fid, '|---:|---|---|\n');
        for r = rows'
            fprintf(fid, '| %d | %s | %s |\n', T.channel_index(r), T.channel_name{r}, T.group{r});
        end
        fprintf(fid, '\n');
    end
end

function labs = local_read_labels(fpath)
% Minimal reader for EEGLAB .locs / .ced text files.
    assert(exist(fpath,'file')==2, 'Channel file not found: %s', fpath);
    [~,~,ext] = fileparts(fpath);
    fid = fopen(fpath,'r'); raw = textscan(fid,'%s','Delimiter','\n','Whitespace',''); fclose(fid);
    lines = raw{1}; lines(cellfun(@isempty,strtrim(lines))) = [];
    labs = {};
    if strcmpi(ext,'.ced')
        hdr = regexp(lines{1}, '\t', 'split');
        col = find(strcmpi(strtrim(hdr),'labels'), 1);
        if isempty(col), col = 2; end
        for i = 2:numel(lines)
            tok = regexp(lines{i}, '\t', 'split');
            if numel(tok) >= col, labs{end+1} = strtrim(tok{col}); end %#ok<AGROW>
        end
    else
        for i = 1:numel(lines)
            tok = regexp(strtrim(lines{i}), '\s+', 'split');
            labs{end+1} = strtrim(tok{end}); %#ok<AGROW>
        end
    end
    labs = labs(:);
end
