# WilledAttention

## Purpose

This repository contains the analysis code written in response to peer review for the
manuscript *Neural Mechanisms of Willed Attention Control* (Xiong, Yang, Kim, Meyyappan,
Bengson, Mangun & Ding; submitted to *eLife*). The study contrasts **willed (free-choice)**
attention against **instructed** attention using combined EEG and fMRI collected at two sites
(University of Florida and UC Davis), and asks:

1. What extra brain machinery supports self-generated attention decisions?
2. Can the direction of attention be decoded from brain activity?
3. Does pre-cue brain state bias the upcoming choice?

The scripts here answer four specific reviewer concerns raised during revision by re-running
or extending pieces of the original decoding analysis:

- **R1 — Posterior electrode decoding** (`01_R1_Posterior_Electrode_Decoding/`): re-runs EEG
  choice-direction decoding restricted to posterior electrodes (vs. whole-scalp) to check that
  the original decoding result isn't driven by frontal/non-posterior artifacts.
- **R1/R2 — Post-cue alpha and arousal** (`02_R1_R2_PostCue_Alpha_Arousal/`): compares post-cue
  occipital alpha power for choice vs. instructed trials, to address whether willed attention
  might simply reflect greater effort/arousal.
- **R2/R3 — Choice-bias permutation testing** (`03_R2_R3_Choice_Bias_Permutation/`): tests
  whether left/right choice imbalance could artificially inflate above-chance decoding, using
  balanced decoding and subject-level permutation baselines.
- **R3 — Decoding inference** (`04_R3_Decoding_Inference/`): EEG (time-resolved, cluster-based)
  and fMRI (ROI-based) statistical inference for the choice-direction decoding results,
  including cluster/permutation correction.

This repository is **not** the full experimental codebase. It contains only the scripts
written for the reviewer-response package, and it contains **no raw or processed data,
figures, or results** — see Dependencies below for what you need to supply yourself.

## Contents

```
README.md
Reviewer_Response/
  README.md                                       folder guide
  Paper_Summary.md                                 plain-language summary of the paper
  Reviewer_Response_Log.md                         dated work log tracking each reviewer comment
  01_R1_Posterior_Electrode_Decoding/
    channel_layouts/                                EEG cap montage definitions (.ced/.locs)
    code/
      PreCueAlpha_PosteriorDecoding.m                runs the posterior/frontal/all-channel decoding
      PlotVoltageOnly_FromSavedResults.m             plots the saved decoding results
    results/channel_groups.md                        electrode-group documentation
  02_R1_R2_PostCue_Alpha_Arousal/
    code/reviewer2_postcue_alpha_choice_vs_instructed.m
  03_R2_R3_Choice_Bias_Permutation/
    code/reviewer3_choice_bias_balanced_permutation.m
  04_R3_Decoding_Inference/
    code/EEG/reviewer4_eeg_cluster_inference.m
    code/fMRI/reviewer4_fmri_roi_inference.m
```

Each numbered folder is a self-contained, complete analysis for one reviewer comment; only
`01_R1_Posterior_Electrode_Decoding/` has an internal two-step order (decode, then plot — see
below). An earlier draft of this repository also included a top-level `code/` folder holding
older, iterative working versions of the same scripts (plus Word-document-generation helper
scripts); those were superseded by the versions kept here and have been removed as
scratch/duplicate content.

## How to Use

None of these scripts include or download data — the underlying EEG/fMRI recordings and
preprocessed derivatives are not part of this repository and must be supplied separately by
placing them in the directory layout each script expects (below). All scripts are standalone
MATLAB `.m` files: open in MATLAB and run, or run non-interactively with
`matlab -batch "run('script_name.m')"` from the script's folder.

### 01 — Posterior electrode decoding (2 steps, run in order)

1. **`PreCueAlpha_PosteriorDecoding.m`**
   - Input: filtered, epoched EEG `.mat` files expected at
     `<DATA_ROOT>/EEG/<site>/EEG processed data/Filtered/<condition folder>/<Choice|Instructed>`
     for both sites (`UF`, `UCD`). Edit the `DATA_ROOT` variable near the top of the script
     (currently set to the original author's path) to point at your local copy of this data.
   - Does: pre-cue (-500 to 0 ms) alpha-power decoding of attend-left vs. attend-right, on
     all channels vs. a posterior ROI vs. a frontal control ROI, for both CSD and voltage
     features, using the montage files in `../channel_layouts/`.
   - Output: `Decoding_<tag>.mat` result files, `channel_groups.csv`/`.md`, and
     `accuracy_summary.csv`, written to a `results_posterior/` folder created next to the
     script.
2. **`PlotVoltageOnly_FromSavedResults.m`**
   - Input: the `Decoding_Voltage__*.mat` files produced by step 1, read from that same
     `results_posterior/` folder — run step 1 first.
   - Does: reproduces the reviewer-facing bar plot (all/posterior/frontal decoding accuracy)
     from the saved voltage-decoding results only.
   - Output: a figure plus the underlying summary table, written to `results_posterior/`.

### 02 — Post-cue alpha, choice vs. instructed

- **`reviewer2_postcue_alpha_choice_vs_instructed.m`**
  - Input: epoched, ICA-cleaned EEG `.mat` files expected at
    `<repo root>/EEG/<UF|UCD>/EEG processed data/Filtered/Epoched ICA/<Choice|Instructed>`,
    where `<repo root>` is computed automatically from the script's own location (three levels
    up from `code/`) — place the `EEG/` data folder there, or edit `cfg.rootDir` to point
    elsewhere.
  - Does: computes induced, baseline-normalized post-cue alpha power in sliding time windows
    for choice vs. instructed trials at both sites, at Fz/Oz and averaged across channels.
  - Output: figures and summary statistics under `Result Figures/Reviewer2_PostCueAlpha/`.

### 03 — Choice-bias permutation testing

- **`reviewer3_choice_bias_balanced_permutation.m`**
  - Input: pre-cue voltage alpha-power `.mat` files expected at
    `<repo root>/EEG/<UF|UCD>/EEG processed data/Filtered/Voltage PreCue Alpha/Choice`,
    containing `Pxx_choice_left_final` / `Pxx_choice_right_final` variables (same convention as
    script 02). `<repo root>` is again derived automatically from the script location.
  - Does: quantifies left/right choice proportions per subject, reruns choice decoding after
    balancing trial counts, and builds an empirical chance distribution via within-subject
    label permutation (`cfg.nPerm`, default 10000 — set the environment variable
    `REVIEWER3_FAST=1` before launching MATLAB for a fast smoke test with fewer permutations).
  - Output: figures and tables under `Result Figures/Reviewer3_ChoiceBias/`.

### 04 — Decoding inference (EEG and fMRI, independent of each other)

- **`code/EEG/reviewer4_eeg_cluster_inference.m`**
  - Input: subject-by-time decoding-accuracy curves saved by the *original* (not included in
    this repository) time-resolved decoding scripts. Point `cfg.rootDir`/the script's expected
    input path at wherever those saved curves live.
  - Does: cluster-based permutation testing over time (FieldTrip if available, otherwise a
    built-in sign-flip max-cluster fallback), plus effect sizes and confidence intervals at
    every time point.
  - Output: corrected cluster figures and statistics under
    `Result Figures/Reviewer4_DecodingInference/EEG/`. Set `REVIEWER4_FAST=1` for a quick,
    low-permutation-count run.
- **`code/fMRI/reviewer4_fmri_roi_inference.m`**
  - Input: subject-by-CV-repeat-by-ROI decoding accuracies saved by the *original* (not
    included in this repository) fMRI ROI decoding scripts.
  - Does: averages CV repeats within subject, tests subject-level accuracy against chance per
    ROI, reports confidence intervals and within-subject effect sizes, and corrects across
    ROIs.
  - Output: figures and statistics under `Result Figures/Reviewer4_DecodingInference/fMRI/`.
    Set `REVIEWER4_FAST=1` for a quick run with fewer permutations/bootstraps.

## Dependencies

- MATLAB (all six scripts are plain MATLAB). Statistics and Machine Learning Toolbox is
  required for SVM decoding (`fitcsvm`) and t-tests; LibSVM is an optional alternative for the
  decoding scripts.
- FieldTrip is used opportunistically by `reviewer4_eeg_cluster_inference.m` for cluster-based
  permutation statistics; if it isn't on the MATLAB path, the script falls back to a built-in
  sign-flip max-cluster implementation, so FieldTrip is optional but recommended for exact
  parity with the manuscript's reported statistics.
- No external data is bundled. Each script above documents the directory layout it expects —
  you must supply the corresponding EEG/fMRI derivatives yourself.
