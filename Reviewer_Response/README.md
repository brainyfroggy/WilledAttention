# Reviewer_Response folder guide

This folder consolidates the reviewer-response analyses for *Neural Mechanisms of Willed
Attention Control*. See the repository root `README.md` for full Purpose/How-to-Use/
Dependencies details; this file is a short structure guide.

## Structure

- `01_R1_Posterior_Electrode_Decoding/`
  - Response to the concern that whole-scalp EEG decoding might reflect non-posterior
    artifacts. Contains posterior/frontal ROI decoding code and channel-layout documentation.
- `02_R1_R2_PostCue_Alpha_Arousal/`
  - Response to the concern that willed attention may involve extra effort/arousal despite
    behavioral parity. Contains the post-cue alpha choice-vs-instructed analysis script.
- `03_R2_R3_Choice_Bias_Permutation/`
  - Response to concerns that left/right choice imbalance could inflate decoding above 50%.
    Contains the balanced-decoding and permutation-baseline script.
- `04_R3_Decoding_Inference/`
  - Response to concerns about fMRI ROI inference and EEG time-resolved multiple comparisons.
    Contains cluster/permutation-based inference scripts for EEG and fMRI.

Each numbered folder is a complete, standalone analysis for one reviewer comment (only
folder 01 has an internal two-step run order — see the root README's How to Use section).

## What was excluded

Per the archival policy for this repository, referee-report/manuscript source documents, the
generated Word response draft, all data and result outputs (figures/tables/`.mat`/`.csv`
under `results/` or `Result Figures/`), and QA-render scratch output were left out — this
repository is source code only. An earlier version of this folder also kept a duplicate,
pre-cleanup `code/` working-scripts folder and several Word-document-generation helper
scripts retained "for continuity with earlier work"; those were superseded by the scripts in
the numbered folders above and have since been removed as scratch/duplicate content.
