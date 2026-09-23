# Reviewer response — work log

Manuscript: *Neural Mechanisms of Willed Attention Control* (eLife)
Referee report received: 2026-02-26 · Log started: 2026-06-06

Status key: ✅ done · 🔄 in progress · ⬜ not started · 💬 discussion/text-only (no new analysis)

---

## Currently working on

### R1 · Comment 2 — Decode using posterior (parieto-occipital) and frontal-control electrodes 🔄

**Reviewer:** EEG decoding used the *entire* electrode topography rather than a biologically motivated posterior ROI. Alpha-mediated spatial attention is traditionally localized to parieto-occipital sensors; the full set risks contamination by micro-saccades / muscle artifact.

**Our plan (from referee notes):** "We can try to use posterior electrodes to do decoding and see what happens."

**What was built (2026-06-06):**

- `PreCueAlpha_PosteriorDecoding.m` — a drop-in variant of `code/EEG/CSD_PreCueAlphaPowerDecoding.m`. Identical pipeline (z-score across channels → linear LibSVM or `fitcsvm`, c=1 → 100×10-fold CV) except features are evaluated as all channels, posterior ROI, and frontal control ROI on the same folds.
- Posterior ROI = P/PO/O/CP/TP-line electrodes. Resolved **per montage** by label:
  - **UF (31 ch, BrainAmp / `brainproduct31.locs`):** 15 posterior channels → P3, P4, O1, O2, P7, P8, Pz, Oz, CP1, CP2, CP5, CP6, TP9, TP10, POz. *(verified against the .locs channel order)*
  - **UCD (58 ch, Neuroscan `NeuroScan_58.ced`):** 26 posterior channels. The montage uses old Neuroscan names (apostrophes); the script strips apostrophes and matches by label. Selected: P3', P4', O1', O2', T5', T6', Pz', PzA, C3P, C4P, Oz', C1P, C2P, TCP1, TCP2, P1', P2', P5', P6', P1P, P2P, P3P, P4P, CB1, CB2, PzP. (T5'/T6' = P7/P8; CB1/CB2 = occipital; P*P/PzP = PO line; C*P/TCP = CP/TP line.) Everything excluded is frontal/central/anterior-temporal — i.e. the artifact-prone sensors the reviewer flagged.
- Frontal ROI = Fp/F/FC/FT-line electrodes plus anterior Neuroscan-specific labels. UF has 11 frontal channels; UCD has 21 frontal channels.
- Runs all four cells: UF-Choice, UF-Instructed, UCD-Choice, UCD-Instructed.
- Outputs per-subject accuracy, group t-tests vs 50%, paired ROI contrasts, and an all-vs-posterior-vs-frontal bar plot.

**Update (2026-06-29):**

- Moved the response folder from `C:\Users\yujunchen\Documents\Changhao\2022-2025\WilledAttention\Reviewer_Response` to `N:\Experimental_Data\yujunchen\projects\WilledAttention\Reviewer_Response`.
- Patched paths so analysis reads EEG data from Changhao's shared data root: `N:\Experimental_Data\Changhao Xiong\2022-2025\WilledAttention`.
- Bundled montage files under `Reviewer_Response/channel_layouts/` and added channel documentation outputs: `results_posterior/channel_groups.md`, `channel_groups.csv`, `channel_groups_uf.csv`, and `channel_groups_ucd.csv`.
- Added channel-count assertions so the script fails fast if a data matrix does not match the expected UF/UCD montage.
- MATLAB R2024a `checkcode` passed with no script issues.

**Open items before this is ✅:**

- [x] Run the script in MATLAB — UF cells run via built-in `fitcsvm` (USE_LIBSVM toggle added; LibSVM optional). Record final numbers below.
- [x] Confirm UCD posterior channel list — resolved from `NeuroScan_58.ced` (26 channels).
- [x] Add frontal control ROI and channel index/name/group documentation.
- [ ] Re-run the updated all/posterior/frontal script after the 2026-06-29 patch so `accuracy_summary.csv` and `all_vs_posterior_vs_frontal.png` reflect the frontal-control comparison.
- [ ] Decide expected/interpretable outcome: if posterior-only accuracy ≈ whole-head, the whole-head result is **not** artifact-driven (answers the reviewer); if it drops, discuss why.
- [ ] Add resulting numbers + figure to the manuscript / response letter.

**Comparison built in:** the script decodes **all channels**, the **posterior ROI**, and a **frontal ROI** on the *same CV folds* (paired), and runs **both CSD and voltage** data types. Outputs: `Decoding_<datatype>__<case>.mat`, `all_vs_posterior_vs_frontal.png` (2 subplots), `accuracy_summary.csv`, and `channel_groups.*`. Group section prints all/posterior/frontal tests vs 50% plus paired posterior-vs-all, posterior-vs-frontal, and frontal-vs-all contrasts; significance stars on the figure are one-sample tests vs 50%.

**KEY FINDING — CSD vs voltage (why UCD looked high).** The paper's Figure 4A reports **scalp-voltage** decoding. Verified with sklearn (= LibSVM):

| | CSD | Voltage | Paper 4A |
|---|---|---|---|
| UF Choice | 55.8% | 55.5% | 55.74% |
| UCD Choice | 63.2% | 57.4% | 57.29% |

So the earlier ~63% UCD number was because the script read the **CSD** folder; **voltage matches the paper exactly**. CSD (surface Laplacian) sharpens topography and adds ~6% on the dense 58-ch UCD montage; on the 31-ch UF montage CSD≈voltage. It is **not** a solver issue (sklearn-LibSVM also gives 63% on CSD) and only ~2% is class imbalance (balanced UCD-CSD = 61%). The script now reports both so the figure shows the voltage panel (paper-matching) and the CSD panel (shows the sharpening effect).

**Results (fill in after running):**

| Dataset | Condition | All-channel acc | Posterior acc | Frontal acc | All p vs 50% | Posterior p vs 50% | Frontal p vs 50% | Δ (post−all), paired p | Δ (post−front), paired p |
|---|---|---|---|---|---|---|---|---|---|
| UF | Choice | | | | | | | | |
| UF | Instructed | | | | | | | | |
| UCD | Choice | | | | | | | | |
| UCD | Instructed | | | | | | | | |

> Note: Populate the table from the updated script's printed GROUP-LEVEL SUMMARY / `accuracy_summary.csv` after re-running the 2026-06-29 version. Expected pattern: Choice cells above 50% for all/posterior, frontal weaker if the effect is posterior alpha rather than frontal artifact; Instructed cells near 50%.

---

## All referee comments — tracker

### Reviewer #1 (public review)
- **C1 — fMRI beta-series lacks temporal precision / spontaneous vs pre-planned.** ⬜ 💬 — *Plan:* discuss with Jesse & Ron; cite Abhijit's theta-timing paper for decision timing.
- **C2 — Posterior-electrode decoding.** 🔄 — *see above.*
- **C3 — Physiological cost: same behavior despite extra frontoparietal recruitment (arousal/effort).** ⬜ — *Plan:* compare post-cue alpha following choice vs instructed cues; discuss arousal/compensation.

**R1 recommendations:** ⬜ add color bar + threshold to Fig 2; ⬜ identical scaling for Fig 5B (UF vs UCD); ⬜ permutation/sensitivity test for neural-efficiency ratio; ⬜ report univariate posterior alpha asymmetry; ⬜ state which TRs used for fMRI decoding + justify no time-resolved; ⬜ discuss behavioral parity; ⬜ clarify site-effect modeling.

### Reviewer #2 (public review)
- **C1 — Is "willed attention" really *will*? ("mental coin flip").** ⬜ 💬 — *Plan:* discuss with Jesse & Ron; discuss other phenomena under willed attention + future directions.
- **C2 — Decision process + WHEN the decision was made.**
  - **2a — what processes the decision engages.** ⬜ 💬 — *Plan:* talk to Jesse & Ron; reuse answer to R1-C1.
  - **2b — pre-cue prediction: causal bias vs subject decided early.** ⬜ 💬 — *Plan:* bring in the pre-cue rhythmic-scanning idea; temper causal claims; cite Abhijit's theta work for timing.
- **C3 — Choice bias L vs R → 50% may not be true chance; do permutation.** ⬜ — *Plan:* report %left vs %right per subject; decode with balanced data; build permutation (label-shuffle ≥1000×) baseline. **(Related to the posterior code — can reuse the same pipeline.)**
- **C4 — Novelty beyond Bengson et al.** ⬜ 💬 — *Plan:* explicitly state contributions (decoding, EEG↔fMRI link, two-site).

**R2 recommendations:** mostly map to C1–C4; **(4)** = Monte-Carlo/permutation MVPA (≥1000 shuffles per subject) — shares work with R2-C3 / R3-C1.

### Reviewer #3 (public review)
- **C1 — How "above chance" is determined; multiple comparisons; FWE for time-resolved EEG; effect sizes + CIs.** ⬜ — *Plan:* check FieldTrip for cluster-based permutation (Maris & Oostenveld 2007); report CIs (Combrisson & Jerbi 2015).
- **C2 — Cross-validation: folds blocked by run vs random; trial counts after rejection / short-ITI removal; %left vs %right; balanced folds.** ⬜ — *Plan:* report sampling into 10 folds + trial counts. **(Overlaps R2-C3.)**
- **C3 — ROI definition / circularity (Kriegeskorte 2009); spheres vs anatomical; radius/voxel count.** ⬜ — *Plan:* table of voxel counts per ROI; redo decoding in standard atlas-defined ROIs.
- **C4 — Neural-efficiency ratio unstable.** ⬜ 💬 — *Plan:* discuss; describe the bootstrap (10,000-iteration) approach already used for the distribution.

**R3 minors:** ⬜ report sex imbalance (UCD) as limitation; ⬜ clarify site-difference moderators (TR, montage, simultaneous vs separate); ⬜ soften causal language → "consistent with".

---

## Shared analyses (do once, answers several comments)

1. **Permutation / Monte-Carlo MVPA baseline** (≥1000 label shuffles per subject) → R2-C3, R2-rec4, R3-C1. *Can be added to `PreCueAlpha_PosteriorDecoding.m` as a label-shuffle loop.*
2. **Choice-bias + balanced decoding** (%L vs %R per subject; balanced folds) → R2-C3, R3-C2.
3. **Trial counts & fold structure reporting** → R3-C2.
4. **Atlas-defined ROI decoding robustness check** → R3-C3.
5. **Neural-efficiency sensitivity/bootstrap** → R1-rec, R2-min2, R3-C4.

## People to consult
- **Jesse Bengson / Ron Mangun:** R1-C1, R2-C1, R2-C2a, R2-C4 (conceptual/framing).
- **Abhijit (theta timing paper):** R1-C1, R2-C2b (decision timing).

---
*Update this log as each item moves ⬜ → 🔄 → ✅.*
