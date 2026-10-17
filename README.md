# BEBO EEG pipeline (rerun 2026-10)

Analysis code for the BEBO encoding EEG study (younger vs older adults):
gamma item specificity (spectral RSA), the alpha/beta subsequent memory effect
(SME), and trial-level models of memory detail.

**The binding analysis plan is `BEBO/analysis_plan_2026-10-17_FINAL.txt`.**
Every setting lives in `bebo_config.m`, so scripts and plan cannot drift
apart. The pipeline reads only the preprocessed data; no earlier output is
used anywhere.

## Folder layout on the grid

```
~/BEBO2/bebo-pipeline/     this repository
~/BEBO2/PreprocData/       <ID>_Enc_allcond.mat (input)
~/BEBO2/Elec/elec.mat      electrode positions
~/BEBO2/Performance_All.txt
~/BEBO2/Results_v2/        all output (created by the pipeline)
```

## Running

On the grid, from inside `~/BEBO2/bebo-pipeline`:

| Command | What it submits |
|---|---|
| `bash submit_pipeline.sh one` | TEST: S1 and S2 for the first person only |
| `bash submit_pipeline.sh` | everything (S1, S2, G1-G4, P1) |
| `bash submit_pipeline.sh group` | G1-G4 only (S1 and S2 done) |
| `bash submit_pipeline.sh g3` | G3, G4, P1 only |
| `bash submit_pipeline.sh figs` | P1 (figures) only |
| `bash submit_pipeline.sh behavior` | B1 (behaviour) only |
| `bash submit_pipeline.sh count` | C1: trial counts of all 96 (93 + 3 excluded) from PreprocData |

Each job starts only when the jobs it needs have finished **without error**.
Logs are in `~/logs/`. Then, locally, on the downloaded `Results_v2/tables`:
`Rscript R/models_results.R`, then `Rscript R/models_plots.R`.

## Steps

| Step | Script | Does | Output (in Results_v2/) |
|---|---|---|---|
| S1 | `s1_gamma_subject.m` | per person: DFT filter 50/100 Hz; DPSS TFR 35-100 Hz; aperiodic removal; Isim / Csim pair similarity; reference matched on lag and session position; t-map | `gamma_subject/<ID>_gamma.mat`, `trial_log/` |
| S2 | `s2_alpha_subject.m` | per person: Morlet 2-30 Hz log power per trial; later 'old' / 'similar' averages | `alpha_subject/<ID>_alpha.mat`, `trial_log/` |
| G1 | `g1_ispec_group.m` | Ispec cluster test (t-maps vs 0); per-person Isim, Csim, Ispec, drift check; group-average matrices | `group/ispec_cluster.mat`, `group/ispec_group_means.mat`, `tables/person_ispec.csv` |
| G2 | `g2_sme_group.m` | SME cluster test (old vs similar); per-person SME | `group/sme_cluster.mat`, `tables/person_sme.csv` |
| G3 | `g3_extract.m` | single-trial tables (Ispec cluster only); split half | `tables/items.csv`, `pairs.csv`, `person_split.csv` |
| G4 | `g4_person_stats.m` | age differences, drift check, correlation families, split half | `stats/person_stats.txt`, `tables/person_all.csv` |
| P1 | `p1_figures.m` | group-level figures incl. Fig3 similarity matrices (no computation) | `figures/` |
| C1 | `c1_count_trials.m` | trial counts per person (no EEG computation); sample check | `trial_log/trial_counts_all96.csv`, `trial_counts_summary.txt` |
| B1 | `b1_behavior.m` | behaviour figure + stats (old behavior_analysis.m, unchanged) | `figures/Fig_Behavior.*`, `stats/Behavior_stats.txt` |
| R | `R/models_results.R` | condition models, pair model, memory model | `stats/models_stats.txt`, `stats/models.rds` |
| R | `R/models_plots.R` | figures of the models (no computation) | `figures/` |

Shared: `bebo_config.m` (all settings), `bebo_ids.m` (sample),
`bebo_trials.m` (trial rule), `bebo_cluster_test.m` (cluster permutation
test), `bebo_init_fieldtrip.m`, `spectral_rsa.m`, `subtract_mean_noise_spectrum.m`.

## Testing

`test/` holds a simulator with known planted effects (item code, category
code, slow drift, line noise, alpha SME) and `test/run_test_pipeline.m`,
which runs every MATLAB step on simulated data. `test/spm_bwlabel.m` is a
stand-in for local testing only (the local FieldTrip has no SPM); the grid
FieldTrip has the real one.
