# Hybrid L1 tracking: where low-pT barrel tracks are lost

Truth tracing of the hybrid tracklet chain for the low-pT (v2) configuration, to
find which step loses barrel tracks at 0.5–1 GeV (efficiency 0.02–0.4, while
forward tracks at the same pT are found at 0.9+ and offline at 0.95+).

## How it works

The emulation is run on the STARlight QED μμ step-2 files
(`Step2_DIGI_RAW_CMSSW_14_0_6_v4_tpminpt0p3`) with Nicole's v2 code from this
repository (commits of 2026-09-24/25) plus a diagnostic patch that writes, for
true particles, one line per step to `truthtrace.txt`:

| record | written in | content |
|---|---|---|
| `E`, `T`, `S` | TrackletEventProcessor | event; truth particles; stubs with their TP |
| `C` | TrackletCalculatorBase (barrel seeds) | a true seed pair reached the tracklet calculator |
| `F` | same | the calculator rejected it (`rz`: rinv/z0 range, `sec`: sector) |
| `P` | same | projection to each layer: status (0 ok, 1 z / 2 φ out of range), exact, Taylor and fixed-point φ and z |
| `A` | same | projection routed to the match memories |
| `M` | MatchProcessor (barrel) | every match candidate: stub TP, r-φ and z residual, window, accepted |
| `B` | FitTrack | track builder: matched layers (≥ 2 needed) |

Files:

- `truthtrace.patch`: against this repo's `TrackletCalculatorBase.cc`, `MatchEngineUnit.cc` and
  `MatchProcessor.cc`, `L1FPGATrackProducer.cc` (truth hand-over cut lowered from
  1.0 to 0.2 GeV, diagnostic only) and `L1TrackHitNtupleMaker.cc` (restores the
  `cluster_*` declarations; the committed file does not compile without them).
- `release_FitTrack.patch`, `release_TrackletEventProcessor.patch`: against the
  CMSSW_15_1_0_patch3 release files (not in this repo).
- `scripts/diag_cfg.py`: wraps `test/rerunL1_trackhitntuple_cfg.py`. Besides input
  and output, it **enables the L1 track finding** (`L1TrackTrigger_step` is
  commented out of the committed schedule, which then re-makes the dummy stubs but
  neither re-runs the track finding nor the stub/track truth association), runs
  only the prompt hybrid tracking, and sets `readMoreMcTruth`.
- `scripts/run_one.sh`: one step-2 file per job.
- `scripts/analyze_trace.py`, `scripts/stage_report.py`: the analysis;
  `scripts/compare_variants.py`: efficiency, fakes and duplicates per cut variant.

Dummy stubs carry the TP of their inner cluster only, so the trace takes a stub's
truth from the inner cluster (the emulation's own `Tracklet::tpseed()` requires
both clusters and is always 0 for dummy stubs).

The rerun reproduces the v2 ntuples (hybrid efficiency, barrel 0.5–0.75 / 0.75–1 /
1–1.5 GeV: 0.03 / 0.45 / 0.75 vs 0.02 / 0.40 / 0.76 in v2; forward 0.5–0.75 GeV:
0.91 vs 0.93).

## Result (8 files, 2556 events, barrel |η| < 0.8)

| pT [GeV] | muons | seed pair formed | tracklet passes calculator | ≥ 1 match | ≥ 2 matched layers | found |
|---|---|---|---|---|---|---|
| 0.5–0.75 | 346 | 0.93 | 0.90 | 0.24 | 0.03 | 0.03 |
| 0.75–1 | 108 | 1.00 | 1.00 | 0.82 | 0.45 | 0.45 |
| 1–1.5 | 60 | 1.00 | 1.00 | 0.97 | 0.78 | 0.75 |

Seeding and the tracklet calculator work; tracks are lost in **matching**,
in two steps that both use fixed φ windows:

1. **MatchEngine** (candidate formation): the true stub becomes a candidate for
   49% / 40% / 17% of L1L2 → L3 / L4 / L5 projections at 0.5–0.75 GeV
   (92% / 79% / 84% at 1–2 GeV). With the bend check bypassed for dummy stubs,
   the selection is `|Δ fine-φ| < 3` (`idphicut = 3` in `MatchEngineUnit.cc`)
   plus the rz window.
2. **MatchProcessor**: z passes 90–100%, but the r-φ window rejects the true
   stub. L1L2 → L3 (window 0.10 cm): 26% pass at 0.5–0.75 GeV, 73% at 0.75–1,
   93% at 1–2; the median residual is 2–4 windows at 0.5–0.75 GeV and 0.4–0.6
   at 1–2 GeV. The windows are per seed and layer, independent of pT.

Secondary: for outer seeds (L5L6) and outer projections, the Taylor /
fixed-point projection deviates from the exact helix at low pT (L5L6 → L1 at
0.5–0.75 GeV: ~5.7 cm vs ~0.5 cm); the reach-layer gate is commented out in
this version.

## Relaxing the two cuts (constant, not pT dependent)

`TT_IDPHICUT` (MatchEngine fine-φ cut, default 3), `TT_REACHGATE` (reach check before barrel projections, default off) and `TT_PHIWIN_SCALE` (factor on
the MatchProcessor r-φ window, default 1) are read from the environment by the
diagnostic build. Hybrid efficiency (`tp_nmatch > 0`), same 2556 QED μμ events:

| variant | barrel 0.5–0.75 | barrel 0.75–1 | barrel 1–2 | forward 0.5–0.75 | stubs/track, barrel 0.5–1 | duplicates |
|---|---|---|---|---|---|---|
| baseline (3, ×1) | 0.03 | 0.45 | 0.80 | 0.91 | 4.09 | 7.1% |
| MatchEngine 6 | 0.06 | 0.57 | 0.82 | 0.92 | 4.11 | 7.4% |
| window ×4 | 0.18 | 0.74 | 0.94 | 0.92 | 4.25 | 4.8% |
| MatchEngine 6, window ×4 | 0.26 | 0.83 | 0.96 | 0.93 | 4.42 | 5.3% |
| MatchEngine 8, window ×4 | 0.27 | 0.82 | 0.96 | 0.93 | 4.41 | 5.2% |

- The MatchEngine cut saturates at 6: beyond it the virtual-module routing, not
  the cut, decides which stubs are seen.
- The window cannot go beyond about ×4: with ×6, ×8 and ×16 accepted residuals
  overflow the 12-bit r-φ residual word (`phiresidbits_ = 12`;
  `Residual.cc: assert(!fpgaphiresid_.atExtreme())` aborts the job).
- No fakes appear in UPC μμ (nothing to make them from). EPOS pPb, 600 events
  (1400 stubs/event), baseline → MatchEngine 6 + window ×4:

  | | barrel 0.5–0.75 | barrel 0.75–1 | barrel 1–2 | forward 0.5–0.75 | tracks/event | not loosely genuine | duplicates |
  |---|---|---|---|---|---|---|---|
  | baseline | 0.02 | 0.27 | 0.68 | 0.63 | 27.5 | 0.4% | 5.0% |
  | ME 6, ×4 | 0.19 | 0.59 | 0.81 | 0.66 | 33.2 | 0.7% | 4.3% |

  Barrel fakes rise from 0.1–0.6% to 0.7–1.8% per pT bin.

QED e⁺e⁻, 40 step-2 files (18 751 events), same events for all columns, all
|η| < 2.4 (hybrid found = `tp_nmatch > 0`; GNN from `gnn_tracking`, trained on the
v2 e⁺e⁻ ntuples with these 18 751 events excluded from training and validation;
statistical uncertainties ≤ 0.01):

| TP pT [GeV] | TPs | hybrid, v2 cuts | hybrid, ME 6 + ×4 | GNN + 1 crossing/layer |
|---|---|---|---|---|
| 0.3–0.5 | 4 376 | 0.02 | 0.02 | 0.31 |
| 0.5–0.75 | 5 824 | 0.15 | 0.24 | 0.67 |
| 0.75–1 | 6 410 | 0.30 | 0.49 | 0.83 |
| 1–2 | 10 780 | 0.47 | 0.63 | 0.88 |
| > 2 | 1 964 | 0.65 | 0.77 | 0.88 |
| fake rate | | 0.22% | 0.26% | 1.8% |

Barrel / forward efficiency of the hybrid, v2 cuts → ME 6 + ×4: barrel 0.75–1 GeV
0.17 → 0.47, 1–2 GeV 0.43 → 0.66; forward (|η| > 1.8) 0.75–1 GeV 0.52 → 0.56,
1–2 GeV 0.60 → 0.64. The relaxed cuts help electrons in the barrel but hardly in
the forward region, where the loss is not the match window.

What remains at 0.5–0.75 GeV (barrel, ME 6 + ×4): the true stub is not a
candidate for 27% / 47% / 77% of L1L2 → L3 / L4 / L5 projections (routing), the
window still rejects 30–50% of the rest, and projections from L3L4 to L5 are
off by ~7 windows (projection approximation at low pT).

## Checking `reachesRadius` before barrel projections

The repository code defines `Settings::reachesRadius(r, rinv)` (0.5 r |rinv| < 0.9) but
the gate before barrel projections in `TrackletCalculatorBase::barrelSeeding` is
commented out; the 2026-09-24 and 2026-09-25 v2 sample sets (identical hybrid
efficiency within 0.02 for QED μμ and e⁺e⁻) are both without it. The diagnostic build
switches it on with `TT_REACHGATE=0.9`. Same events, hybrid efficiency (barrel |η| < 0.8):

| | μμ 0.5–0.75 / 0.75–1 / 1–2 | e⁺e⁻ 0.5–0.75 / 0.75–1 / 1–2 | pPb 0.5–0.75 / 0.75–1 / 1–2 | pPb fakes, tracks/event |
|---|---|---|---|---|
| v2 cuts | 0.03 / 0.46 / 0.80 | 0.01 / 0.17 / 0.43 | 0.02 / 0.27 / 0.68 | 0.4%, 27.5 |
| v2 cuts + reach check | 0.03 / 0.49 / 0.80 | 0.01 / 0.17 / 0.43 | 0.02 / 0.27 / 0.68 | 0.4%, 27.5 |
| ME 6 + ×4 | 0.26 / 0.84 / 0.96 | 0.16 / 0.47 / 0.66 | 0.19 / 0.58 / 0.81 | 0.7%, 33.2 |
| ME 6 + ×4 + reach check | 0.27 / 0.83 / 0.96 | 0.16 / 0.47 / 0.66 | 0.19 / 0.59 / 0.81 | 0.7%, 33.2 |

(μμ: 2556 events; e⁺e⁻: 3759; pPb: 600.) No change beyond statistics. The check removes
26–38% of the L1L2 → L5/L6 and 30% of the L3L4 → L6 projections of true muon tracklets;
without it, about half of those are rejected later as out of φ/z range and the rest are
routed to matching but never match the true stub. It saves matching bandwidth (relevant
for truncation at high occupancy) and keeps projections finite; it does not recover
tracks. With the check on, the code also tries the disk projections those tracklets
skip today, which did not add efficiency here.

## Rerun

On an lxplus node (the area lives in node-local `/tmp`):

```bash
scram project CMSSW CMSSW_15_1_0_patch3   # inside cmssw-el8
# copy L1Trigger/{TrackFindingTracklet,TrackFindingTMTT,TrackTrigger,TrackerDTC} and
# SimTracker/TrackTriggerAssociation from the release, overlay this repo's files,
# apply the patches, scram b
bash run_one.sh <step2 index> <nevents>   # per job directory: truthtrace.txt, diag_ntuple.root
python3 analyze_trace.py <runs dir> <out dir>
python3 stage_report.py <out dir>
```
