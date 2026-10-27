# GNN edge classification for L1 track finding

An interaction-network (IN) edge classifier on outer-tracker stubs, following
DeZoort et al., *Charged particle tracking via edge-classifying interaction
networks*, Comput. Softw. Big Sci. 5 (2021) 26,
[arXiv:2103.16701](https://arxiv.org/abs/2103.16701), trained and tested on the
`L1TrackHitNtuple` files (`allstub_*` branches). Tracks are built from the
classified edges and matched to TrackingParticles.

Plain PyTorch, no CMSSW or ROOT needed (the ntuples are read with uproot).

## Pipeline

| step | script | what it does |
|---|---|---|
| 1 | `build_graphs.py` | ntuple → graphs: nodes are stubs, candidate edges join stubs on allowed layer pairs inside Δφ/Δη windows derived from truth |
| 2 | `train.py` | trains the IN (BCE loss), sets the working point δ where TPR = TNR on the validation set |
| 3 | `evaluate.py` | edge metrics, then track building and matching; efficiency / fake / duplicate rate vs pT, η, threshold |
| 4 | `compare_hybrid.py` | GNN vs the hybrid L1 tracks (`trk_*`) in the same ntuple events, on the same particles |

`gnn_common.py` (ntuple reading, graph shards) and `model.py` (the IN) are shared.

## Setup

```bash
python3 -m venv ~/venvs/gnn
~/venvs/gnn/bin/pip install -r requirements.txt
```

## Usage

Run from `scripts/`. Inputs follow `analyses/common/InputFiles.h`: a file, a
directory (every `*.root` in it) or a comma-separated list; `--nfiles N` takes
the first N files of a directory.

```bash
cd analyses/gnn_tracking/scripts
python3 build_graphs.py --input /my/ntuples/DummyStub/EPOS_pPb --nfiles 2 --out ../output/graphs_ppb
python3 train.py --graphs ../output/graphs_ppb --out ../output/model_ppb --epochs 20
python3 evaluate.py --checkpoint ../output/model_ppb/model.pt
```

Test a trained model on another sample by building its graphs with the same
windows and passing them to `evaluate.py` (all events are used, since none were
seen in training):

```bash
python3 build_graphs.py --input /my/ntuples/DummyStub/HYDJet_PbPb --out ../output/graphs_hydjet \
    --config ../output/graphs_ppb/graph_config.json
python3 evaluate.py --checkpoint ../output/model_ppb/model.pt --graphs ../output/graphs_hydjet
```

Compare with the hybrid L1 tracks of the same events (and, in a third panel,
the hybrid efficiency of another sample with the standard `tp_nmatch`
definition of `compare_eff.C`):

```bash
python3 compare_hybrid.py --checkpoint ../output/model_ppb/model.pt --threshold 0.9 \
    --reference /my/ntuples/DefaultStub/EPOS_pPb
```

Outputs (`output/` is not tracked): `graph_config.json` with the layer pairs and
windows, `graphs_*.npz` shards, `model.pt`, `history.json`, `training.png`, and
per evaluated sample `eval_<graphs>/summary.json` and `performance.png`.

## Details

### Graph construction (`build_graphs.py`)

- **Nodes:** stubs, with features r, φ, z (default), and optionally bend, isPS,
  isBarrel (`train.py --node-features r,phi,z,bend`). DummyStub samples have no
  bend (`trigBend = 999999`), so bend is only useful for DefaultStub samples.
- **Layers:** barrel 1–6, endcap disks 11–15 (+z) and 21–25 (−z). The ntuple's
  `allstub_layer` does not encode the disk side, so the sign of z is used.
- **Particle identity:** the ntuple stores no TP index per stub, so a particle
  is the unique `(pdgid, pt, eta, phi)` of the matched TP.
- **Candidate edges:** all stub pairs on an allowed layer pair with
  |Δφ| < Δφ_max and |Δη| < Δη_max, oriented inner → outer. Edge features:
  Δr, Δφ, Δz, Δη, ΔR (the paper uses Δr, Δφ, Δz, ΔR).
- **Windows from truth, not hand-tuned.** Each particle's stubs are ordered by
  distance from the origin; pairs on consecutive layers and one layer apart
  define which layer pairs are allowed (`--min-pair-count`). The windows are
  the 99.5% quantile (`--quantile`) × 1.1 (`--safety`) of each particle's
  closest pair on that layer pair. Two choices matter:
  - *Closest pair per particle*, not every pair: one surviving pair keeps a
    particle connected, and every-pair statistics are dominated by electrons
    with tens of stubs on one TP.
  - *On-helix stubs only* (`--helix-tol`, 0.05 rad in φ from the production
    helix): 98% of muon stubs pass but only 78% of electron stubs; the rest are
    bremsstrahlung or looping segments of the same TP. Without this the L1–L2
    window is ~2 rad instead of 0.13 rad.
- **Physics preconditions** (after the windows), each with a per-layer-pair
  tolerance measured on prompt particles:
  - `--pt-min` (default 0.5 GeV): the circle through the beamline (d0 = 0)
    and the two stubs must have |κ| ≤ 0.3·B / pT_min. On true pairs this pair
    pT has a −1% median bias; 90% of barrel pairs are within 6% of the TP pT
    (19% when a disk is involved).
  - `--z0-max` (default 15 cm): stub b must lie inside the cone of lines from
    the luminous region |z0| ≤ Z through stub a, in the (s, z) plane with s the
    arc length on that circle. This is used instead of extrapolating z0,
    which blows up for disk-to-disk pairs where Δs is a few cm, comparable to
    the 5 cm 2S strip length along r. 97% of true prompt barrel PS pairs sit
    exactly inside the cone (tolerance 0); disk pairs need 20–400 cm.

  The reference z0 of a particle is the median of z − s·sinh η_TP over its PS
  stubs (the DummyStub QED files have no `tp_*` branches to join).
  With `--config`, the preconditions of the config are used unless given
  explicitly. The pair's z0 (clipped to ±100 cm) and q/pT are also stored as
  edge features `z0`, `curv` (`train.py --edge-features ...,z0,curv`).
- **Labels:** an edge is true when both stubs belong to the same particle. Every
  such pair is true (not just consecutive stubs), so connected components of
  the true edges are exactly the particles. This matters because DummyStub
  gives several stubs per layer per particle.
- **Unmatched stubs are dropped by default** (the paper's idealized hit
  filtering). With `ptMinTP = 1 GeV` in step 2 (see the top-level README), a
  stub with no TP is often a real particle below 1 GeV; keeping it
  (`--keep-unmatched`) would label real low-pT tracks as noise.

### Model (`model.py`)

The paper's IN (Fig. 6): relational model φ_R1 → sum over incoming edges →
object model φ_O → edge-weight model φ_R2, each an MLP with two hidden layers
of 40 units, 4-dimensional messages and a 3-dimensional node update. 6488
parameters (the paper's 6448 plus the extra Δη input). Graphs are made
bidirectional, and an edge's weight is the minimum of its two directions, as in
the paper's track-building code. `--iterations N` repeats the message-passing
step (the paper uses 1).

### Training (`train.py`)

Adam, learning rate 5×10⁻³ decayed by 0.95 every 10 epochs, BCE on both edge
directions, 64 events per batch. Events are split 80/10/10 (train /
validation / test) with a fixed seed; the split is saved in the checkpoint.
The checkpoint with the lowest validation loss is kept. δ is the threshold at
TPR = TNR on the validation set, as in the paper.

### Track building and matching (`evaluate.py`)

1. Keep edges with weight > δ (`--threshold` to override).
2. Connected components of the kept edges (union-find).
3. A component with ≥ 3 stubs (`--min-nodes`) is a track candidate.
4. A candidate is a true track if ≥ 2/3 of its stubs (`--min-purity`) come from
   one particle; that particle is found. A second true track of an already
   found particle is a duplicate.

**One crossing per layer** (`--crossing-dist 5`, `--max-layers 8`): a track
that does not loop crosses each layer once and, in the barrel–endcap transition,
spans at most 8 layers (2 barrel + 5 disks + one more; 11% of non-looping
muons at 1.5 < |η| < 1.8 have 8). Stubs on one layer within 5 cm of each other
(single linkage) are one crossing: on-helix stubs of one muon on a layer are
within 5 cm in 64% of pairs, stubs of different muons in 1%. Kept edges are
added by decreasing weight, and a union that would give two crossings on a
layer or more than 8 layers is refused, so merged tracks are split instead of
lost.

Efficiency = found / particles with ≥ 3 stubs in the graph. Fake rate =
candidates that are not true tracks / all candidates.

The same track building is run with the truth labels as edge weights ("ideal
classifier"); the gap to 100% is what graph construction loses. A scan of the
edge-weight threshold is also printed and plotted.

## First results

EPOS pPb, DummyStub, 2 files (10k events: 140 stubs and 715 candidate edges per
event after dropping unmatched stubs, 49% of edges true), without the physics
preconditions. Paper IN, 20 epochs, test split of 1000 events, 18.6k particles
with ≥ 3 stubs:

| | efficiency | fake rate | duplicate rate |
|---|---|---|---|
| GNN at δ = 0.49 (TPR = TNR) | 0.81 | 0.072 | 0.008 |
| GNN at threshold 0.8 | 0.87 | 0.032 | 0.010 |
| GNN at threshold 0.9 | 0.87 | 0.018 | 0.012 |
| ideal classifier | 0.94 | 0 | 0.004 |

Edge AUC 0.990, accuracy 95.6% at δ (the paper reaches 99.9% with one pixel
hit per layer). The TPR = TNR point is not the best one for track building:
a false edge merges two tracks and both fail the purity cut, so a tighter
threshold raises the efficiency and lowers the fake rate together. Efficiency
is flat in pT above 1 GeV and lowest in the barrel (|η| < 1), where the paper
also finds most errors.

Three message-passing steps (`--iterations 3`, 15k parameters) help a little:
AUC 0.9929, and at threshold 0.9 an efficiency of 0.878 with a 0.015 fake rate.

### Effect of the physics preconditions

Same events, graphs built with the same windows, one training each (20 epochs;
retraining the same configuration moves the efficiency by about 0.3%):

| preconditions | edges / event | graph purity | ideal eff. | GNN eff. / fake at δ | GNN eff. / fake at 0.9 |
|---|---|---|---|---|---|
| none | 715 | 0.487 | 0.942 | 0.811 / 0.070 | 0.865 / 0.020 |
| pT > 0.5 GeV, \|z0\| < 15 cm | 649 | 0.522 | 0.935 | 0.804 / 0.070 | 0.866 / 0.023 |
| pT > 1.0 GeV, \|z0\| < 15 cm | 609 | 0.552 | 0.932 | 0.753 / 0.094 | 0.839 / 0.029 |
| pT > 0.5, \|z0\| < 15, + `z0,curv` edge features | 649 | 0.522 | 0.935 | 0.822 / 0.061 | 0.867 / 0.021 |

With all stubs kept (`--keep-unmatched`, 1395 stubs per event, 1000 events),
the graph has 22.7k edges per event and 1.5% purity; the preconditions remove
20% (pT > 0.5) or 32% (pT > 1) of the edges and 2–3% of the true ones.

The preconditions shrink the graph (what the FPGA must hold) but do not improve
the tracking: the truth-derived windows already contain most of the pT and
vertex constraint, and the IN learns the rest from Δr, Δφ, Δz. The edges it gets
wrong are pairs from nearby particles that *are* compatible with a prompt
track; no two-stub cut can remove those. pT > 1 GeV loses more, because 2.6% of
true edges near threshold fall outside the tolerance and each lost edge can
split a track.

### Comparison with the hybrid L1 tracks

Same 1000 pPb test events and 18.6k particles (≥ 3 DummyStub stubs), GNN with
the default preconditions at an edge cut of 0.9. A particle counts as found by
the hybrid if a `trk_*` track is loosely genuine to it (≥ 75% purity for its
≥ 4 stubs); by the GNN if a candidate has ≥ 2/3 purity.

| efficiency | particles | GNN, ≥ 3 stubs | GNN, ≥ 4 stubs | hybrid |
|---|---|---|---|---|
| pT 1–2 GeV | 14 959 | 0.870 | 0.829 | 0.020 |
| pT > 2 GeV | 3 433 | 0.879 | 0.842 | 0.769 |

With the standard definition (TPs with `tp_nmatch > 0`, each sample on its own
TPs), the hybrid reaches 0.829 above 2 GeV on these DummyStub files and 0.867
on DefaultStub, so the DummyStub hybrid is about 4% below the default L1
tracking. Read the comparison with these differences in mind:

- **The GNN sees only truth-matched stubs** (140 per event); the hybrid sees all
  1424. Its fake rates are therefore not comparable (GNN 2.3% of candidates,
  hybrid 0.7% of tracks), and its efficiency is an upper bound for a GNN on
  real inputs.
- Below 2 GeV the hybrid is off by configuration (its pT threshold is
  1.95 GeV), not by design; the fair low-pT baseline is the hybrid with that
  threshold lowered.
- The hybrid needs ≥ 4 stubs; the ≥ 4-stub GNN column is the like-for-like one.
- The 0.9 edge cut was chosen from the threshold scan on these test events.

### Low-pT QED μμ with the extended hybrid (v2 samples)

STARlight QED μμ, `L1TrackNtuplev2` production (truth down to 0.3 GeV, hybrid
tracking extended to low pT, DummyStub). 50 files, 155k events; 20 stubs per
event, 98% of them truth-matched, so all stubs are kept (`--keep-unmatched`)
and the GNN sees the same stubs as the hybrid.

```bash
python3 build_graphs.py --input /my/ntuples/QED_mumu_v2 --out ../output/graphs_mumu_v2 \
    --keep-unmatched --pt-min 0.3
python3 train.py --graphs ../output/graphs_mumu_v2 --out ../output/model_mumu_v2
python3 compare_hybrid.py --checkpoint ../output/model_mumu_v2/model.pt --threshold 0.9 \
    --count layers --denominator tp --eta-pt-min 1.0 \
    --offline ../../offline_tracking_performance/output/offperf_qedmumu_hp_prim_mu.root \
    --offline-ntuple /my/ntuples/OfflineTrackNtuple_QEDMuMu_v1
```

`gnn_vs_hybrid.png` has efficiency vs pT, η and φ (the last two above
`--eta-pt-min`) with the offline efficiency as a line; `eff2d.png` has the
efficiency maps in (η, pT) for offline (from `--offline-ntuple`, the
`offline_perf.C` definition for primary muons), hybrid and GNN, and their
differences.

Efficiency on the 15.5k test events, denominator all TPs with |η| < 2.4
(`compare_eff.C` definition; hybrid found = `tp_nmatch > 0`). GNN candidates
are counted in distinct layers, since the hybrid needs 4 stubs on 4 layers
and DummyStub gives several stubs per layer:

| TP pT [GeV] | TPs | GNN ≥ 3 layers | GNN ≥ 4 layers | no classifier, ≥ 4 layers | hybrid |
|---|---|---|---|---|---|
| 0.3–0.5 | 21 390 | 0.896 | 0.669 | 0.446 | 0.047 |
| 0.5–0.75 | 5 395 | 0.984 | 0.949 | 0.776 | 0.382 |
| 0.75–1 | 1 534 | 0.998 | 0.993 | 0.988 | 0.665 |
| 1–2 | 1 187 | 0.997 | 0.991 | 0.996 | 0.901 |
| > 2 | 196 | 1.000 | 1.000 | 0.995 | 0.995 |

Fake rate: GNN 1.1% of candidates at an edge cut of 0.9 (22% with no
classifier); hybrid 0%. Edge AUC 0.945.

- **Above 0.75 GeV the gain over the hybrid is not the classifier.** The same
  connected components with every edge kept ("no classifier") already find
  99%: the hybrid loses these tracks in its own windows, seeding or fit, not to
  combinatorics.
- **Below 0.75 GeV the classifier matters.** The windows needed at 0.4 GeV
  connect the two muons: a back-to-back μ⁺μ⁻ pair with equal pT (STARlight
  pairs have almost no pair pT) lies on the same circle in x–y, so only z/η
  separates them. The IN raises the efficiency from 0.78 to 0.95
  (0.5–0.75 GeV) and cuts the fake rate from 22% to 1%.
- **The extended hybrid fails in the barrel.** For pT > 0.5 GeV it is ~95%
  efficient for |η| > 1.5 but ~25% for |η| < 1, although nearly every muon
  leaves stubs on ≥ 4 layers there. That points inside the hybrid chain
  (barrel seeding, projection windows or fit at low pT).
- The GNN candidates are stub clusters without a track fit; a 7% duplicate
  rate comes from looping muons whose outgoing and returning arcs form separate
  pure components.
- **41% of the stubs attributed to the muons are off their helix**, 70% of
  those in PS modules: looping low-pT muons (a central 0.35 GeV muon turns at
  r = 61 cm and advances ~40 cm in z per turn), and for muons that do not loop,
  most likely delta rays whose hits are assigned to the parent muon.

With one crossing per layer (≥ 4 layers; GNN at an edge cut of 0.5, the
no-classifier variant adds edges by increasing ΔR), same events:

| TP pT [GeV] | GNN | GNN + 1 crossing/layer | no classifier | no classifier + 1 crossing/layer | hybrid |
|---|---|---|---|---|---|
| 0.3–0.5 | 0.669 | 0.765 | 0.446 | 0.686 | 0.047 |
| 0.5–0.75 | 0.949 | 0.968 | 0.776 | 0.948 | 0.382 |
| 0.75–1 | 0.993 | 0.995 | 0.988 | 0.996 | 0.665 |
| fake rate | 0.9% | 1.5% | 24% | 9.6% | 0 |
| duplicate rate | 0.8% | 4.8% | 0.7% | 6.4% | – |

The constraint alone recovers most of the merged tracks (no classifier:
0.78 → 0.95 at 0.5–0.75 GeV) but leaves a 9.6% fake rate; with the GNN it gives
the best efficiency at every pT with 1.5% fakes, close to offline
(`generalTracks` highPurity, different events) down to 0.35 GeV. Split looper
arcs and delta-ray stubs raise the duplicate rate to ~5%.

### Low-pT QED e⁺e⁻ (v2)

Same procedure as for μμ (`--keep-unmatched --pt-min 0.3`; 234k events, 26 stubs per
event of which 58% are truth-matched; 44% of candidate edges are true, edge AUC 0.945).
Test split of 23.4k events, all TPs with |η| < 2.4, ≥ 4 layers; GNN with one crossing
per layer at an edge cut of 0.7; offline = highPurity `generalTracks`, primary
electrons, from the offline ntuples (different events):

| TP pT [GeV] | TPs | offline | GNN + 1 crossing/layer | no classifier + 1 crossing/layer | hybrid |
|---|---|---|---|---|---|
| 0.3–0.5 | 5 230 | ~0.6 | 0.35 | 0.50 | 0.03 |
| 0.5–0.75 | 7 185 | ~0.8 | 0.67 | 0.70 | 0.17 |
| 0.75–1 | 7 897 | ~0.87 | 0.82 | 0.80 | 0.31 |
| 1–2 | 13 809 | ~0.9 | 0.88 | 0.82 | 0.46 |
| > 2 | 2 389 | ~0.93 | 0.88 | 0.82 | 0.64 |
| fake rate | | | 1.8% | 36% | 0.2% |

- The hybrid loses electrons at every pT (0.64 above 2 GeV, against 0.93 offline),
  worst in the barrel: bremsstrahlung on top of the low-pT matching losses.
- The GNN is within a few percent of offline above 0.75 GeV everywhere in η, and
  above offline at 0.3–0.6 GeV for 0.5 < |η| < 2. It falls short at |η| < 0.4 below
  0.6 GeV, where electrons reach only L1–L3.
- Unlike μμ, the classifier is essential here: without it the same constrained
  linking has 36% fakes.

## Known limitations

- **No truth below 1 GeV** in the first samples (`ptMinTP = 1 GeV`); the v2
  productions (`L1TrackNtuplev2`) go down to 0.3 GeV.
- **QED e⁺e⁻ (DummyStub) is nearly trivial** after dropping unmatched stubs:
  4.3 stubs per event and 99.5% of candidate edges true.
- **DefaultStub EPOS pPb has no stub truth** (every `allstub_matchTP_pdgid`
  is −999), so it cannot be used for training.
- Muons in pPb are few (69 in the test split, mostly decay products), so their
  efficiency is not meaningful there; use the QED μμ sample.
