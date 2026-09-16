# Offline tracking performance

Efficiency, fake rate and duplicate rate for **offline** tracks (`generalTracks`)
matched to TrackingParticles, plus the event trigger efficiency.

Offline reconstruction is a single algorithm, so nothing is compared against
anything: within a panel the curves differ only by pT threshold. The L1
equivalent is `../L1_tracking_performance/`, which compares Default vs Dummy stub.

```
efficiency     = N(TP matched by >=1 selected track) / N(tp)
fake rate      = N(track not matched to a TP)        / N(selected tracks)
duplicate rate = N(TP matched by  >1 selected track) / N(matched TP)
```

## 1. Make the ntuple

```
cd $CMSSW_BASE/src/<path>/configurations
cmsRun offlineTrackNtuple_cfg.py
```

Edit `process.source.fileNames` to point at your step3 file. It must have been
produced with the extra keeps in `simulations/commands_cmsDriver` -- the matching
is hit-based and needs `generalTracks` + trackExtras + rec hits, plus
`siPixelClusters` and `siPhase2Clusters`. Check first:

```
edmDumpEventContent step3.root | grep -E "generalTracks|Cluster|TrackingParticle"
```

If the keeps are missing the job still exits 0 and every `trk_*` branch is empty.

No L1 re-emulation is involved, so this runs on step3 alone.

### On the grid

For the published step3 datasets use `configurations/crab_offlineTrackNtuple_cfg.py`:
pick a sample at the top, then `crab submit -c crab_offlineTrackNtuple_cfg.py`.

Submit from **CMSSW_15_1_0_patch3**, where the plugin is built -- CRAB ships the
local release area with the job, so the plugin must exist where you submit from.

## 2. Fill the histograms

Run from `scripts/`. **Give the path to your own ntuple** -- the defaults are
placeholders.

```
cd scripts
root -l -b -q 'offline_perf.C("../output/OfflineTrackNtuple.root")'
```

For many files pass the **directory** -- every `*.root` in it is read, whatever
the names -- how many to read, and the rest:

```
# input, outname, nfiles, muonsOnly, hpOnly, nchMax, tpClass
root -l -b -q 'offline_perf.C("/eos/.../0000","../output/offperf_hydjet_hp.root",900,false,true,8000,1)'
```

| argument | meaning |
|---|---|
| `nfiles` | how many files of a directory to read, 0 = all; ignored for a single file |
| `muonsOnly` | keep only `\|pdgid\| == 13` |
| `hpOnly` | keep only highPurity tracks |
| `nchMax` | upper edge of the Nch axis -- **sample dependent**, see below |
| `tpClass` | 0 all, 1 primary (`tp_ngenpart > 0`), 2 GEANT secondary |

`nchMax` measured on the current samples: QED µµ **4**, EPOS pPb **250**,
HYDJET PbPb **8000**. Entries above it go to overflow and vanish from the plot,
so set it generously. Below 40 the macro gives one bin per unit.

## 3. Plot

```
root -l -b -q 'plot_offline_perf.C("../output/offperf_hydjet_hp.root","_hydjet_hp","HYDJET PbPb, offline highPurity")'
```

Arguments are the histogram file, a filename tag, and the label drawn on each
plot. Writes PDF + PNG to `../figures/` and dumps binned numbers as text.

Panels: `off_{eff,fake,dup}_vs_{pt,eta,phi,nch}`. The pT panels use a **log x
axis** from 0.3 to 10 GeV; the eta, phi and Nch panels **overlay four pT
thresholds** (0.3 / 0.6 / 1 / 2 GeV).

## 4. Event trigger efficiency

Fraction of events with at least one object above a pT threshold -- primary truth
particles against offline highPurity tracks, so the gap between the curves is
what reconstruction costs relative to a perfect tracker.

```
# input, outname, nfiles, nchMax
root -l -b -q 'mbeff_offline.C("/eos/.../0000","../output/mbeff_hydjet.root",900,6000)'
root -l -b -q 'plot_mbeff_offline.C("../output/mbeff_hydjet.root","_hydjet","HYDJET PbPb")'
```

The denominator is every event and is the same at every threshold, so events with
nothing above it stay in the denominator. Writes `mbeff_vs_ptmin<tag>.pdf`.

## 5. Trigger rate

The same histograms, scaled by an assumed total collision rate:

```
root -l -b -q 'plot_rate_offline.C("../output/mbeff_hydjet.root","_hydjet","HYDJET PbPb",50e3)'
root -l -b -q 'plot_rate_offline.C("../output/mbeff_epospb.root","_epospb","EPOS pPb",3e6)'
```

The last argument is the total rate in **Hz** and is an assumption you supply --
the shape comes from the MC, the normalisation from you. Writes
`rate_vs_ptmin<tag>.pdf` (require >= 1 track above a pT threshold) and
`rate_vs_nchmin<tag>.pdf` (require Nch >= N), both with a log rate axis.

Note the reco curve sits **below** truth on the Nch plot: reconstruction misses
tracks, so measured multiplicity is systematically lower and a given threshold is
harder to reach. Selecting the same events needs a lower cut on reconstructed
multiplicity than on true multiplicity.

## Conventions

Everything is inside **|eta| < 2.4**, the pT spectra included.

The efficiency denominator carries **no hit or stub requirement** -- the TP
selection in the ntuplizer is deliberately loose, so efficiency is not biased by
the reconstruction it measures. The low-pT reach is set upstream by `ptMinTP` in
the step-2 mixing module (0.3 GeV in the current samples), not here.

**Primaries vs secondaries matters.** About a quarter of the stored TPs in the
heavy-ion samples are GEANT secondaries, reconstructed at ~14% against ~84% for
primaries, so the all-particle efficiency sits far below the primary one. Run
`tpClass = 1` and `2` into separate files to see the two populations.

**Every truth quantity is primaries-only**, regardless of `tpClass`: the Nch axis,
the truth curve of the trigger efficiency, and the truth curve of the rate plots
all require `tp_ngenpart > 0`. `tpClass` selects which particles the efficiency
and duplicate rate are measured *for*; it does not change how Nch is counted.
Multiplicity and the trigger ceiling should describe the collision, not how much
material a secondary happened to traverse.

With `tpClass = 1` a track matching only a secondary counts as **fake** -- the
usual convention when the denominator is primaries. The primary fake rate is
therefore higher than the inclusive one, and the difference is the secondary
contamination. In the `tpClass = 2` output the fake-rate panels are meaningless
by construction (a track is "fake" unless it matches a secondary, and most match
primaries); only efficiency and duplicate rate mean anything there.

Because `matchtrk_*` has no highPurity field, `tp_nmatch` cannot be used once a
track selection is applied -- it counts every match regardless of quality. The
macro instead rebuilds the per-TP count by joining tracks to TPs on their truth
values, which are identical floats on both sides. Validated against `tp_nmatch`
with no selection: identical over 31301 TPs, duplicates included.

The y axis is 0-1.35 for all three quantities, as in the L1 plots.

`output/` and `figures/` are gitignored.
