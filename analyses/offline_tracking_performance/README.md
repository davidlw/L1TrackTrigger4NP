# Offline tracking performance

Efficiency, fake rate and duplicate rate for **offline** tracks (`generalTracks`)
matched to TrackingParticles.

Offline reconstruction is a single algorithm, so there is nothing to compare
against and every plot has one curve. The L1 equivalent lives in
`../L1_tracking_performance/`, which compares Default vs Dummy stub.

```
efficiency     = N(tp_nmatch > 0)   / N(tp)
fake rate      = N(trk_isTrue == 0) / N(trk)
duplicate rate = N(tp_nmatch > 1)   / N(tp_nmatch > 0)
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

## 2. Fill the histograms

Run from `scripts/`. **Give the path to your own ntuple** -- the defaults are
only placeholders.

```
cd scripts
root -l -b -q 'offline_perf.C("../output/OfflineTrackNtuple.root")'
```

One file is the common case. For many files, pass the **directory** plus how many
to read and the filename pattern:

```
root -l -b -q 'offline_perf.C("/path/to/dir","../output/offperf.root",50)'
root -l -b -q 'offline_perf.C("/path/to/dir","../output/offperf.root",50,false,"myntuple_%d.root")'
```

Arguments: `input`, `outname`, `nfiles` (0 = `input` is a single file),
`muonsOnly` (keep only `|pdgid| == 13`), `pattern`.

Integrated rates are printed to the terminal; the binned histograms go to
`outname`.

## 3. Plot

```
root -l -b -q 'plot_offline_perf.C("../output/offperf_qed_mumu.root")'
```

Writes PDF + PNG to `../figures/` and dumps the binned numbers as text. Optional
2nd/3rd arguments are a filename tag and the sample label drawn on each plot:

```
root -l -b -q 'plot_offline_perf.C("../output/offperf_hydjet.root","_hydjet","HYDJET PbPb, offline tracks")'
```

Panels: `off_{eff,fake,dup}_vs_{pt,eta,phi}` plus `_vs_z0` for efficiency and fake
rate, and `_pt1` variants of the eta/phi panels.

## Conventions

Binning, the |eta| < 2.4 acceptance and the pT > 1 / pT > 2 variants are identical
to the L1 macros, so an offline panel overlays the corresponding L1 one without
rescaling. The y axis is 0-1.35 for all three quantities, again as in the L1 plots.

The efficiency denominator carries **no hit or stub requirement** -- the TP
selection in the ntuplizer is deliberately loose, so efficiency is not biased by
the reconstruction it measures. The low-pT reach is set upstream by `ptMinTP` in
the step-2 mixing module, not here.

`output/` and `figures/` are gitignored.
