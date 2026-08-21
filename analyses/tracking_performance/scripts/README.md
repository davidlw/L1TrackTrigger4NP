# Tracking performance macros: Default Stub vs Dummy Stub

Two passes: `compare_*.C` reads the ntuples and writes a small `.root`; `plot_*.C`
reads that and draws figures. The split means you can re-plot without re-reading the
input files.

**Run everything from `../output/`** — figures go to `../figures/`, so the relative
paths only resolve from there.

## Quick start

```bash
cd analyses/tracking_performance/output

root -l -b -q '../scripts/compare_eff.C(50)'                    # efficiency + duplicate rate
root -l -b -q '../scripts/plot_eff.C("eff_qed_mumu.root")'
root -l -b -q '../scripts/compare_fake.C(50)'                   # fake rate
root -l -b -q '../scripts/plot_fake.C("fake_qed_mumu.root")'
```

That uses the built-in sample. **For your own files, set the paths first.**

## Setting the input paths

Each `compare_*` macro reads two directories, one per production. Pass them as
arguments — every macro takes `dirDef` / `dirDum`, and an empty string keeps the
default:

```bash
root -l -b -q '../scripts/compare_eff.C(100,false,"eff_mine.root",\
  "/my/path/DefaultStub/HYDJet_PbPb",\
  "/my/path/DummyStub/HYDJet_PbPb")'
```

Or edit `kDirDefault` / `kDirDummy` at the top of the macro. The built-in defaults are
machine-specific and will not work elsewhere.

```
efficiency     = N(tp_nmatch > 0)    / N(tracking particles)
duplicate rate = N(tp_nmatch > 1)    / N(tp_nmatch > 0)
fake rate      = N(trk_fake == 0)    / N(tracks)
non-genuine    = N(trk_genuine == 0) / N(tracks)      (looser than fake)
```
