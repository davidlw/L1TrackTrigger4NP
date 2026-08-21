# Tracking performance: Default Stub vs Dummy Stub

ROOT macros comparing L1 track-finding performance between the two ntuple
productions — efficiency, fake rate, duplicate rate, momentum resolution, track
quality, and stub- and trigger-level distributions.

Works on any sample with these ntuples: STARlight QED mumu, HYDJet PbPb, EPOS pPb,
pythia pp.

## Layout

| path | contents |
| --- | --- |
| `scripts/` | the macros, plus [`scripts/README.md`](scripts/README.md) — **start there** |
| `output/` | `.root` files from the `compare_*` pass, and run logs. Run the macros from here. |
| `figures/` | every figure, written as both PDF and PNG |

## Getting started

```bash
cd output
root -l -b -q '../scripts/compare_eff.C(50)'                  # efficiency + duplicate rate
root -l -b -q '../scripts/plot_eff.C("eff_qed_mumu.root")'
root -l -b -q '../scripts/compare_fake.C(50)'                 # fake rate
root -l -b -q '../scripts/plot_fake.C("fake_qed_mumu.root")'
```

The built-in input paths are machine-specific — pass `dirDef` / `dirDum` for your own
files. [`scripts/README.md`](scripts/README.md) covers that, every performance
variable with its command, and the branch pitfalls worth knowing before you trust a
number.

## Results are deliberately not kept here

Every number depends on which sample, which production and how many files you run
over, so quoting results in a README would go stale and mislead. Run the macros and
read your own output.
