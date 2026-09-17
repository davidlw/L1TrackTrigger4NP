# Tracking performance macros: Default Stub vs Dummy Stub

Two passes: `compare_*.C` reads the ntuples and writes a small `.root`; `plot_*.C`
reads that and draws figures. The split means you can re-plot without re-reading the
input files.

**Run everything from `../output/`** — figures go to `../figures/`, so the relative
paths only resolve from there.

## Quick start

```bash
cd analyses/L1_tracking_performance/output

D=/my/path/DefaultStub/QED_mumu; U=/my/path/DummyStub/QED_mumu    # directories or single files
root -l -b -q "../scripts/compare_eff.C(50,true,\"eff_qed_mumu.root\",\"$D\",\"$U\")"    # efficiency + duplicate rate
root -l -b -q '../scripts/plot_eff.C("eff_qed_mumu.root")'
root -l -b -q "../scripts/compare_fake.C(50,\"fake_qed_mumu.root\",\"$D\",\"$U\")"        # fake rate
root -l -b -q '../scripts/plot_fake.C("fake_qed_mumu.root")'
```

There is no built-in sample: **the input locations are always arguments.** A
macro run without them stops with a message.

## Setting the input paths

Each `compare_*` macro reads two productions, Default stub and Dummy stub. Each
location is either a **directory** -- every `*.root` in it is read, sorted by
name, and `nfiles` caps how many (0 = all) -- or a **single file**. File names do
not matter, so three files called `_4`, `_17`, `_22` copied to a laptop work as
well as a full `/eos/.../0000`. What was found is printed at the start of the run.

```bash
root -l -b -q '../scripts/compare_eff.C(100,false,"eff_mine.root",\
  "/my/path/DefaultStub/HYDJet_PbPb",\
  "/my/path/DummyStub/HYDJet_PbPb")'
```

`compare_eff` takes a trailing `pdgSel` (11 electrons, 13 muons, 211 pions; 0 = all)
after `dirDum`; `muonsOnly` is the old spelling of 13. On the STARlight QED samples
with default stubs the L1 efficiency at 2-2.5 GeV is 0.52 for electrons against 0.96
for muons -- bremsstrahlung breaks the single-helix assumption -- so e+e- is not
"mu+mu- with another label".

Macros that pair the two productions event by event (`compare_res`, and
`compare_stubmult` with `commonIndicesOnly`) keep only the file names present on
both sides. The file lister is shared: [`../../common/InputFiles.h`](../../common/InputFiles.h).

**Assume the ntuples live on lxplus**, under `/eos/cms/store/group/phys_heavyions/`.
Run the `compare_*` pass there; it is the only step that touches the samples. The
`.root` it writes is tens of kB, so copy that back and run `plot_*` locally --
axis ranges and styling can then be reworked in seconds without re-reading TB of
input. `compare_stubs.C` already defaults to EOS paths; the older macros still
default to local copies and need their paths passed.

```
efficiency     = N(tp_nmatch > 0)    / N(tracking particles)
duplicate rate = N(tp_nmatch > 1)    / N(tp_nmatch > 0)
fake rate      = N(trk_fake == 0)    / N(tracks)
non-genuine    = N(trk_genuine == 0) / N(tracks)      (looser than fake)
```


## Stub distributions

`compare_stubs.C` / `plot_stubs.C` compare the two productions at stub level,
splitting ACCEPTED from REJECTED via `allstub_isRejected`.

```bash
# on lxplus -- the only step that reads the samples
root -l -b -q 'compare_stubs.C(5,"<defDir>","<dumDir>","../output/stub_epos.root")'
# anywhere -- instant to redraw
root -l -b -q 'plot_stubs.C("../output/stub_epos.root","_epos","EPOS pPb","epos")'
```

Panels: stubs per event (the main one, log-log, four curves), plus r, z,
layer/disk and trigBend. The multiplicity is stored at one bin per stub, so the
plot step can pick any binning or upper limit.

`plot_stubs.C(file, tag, label, pfx, xmax, binsPerDecade, minWidth, logx)` -- `xmax` 0
takes 25% above the largest count seen; `binsPerDecade` (12) and `minWidth` (1)
control the multiplicity binning, which is geometric with a floor on the width.
**Raise `minWidth` when the low-N end looks spiky.** In the heavy-ion samples
almost every event carries thousands of stubs, so that region holds only a
handful of events and width-1 bins there show single counts. HYDJET reads well
at `minWidth = 4`. Redrawing is instant, so tune it.

`logx` (default true) switches the multiplicity panel to a linear axis and to
uniform bins. Log is the better default here: the four curves span a large range
(HYDJET ⟨N⟩ runs 33 / 1020 / 5454 / 38763), so on a linear axis wide enough for
the dummy tail everything else collapses into the first bin or two. Linear is
worth it for the *shape* of one high-multiplicity distribution -- the dummy
accepted curve is a broad plateau with a falling edge, which log compresses away
-- usually together with an `xmax` to zoom.

Measured, EPOS pPb / HYDJET PbPb: dummy produces **8.9x / 7.3x** more accepted
stubs than default (1408 vs 158 per event, and 39874 vs 5454).

### Which production is which

A dummy production always has `allstub_trigBend` at the 999999 sentinel. What
varies between productions is whether it rejects anything at all:

| sample | dummy | rejected/evt |
|---|---|---|
| EPOS pPb | `260825_213152` | 0 |
| EPOS pPb | `260828_193234` | 9.5 (0.68%) |
| HYDJET PbPb | `260825_213428` | 0 |
| HYDJET PbPb | `260828_193408` | 1020 (2.6%) |

**An empty rejected curve means that production applied no stub selection**, not
that dummy stubs cannot be rejected. In the `260825` EPOS dummy every inner-sensor
cluster became exactly one stub (`dummy stubs / default sensor-0 clusters =
1.0000` over 25000 events), so nothing was paired and no window was applied and
no selection existed to fail. The `260828` productions do apply one. Prefer the
newer ones unless you specifically want the no-selection variant.

## Offline reference on the efficiency plots

`plot_eff.C` takes an optional 4th argument: a histogram file from
`../../offline_tracking_performance`. Given one, it overlays the offline
efficiency as a dashed black line -- useful as the ceiling L1 is working toward.

```bash
root -l -b -q 'plot_eff.C("../output/eff_qed_mumu_v2.root","_v2","STARlight QED #mu#mu",\
  "../../offline_tracking_performance/output/offperf_qedmumu_hp_prim.root")'
```

eta and phi overlay bin-for-bin; the pT binning differs, which is why the
reference is drawn with `HIST` rather than interpolated. The z0 panel gets no
reference, since offline has no matching pT > 2 selection. Check whether the two
files come from the same production before reading the gap quantitatively.
