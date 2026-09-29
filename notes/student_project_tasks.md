# Senior thesis project — task list

Building on Nicole's low-pT configuration. Ordered by priority, and grouped so
that each phase produces something usable before the next one starts.

**Suggested thesis spine:** *Fixed-point precision and occupancy limits of a
low-pT Level-1 track trigger for heavy-ion collisions at the HL-LHC.* Tasks
T2, T4 and T5 are the core; T3 and T6 make it a systems result rather than a
numerical one.

Two things gate much of this and are **not** the student's to fix:

- samples regenerated with `ptMinTP = 0.3` (the current standard-configuration
  productions were made at 1.0, so below 1 GeV there is no denominator)
- one frozen low-pT configuration that everyone agrees to call final

Tasks are marked **[free]** if they need neither.

---

## Phase 0 — onboarding (weeks 1–2)

### T0. Reproduce an existing result  **[free]**

Install ROOT locally, clone the repository, and remake one figure that already
exists, without changing anything.

```
analyses/offline_tracking_performance/scripts/offline_perf.C
analyses/offline_tracking_performance/scripts/plot_eff_electron_quality.C
```

Ntuples are on this Mac under `~/Documents/OfflineNtuple/` (QED ee, QED μμ,
J/ψ μμ — 44–70 MB each) and the L1 ones under `~/Documents/{Default,Dummy}Stub/`.
The macros take the input directory as an argument and read every `*.root` in
it, so nothing needs to be renamed.

*Deliverable:* the electron-vs-muon efficiency figure, regenerated.
*Learns:* the repository layout, the efficiency/fake/duplicate definitions,
that "primary" and "highPurity" are choices and not details.

---

## Phase 1 — the projection arithmetic (weeks 2–6)

### T1. Layer reach versus pT  **[free]**

Purely analytic, no samples. A charged particle of transverse momentum pT in
the 3.8 T field has radius R = 87.5·pT cm and reaches a maximum distance 2R
from the beam line. Plot that against the six barrel layer radii (24.9, 37.2,
52.3, 68.7, 86.0, 108.3 cm).

*Deliverable:* figure F2 of the DP note list. One afternoon.
*Learns:* why below 0.62 GeV a track cannot produce six stubs at all, which is
the fact the whole project rests on.

### T2. How wrong is the projection at low pT?  **[free]**

**Do this first. It is arithmetic, needs no data, and the answer decides whether
the rest of the project is tuning or repair.**

The emulation projects a track to a layer with a third-order expansion of the
arcsine (`IMATH_TrackletCalculator.h:250-258`):

```
phi(r) = phi0 - x·(1 + x²/6)      where x = r·rinv/2      [exact: phi0 - asin(x)]
dphi/dr = -rinv/2                                         [exact: -rinv/2 / sqrt(1-x²)]
```

Both are excellent at 2 GeV and both fail as x → 1. Expressed in millimetres at
the layer, against match windows of 1–3 mm:

| pT [GeV] | L1 | L2 | L3 | L4 | L5 | L6 |
|---|---|---|---|---|---|---|
| 2.0 | 0.00 | 0.00 | 0.00 | 0.02 | 0.06 | 0.24 |
| 1.5 | 0.00 | 0.00 | 0.01 | 0.07 | 0.26 | 1.07 |
| 1.0 | 0.00 | 0.01 | 0.10 | 0.52 | 2.13 | **9.5** |
| 0.8 | 0.00 | 0.04 | 0.31 | 1.69 | **7.2** | **36** |
| 0.6 | 0.01 | 0.17 | 1.40 | **8.3** | **41** | — |
| 0.5 | 0.04 | 0.43 | **3.8** | **25** | **201** | — |
| 0.4 | 0.11 | 1.41 | **14** | **157** | — | — |

Reproduce this table, then check it against the emulation (T3), then fix it.
The limits of the current expansion, at 1 mm accuracy:

| layer | valid for | layer | valid for |
|---|---|---|---|
| L1 | pT > 0.26 GeV | L4 | pT > 0.88 GeV |
| L2 | pT > 0.43 GeV | L5 | pT > 1.15 GeV |
| L3 | pT > 0.64 GeV | L6 | **pT > 1.52 GeV** |

The repair is to extend the series: `x10_i = 1 + x²/6` → `1 + x²/6 + 3x⁴/40`
(one extra multiply and add per projection, reusing `x12A_i = x²`), and
`der_phiL_i = der_phiL·(1 + x²/2)` for the derivative. Both need their bit
widths rechecked — `x10`'s range grows from 1.1 to about 1.3 — and the HLS side
carries the same expansion, so it is a firmware change as well as an emulation
one.

*Deliverable:* the table above, measured; then the same table after extending
the expansion, showing what order is needed to reach the match window at
0.5 GeV.
*Learns:* series expansions in fixed point; that "the algorithm is exact" is
never true; how to compare an arithmetic error to a physical cut.
*Consequence to state plainly in the thesis:* below roughly 1.5 GeV the
dominant error on a projection to L6 is **not** multiple scattering or
detector resolution, it is the arithmetic. Widening the match windows cannot
compensate for a systematic offset an order of magnitude larger than the
window.

### T3. Does the emulation agree with itself?

Verify T2 against the emulation rather than only on paper. `MatchProcessor.cc`
computes three residuals per stub: `dphi` against the exact helix,
`dphiapprox` against the same arithmetic the firmware does, and `ideltaphi` in
integers. Only the last one decides the match.

`writeMonitorData("Residuals")` is already enabled and writes
`layerresiduals.txt` with columns 4 and 5 (integer and approximate), but not
`dphi`. Add one column at `MatchProcessor.cc:574`, rerun, and plot
`dphi − dphiapprox` versus pT per layer.

*Deliverable:* the measured truncation error, to compare with T2's prediction.
*Learns:* reading emulation internals; that a quantity can be "in the code"
three times with three different meanings.

---

## Phase 2 — the output format (months 2–3)

### T4. The 1.9 GeV wall in the track word

The 96-bit output track word encodes curvature with `minRinv = -0.006` in
`DataFormats/L1TrackTrigger/interface/TTTrack_TrackWord.h:87` — exactly
1.9 GeV — and `digitizeSignedValue` **clamps silently**. The tracks are found
correctly and the ntuple reads the floating-point momentum, so nothing in our
plots shows it; but every downstream L1 client (GTT, vertexing, jets, MET,
particle flow) decodes the word.

Steps:

1. Add `trk_rinv_word` and `trk_pt_word` branches to the ntuplizer, filled from
   `TTTrack::getRinvWord()` / `getRinv()`.
2. Run on a low-pT sample and plot word-pT against float-pT. Expect every track
   below 1.9 GeV to pile up in one bin.
3. Change `minRinv` to −0.024 in a local checkout, rerun, and show the pile-up
   is gone. Quantify what it costs: the curvature LSB coarsens 4×, so δpT/pT at
   2 GeV goes from 0.006% to 0.026% — far below the fit resolution.
4. Look one step further downstream: GTT re-encodes pT as `ap_ufixed<10,7>`
   ([`L1GTTInputProducer.cc:66`](../L1Trigger)) — LSB 0.125 GeV, so 0.5 GeV is
   four counts, 25% granularity. Work out whether reallocating integer bits to
   fractional bits at constant width would fix it.

*Deliverable:* one figure and a one-constant recommendation. This is the
evidence for a request to the L1 upgrade group, and it has the longest lead
time of anything in the project because it needs external agreement.
*Learns:* CMSSW plugin editing, DataFormats, and that a correct algorithm can
still be unusable if the output format cannot express it.

---

## Phase 3 — the thesis core (months 3–6)

### T5. Match windows from the residuals

**The freest knob in the system, and nobody has touched it.** The cuts
`rphimatchcut_` and `zmatchcut_` in `Settings.h` are still the values tuned for
2 GeV tracks, while multiple scattering scales as 1/pT — a 0.5 GeV track's
residual is roughly four times larger against an unchanged window.

Changing them costs **nothing**: no bits, no memory, no firmware interface, so
long as the LUT entries stay under the 10-bit ceiling of 1023
(`TrackletLUT::initmatchcut`). With full duplication the largest entry drops
from 741 to about 565, leaving ~1.8× flat headroom, more per layer.

Do not scan blindly. **Measure the residual distribution and set the cut from
it.** The machinery exists: `MatchProcessor.cc` computes the r-φ and z residuals
for every candidate and `FillLayerResidual` already receives a `truthmatch`
flag, so the distribution for genuine matches can be separated from
combinatorics.

1. **Open the windows wide first** — 4–5× — and only then measure. This is the
   trap in this task: with the current cuts you can only observe residuals that
   already passed them, so the measured width is biased low by construction and
   the answer will look fine when it is not.
2. Plot the r-φ and z residual for true matches versus pT, per layer and per
   seed. Fit the width.
3. The r-φ width should go as 1/pT (scattering). If instead it shows a
   pT-dependent **shift**, that is T2's arithmetic error and must be fixed
   before any tuning is meaningful.
4. Set each cut to cover the measured distribution — the current values
   correspond to some number of sigma at 2 GeV; keep that convention and let
   the widths scale.
5. Check the 1023 ceiling. L1 is the layer near it; L2–L6 have 2–3× more room,
   so a per-layer result is the natural outcome, not a flat factor.
6. Measure what it costs: fake rate and duplicate rate versus pT, and the
   occupancy effect — wider windows mean more candidates per match engine,
   which feeds straight into T7.
7. The physically right answer is a window that depends on pT, since that is how
   the scattering scales. There is precedent for an rinv-binned table at
   `TrackletLUT.cc:928`. This is the ambitious version and would be a strong
   thesis result.

z deserves its own treatment: the z residual is driven by the z0 spread and the
tilted-module geometry rather than by scattering alone, so it will not scale the
same way as r-φ.

*Blocked by:* low-pT samples, frozen configuration, and T2 — the cut cannot be
set from a residual distribution that is dominated by an arithmetic offset.
*Learns:* setting a selection from a measured distribution rather than by
scanning; selection-bias traps; that "more efficiency" is never free.

### T5b. The `reachesRadius` threshold

Nicole's guard is

```cpp
bool reachesRadius(double r, double rinv) { return 0.5 * r * std::abs(rinv) < 0.9; }
```

The geometric limit is 1.0 — a track whose turning diameter just reaches the
layer, arriving tangentially. The 0.9 is a tunable that nobody has scanned, and
it does two different jobs that should be separated:

- **preventing NaN and overflow.** `asin` is undefined past 1, and the exact
  derivative carries 1/√(1−x²), which diverges. This is what the guard was
  written for.
- **rejecting tracks whose projection is worthless.** A tangential crossing
  gives a poor stub and an ill-conditioned derivative even when the arithmetic
  is finite.

It does **not** protect the arithmetic: at x = 0.9 the truncated expansion is
already off by about 10 cm at every layer (T2). So today the guard is admitting
tracks whose projections are meaningless.

1. Work out the per-layer pT floor the threshold implies — at 0.9, L6 needs
   0.685 GeV where the geometric limit is 0.617.
2. Scan it — 0.99, 0.95, 0.9, 0.8, 0.7 — against efficiency, fake rate,
   duplicate rate, and numerical health (how many projections hit `atExtreme`,
   how many derivative words saturate).
3. Report the interaction with T2: with an extended expansion the threshold can
   be raised toward the geometric limit, recovering real tracks. Without it, the
   threshold is doing the expansion's job badly.

*Blocked by:* low-pT samples. Cheap once unblocked — a one-line change and a
scan.

### T6. Does the low-pT algorithm need wider memories at all?

Nicole widened three internal words: `nbitsrinv_` 14→17, `nbitsrinvfit_` 15→16,
`nbitsphiprojder*` 10→12. Every fixed-point field satisfies
`range = 2^bits × LSB`, so the same range is reachable by coarsening the LSB
instead — and the LSB route needs no memory change at all.

The arithmetic says a constants-only configuration exists:

| Nicole | constants-only alternative | cost |
|---|---|---|
| `nbitsrinv_` 17 | keep 14, `rinv_shift_` −8 → −6 | δpT/pT 0.10% at 2 GeV |
| `nbitsrinvfit_` 16 | keep 15, follows automatically | none |
| `nbitsphiprojder*` 12 | keep 10, `SS_phiderL_shift_` −5 → −3 | projection quantisation at L6: 0.09 → 0.18 mm |

Run both configurations on the same events and compare efficiency, fake rate
and resolution. If they are indistinguishable — which the numbers suggest — the
result is that **the low-pT algorithm requires no rewidening of any memory**,
only changed constants. That is a much easier change to get into firmware, and
it is a genuinely useful conclusion for the note.

Note the coupling: coarsening the projection derivative eats into the match
window, which T5 is simultaneously trying to widen. The two tasks have to be
done in that order and the interaction reported.

*Blocked by:* frozen configuration. Four numbers in `Settings.h`, so cheap once
unblocked.

---

## Phase 4 — system limits (months 6–8)

### T7. Occupancy, truncation, and whether this survives PbPb

The question a trigger reviewer asks first, and the one nobody has answered.
Each processing module gets ~108 clock cycles (`maxstep_`, fixed by the L1
latency budget) and each virtual-module bin holds at most 15 stubs
(`maxStubsPerBin_`, "16 causes overflow"). Both were sized for pp with 200
pileup. Central HYDJET has far more stubs per module, and lowering `ptcutte_`
to 0.3 multiplies the candidates each engine must walk through.

There is a second effect: full duplication sets `regions_.set(0); regions_.set(1)`
unconditionally in `TrackerDTC/src/Stub.cc`, so **every stub now goes to both
track-finding processors** where the boundary logic previously sent most to one.
That is close to a doubling of the input rate per processor.

1. Add counters for modules that hit `maxstep_` and bins that overflow
   `maxStubsPerBin_`. They exist as internal state but are not written out.
2. **Virtual-module occupancy as a function of multiplicity**, separately for
   pPb and PbPb — they differ by more than an order of magnitude in N_ch
   (~250 versus ~8000), so they are two different regimes and must be plotted
   on their own axes, not overlaid. For each: the mean and the tail of the
   stub count per VM bin against the limit of 15, and the fraction of bins
   that overflow.
3. Break the loss down by layer and by seed. L1 has the highest stub density
   and the most virtual modules (`nbitsallstubs_` = 3 there against 2
   elsewhere), so it will saturate first.
4. Plot the fraction of events affected versus multiplicity, and efficiency
   versus multiplicity alongside it, to show whether the efficiency loss at
   high occupancy is truncation or something else.
5. Separate the two causes of the increase: the lower `ptcutte_`, which
   multiplies the candidates each engine walks through, and full duplication,
   which roughly doubles the stubs arriving at each processor. Running with one
   and not the other separates them.
6. Measure the stub-rate increase from full duplication; the DTC's own `lost_`
   counters cover the input side.
7. If the limits are exceeded, the response is not to raise them — `maxstep_`
   is fixed by the L1 latency budget and `maxStubsPerBin_` by memory. The
   response is more units in parallel (`teunits_`), which costs resources. State
   what it would cost.

*Deliverable:* figures F6 and F10 of the DP note list. This is the difference
between "works in simulation" and "works in the trigger".
*Learns:* that an FPGA algorithm has a time budget, and that exceeding it fails
silently as dropped data rather than as an error.

---

## Lower priority — good extensions, not the spine

**T8. Bremsstrahlung and electrons.** L1 efficiency for electrons is 0.52 at
2–2.5 GeV against 0.96 for muons; offline reaches 0.91 against 0.99. The stubs
are there — 95.5% of electrons have ≥ 4 — so it is the single-helix assumption
failing, not missing hits. Understanding whether the low-pT changes recover any
of it is a real question, but it opens a second front before the main result is
established. `plot_eff_species.C` already does the comparison.

**T9. Track parameter resolution below 2 GeV.** σ(pT)/pT, σ(z0), σ(φ0) versus
pT extended into the new range, to confirm the fixed-point granularity is
adequate. `compare_res.C` exists; needs the low-pT samples.

**T10. Sample production.** Mechanical but necessary, and learning to submit
CRAB jobs is worth it. Better shared than owned.

**T11. Firmware.** Only if the student is strong and time remains. The HLS
project needs Vivado 2020 and the test-vector recipe in `firmware/`.

---

## Ordering rationale

| | task | value | blocked? | difficulty |
|---|---|---|---|---|
| 1 | T0 reproduce | onboarding | no | low |
| 2 | T1 layer reach | note figure | no | low |
| 3 | **T2 projection error** | **highest — decides whether T5 is tuning or repair** | **no** | medium |
| 4 | T3 residual cross-check | validates T2 | partly | medium |
| 5 | T4 track word | high — longest external lead time | partly | medium |
| 6 | **T5 match windows from residuals** | **high — free and untouched** | yes, incl. T2 | medium |
| 7 | T5b `reachesRadius` threshold | medium — one line, unscanned | yes | low |
| 8 | T6 constants-only | high — firmware adoption | yes | low |
| 9 | **T7 VM occupancy vs multiplicity** | **highest for FPGA credibility** | yes | high |

T2 sits third on purpose, and the numbers in it are the reason. At 0.5 GeV the
projection into L4 is wrong by 25 mm against a window of 1–3 mm, and the
`reachesRadius` guard at 0.9 admits tracks whose projections are off by 10 cm.
Until the expansion is extended, T5 is measuring an arithmetic offset rather
than a physical residual, and T5b is compensating for it with a blunt cut.
Better to establish that in month one than in month five.

The three tasks are also coupled and should be reported together: the
expansion order (T2) sets how tight the windows can be (T5), the windows set
how many candidates each engine must walk (T7), and the `reachesRadius`
threshold (T5b) trades all three against efficiency.
