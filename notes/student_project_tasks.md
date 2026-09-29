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

**The most valuable thing a new student can do here, and it needs no data.**

The emulation projects a track to a layer with a third-order expansion of the
arcsine ([`IMATH_TrackletCalculator.h:250-258`](../L1Trigger/TrackFindingTracklet)):

```
phi(r) = phi0 - x·(1 + x²/6)          where x = r·rinv/2      [exact: phi0 - asin(x)]
dphi/dr = -rinv/2                                             [exact: -rinv/2 / sqrt(1-x²)]
```

Both are excellent at 2 GeV and both degrade as x → 1. Quantify the error in
**millimetres at the layer** (multiply the angular error by the layer radius),
as a function of pT and layer, and compare it to the match window
`rphimatchcut_` for that layer and seed.

Then extend: add the fifth-order term 3x⁵/40 and the 1/√(1−x²) factor to the
same expressions and show how much of the error is recovered.

*Deliverable:* a two-panel figure (error vs pT per layer, with the match window
overlaid) and a table. Goes straight into the internal note.
*Learns:* series expansions in fixed point, why "the algorithm is exact" is
never true, how to compare an arithmetic error to a physical cut.
*Watch for:* the answer should be a few centimetres at L4 for 0.5 GeV against
windows of 1–3 mm. If it is, that is a result: **no widening of the match
windows can compensate for it**, and the expansion has to be extended.

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

### T5. Match-window optimisation

**The freest knob in the whole system, and nobody has touched it.** The match
cuts `rphimatchcut_` / `zmatchcut_` in `Settings.h` are still the values tuned
for 2 GeV tracks, while multiple scattering scales as 1/pT — a 0.5 GeV track's
residual is roughly four times larger against an unchanged window.

Changing them costs **nothing**: no bits, no memory, no firmware interface, as
long as the LUT entries stay under the 10-bit ceiling of 1023
(`TrackletLUT::initmatchcut`). With full duplication the largest entry drops
from 741 to ~565, leaving about 1.8× flat headroom, more per layer.

1. Add a configurable scale factor rather than editing the array by hand.
2. Scan it and measure efficiency, fake rate and duplicate rate versus pT.
3. Find where the efficiency gain stops paying for the fake-rate cost.
4. Then stop scaling flatly: L1 is the layer near the ceiling, L2–L6 have 2–3×
   more room, so a per-layer scan is the natural next step.
5. The physically right answer is a window that depends on pT, since that is how
   the scattering scales. There is precedent for an rinv-binned table at
   `TrackletLUT.cc:928`. This is the ambitious version and would be a strong
   thesis result.

*Blocked by:* low-pT samples, frozen configuration.
*Learns:* parameter optimisation against a figure of merit; that "more
efficiency" is never free.

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
2. Plot the fraction of events affected versus charged multiplicity, for EPOS
   pPb and HYDJET PbPb.
3. Plot efficiency versus multiplicity alongside it, to show whether the
   efficiency loss at high occupancy is truncation or something else.
4. Measure the stub-rate increase from full duplication; the DTC's own `lost_`
   counters cover the input side.

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
| 3 | **T2 projection error** | **high — may invalidate the window approach** | **no** | medium |
| 4 | T4 track word | high — external lead time | partly | medium |
| 5 | T3 residual check | validates T2 | partly | medium |
| 6 | **T5 match windows** | **highest — free and untouched** | yes | medium |
| 7 | T6 constants-only | high — firmware adoption | yes | low |
| 8 | T7 occupancy | highest for FPGA credibility | yes | high |

T2 is placed early on purpose. If the projection error at L4 really is
centimetres against millimetre windows, then T5 cannot succeed on its own and
the expansion has to be fixed first — better to know that in month one than in
month five.
