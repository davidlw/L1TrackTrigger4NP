# DP note — figure list

Working document for a CMS Detector Performance note on the **low-pT L1 track
trigger for heavy-ion and UPC collisions**. One entry per candidate figure:
what it shows, the sample and selection, the macro that makes it, a draft
caption, and what is still missing.

Target: Quark Matter 2027 (21–27 Mar 2027, Jeju). Abstracts historically close
3–4 months ahead, so the approval chain has to start in **October 2026**.

Status legend — **R** ready (figure exists, numbers final), **S** macro ready,
sample missing, **A** waiting on the new algorithm, **N** not written yet.

---

## 1. The argument

Each figure answers the next objection, in order:

1. The standard Phase-2 track finder stops at ~2 GeV — that is a design choice,
   not a limit of the detector. → **F1, F2**
2. A low-pT version reaches 0.5 GeV with efficiency approaching offline. → **F3, F4**
3. It does not fall apart in central PbPb. → **F5, F6**
4. The purity cost is bounded. → **F7, F8**
5. It fits in the firmware: data volume, truncation, resources. → **F9, F10, F11**
6. It buys real physics. → **F12, F13**

A first note should carry **F3, F5, F6, F7, F12, F13 and one of F9/F10** as the
approved set; the rest are supporting material for the internal note, or a
second DP note once the FPGA work lands.

---

## 2. Conventions

Apply to every figure, so they read as one set.

| item | value |
|---|---|
| label, upper left | `CMS Phase-2 Simulation Preliminary` (check against current CMS guidance before approval) |
| label, upper right | collision system and energy, e.g. `PbPb 5.36 TeV`, `pPb 8.8 TeV`, `PbPb 5.36 TeV (UPC)` |
| pileup | state it even when zero — these are HI samples, `PU 0` |
| acceptance | `|eta| < 2.4` everywhere; never plot outside it |
| truth | **primary** charged particles only (`tp_ngenpart > 0`); say so in every caption |
| default algorithm | "standard Phase-2 L1 track finding (tracklet + KF)" — never "old" |
| new algorithm | "low-pT L1 track finding" — never "dummy stub", which is jargon |
| offline | "offline tracking (`generalTracks`)", with the quality requirement named |
| error bars | Clopper–Pearson, 68.3% |
| x axis | log for pT whenever the range spans 0.3–10 GeV |

**Definitions to repeat verbatim in captions:**

```
efficiency     = N(truth particle matched by >= 1 track) / N(truth particle)
fake rate      = N(track not matched to a truth particle) / N(tracks)
duplicate rate = N(truth particle matched by > 1 track) / N(matched truth particles)
```

**Sample table** (to appear once, in the note text):

| sample | generator | what it is for |
|---|---|---|
| γγ → μ⁺μ⁻ | STARlight | clean two-track UPC benchmark |
| γγ → e⁺e⁻ | STARlight | same, with bremsstrahlung |
| γA → J/ψ → μ⁺μ⁻ | STARlight | photonuclear reference |
| pPb min. bias | EPOS | intermediate occupancy |
| PbPb min. bias | HYDJET | maximum occupancy |

---

## 3. Tier 1 — the core claim

### F1. Truth pT spectrum and the fraction below threshold  — **S**

*What it shows:* the charged-particle spectrum in UPC and MB heavy-ion events,
with the fraction of particles below 2 GeV marked. Establishes the motivation
before any algorithm is discussed.

- **x** truth pT, 0.3–10 GeV, log · **y** dN/dpT per event, and a second panel
  with the cumulative fraction above pT
- **curves** γγ→μμ, γA→J/ψ, EPOS pPb, HYDJET PbPb
- **macro** `analyses/L1_tracking_performance/scripts/plot_tpspectrum.C`
- **needs** samples with `ptMinTP = 0.3` (see §6)

> **Draft caption.** Transverse-momentum spectrum of primary charged particles
> within |η| < 2.4 for ultraperipheral and minimum-bias heavy-ion collisions
> simulated with STARlight, EPOS and HYDJET. The shaded band marks pT < 2 GeV,
> the region inaccessible to the standard Phase-2 Level-1 track finding, which
> contains XX% (γγ→μμ), XX% (pPb) and XX% (PbPb) of all primary charged
> particles in the acceptance.

### F2. Layer reach vs pT  — **N**

*What it shows:* the radius a helix reaches as a function of pT, with the six
barrel layers marked — 0.14, 0.21, 0.30, 0.39, 0.49, 0.62 GeV for L1…L6. A
single analytic figure that explains why the low-pT regime is qualitatively
different: below 0.62 GeV a track cannot produce six stubs at all.

- **x** pT 0.1–3 GeV · **y** maximum radius reached [cm], with horizontal lines
  at the layer radii (24.9, 37.2, 52.3, 68.7, 86.0, 108.3 cm)
- **macro** to write — a few lines, no sample needed

> **Draft caption.** Maximum radial distance reached by a charged particle in
> the 3.8 T field as a function of transverse momentum. Horizontal lines mark
> the mean radii of the six Phase-2 outer-tracker barrel layers. Below
> 0.62 GeV a track originating at the beam line cannot reach the outermost
> layer, and below 0.30 GeV it reaches at most the third; the number of stubs
> available to the trigger is therefore a function of pT.

### F3. Tracking efficiency vs pT, UPC  — **A**

**The headline figure.** Three curves on one pad, log x from 0.3 GeV.

- **sample** STARlight γγ → μ⁺μ⁻
- **selection** primary truth muons, |η| < 2.4
- **curves** standard L1 · low-pT L1 · offline (highPurity)
- **macro** `compare_eff.C` → `plot_eff.C` (offline overlay already supported)
- **needs** the new algorithm run on this sample with `ptMinTP = 0.3`

> **Draft caption.** Level-1 track-finding efficiency as a function of
> transverse momentum for primary muons from γγ → μ⁺μ⁻ in ultraperipheral PbPb
> collisions, for the standard Phase-2 algorithm and for the low-pT
> configuration described in the text. Offline tracking is shown for reference.
> Efficiency is the fraction of primary truth particles within |η| < 2.4
> matched to at least one reconstructed track. Vertical bars are 68.3%
> Clopper–Pearson intervals.

### F4. Tracking efficiency vs η and φ at fixed pT thresholds  — **A**

Shows the gain is uniform in acceptance and exposes any sector-boundary
artefact from the nonant duplication — a change a reviewer will ask about.

- **panels** η and φ, for pT > 0.5, 1, 2 GeV
- **curves** standard vs low-pT, same sample as F3
- **macro** `compare_eff.C` → `plot_eff.C` (η/φ panels exist)

> **Draft caption.** Level-1 track-finding efficiency versus pseudorapidity
> (left) and azimuth (right) for primary muons above three transverse-momentum
> thresholds, for the standard and low-pT configurations. The φ distribution is
> shown over one processing nonant to expose any dependence on the sector
> boundary.

---

## 4. Tier 1 — occupancy, the real test

### F5. Tracking efficiency vs pT, high occupancy  — **A**

Same axes and curves as F3, on EPOS pPb and HYDJET PbPb. If the low-pT
efficiency degrades, **publish the degradation** — an honest curve is worth
more than a UPC-only claim and pre-empts the obvious question.

- **macro** `compare_eff.C` → `plot_eff.C`
- **note** HYDJET statistics are currently thin (1 800 events in the existing
  offline scan); a larger production is needed for a publishable curve

> **Draft caption.** Level-1 track-finding efficiency versus transverse
> momentum for primary charged particles in minimum-bias pPb (EPOS, left) and
> PbPb (HYDJET, right) collisions, for the standard and low-pT configurations,
> with offline tracking for reference.

### F6. Tracking efficiency vs charged multiplicity  — **A**

*The most important single plot for a trigger audience.* It is the direct test
of truncation: each processing module has ~108 clock cycles (`maxstep`) and
each virtual-module bin holds at most 15 stubs (`maxStubsPerBin`), so
efficiency must be shown as a function of event occupancy, not only averaged
over it.

- **x** N_ch (primary truth, pT > 0.4, |η| < 2.4) · **y** efficiency
- **curves** standard vs low-pT, at a fixed pT > 0.5 GeV; EPOS pPb and HYDJET
  PbPb on separate pads (very different N_ch ranges: ~250 vs ~8000)
- **macro** `compare_effvsnch.C`

> **Draft caption.** Level-1 track-finding efficiency for primary charged
> particles with pT > 0.5 GeV as a function of the primary charged-particle
> multiplicity of the event, for the standard and low-pT configurations. The
> multiplicity is counted for primary charged particles with pT > 0.4 GeV and
> |η| < 2.4. The decrease at high multiplicity reflects the finite processing
> time and memory depth of the pattern-recognition step.

---

## 5. Tier 1 — the cost

### F7. Fake rate vs pT and vs N_ch  — **A**

- **curves** standard vs low-pT; both UPC and MB samples
- **macro** `compare_fake.C` → `plot_fake.C`

> **Draft caption.** Fraction of Level-1 tracks not matched to a primary
> charged particle, as a function of track transverse momentum (left) and event
> multiplicity (right), for the standard and low-pT configurations in
> minimum-bias PbPb collisions.

### F8. Duplicate rate vs pT  — **A**

Relevant because the low-pT configuration widens the match windows and, with
full nonant duplication, allows the same track to be found in two sectors.

- **macro** `compare_eff.C` (duplicate rate is computed alongside efficiency)

> **Draft caption.** Fraction of primary charged particles matched by more than
> one Level-1 track, for the standard and low-pT configurations.

---

## 6. Tier 2 — firmware feasibility

### F9. Stub multiplicity per event  — **R**

Already made. Four curves per sample: accepted and rejected stubs, standard vs
low-pT stub definition. Quantifies the data volume the DTC links must carry —
the first hardware objection.

- **samples** EPOS pPb, HYDJET PbPb
- **macro** `compare_stubs.C` / `compare_stubmult.C` → `plot_stubs.C`
- **note** relabel "dummy" → "low-pT stub definition" before approval

> **Draft caption.** Distribution of the number of stubs per event accepted and
> rejected by the front-end pT selection, for the standard and low-pT stub
> definitions, in minimum-bias pPb (left) and PbPb (right) collisions. Removing
> the front-end bend requirement increases the accepted stub multiplicity by a
> factor XX, which sets the input bandwidth the track-finding system must
> absorb.

### F10. Truncation loss vs occupancy  — **N**

The companion to F6, showing *where* the loss occurs rather than only that it
occurs: fraction of events in which at least one module hits `maxstep`, or one
VM bin overflows `maxStubsPerBin`, as a function of N_ch.

- **needs** counters added to the emulation (they exist as internal state; not
  currently written out)

> **Draft caption.** Fraction of events in which the processing-time or
> memory-depth limit of at least one pattern-recognition module is reached, as
> a function of charged-particle multiplicity, for the low-pT configuration.

### F11. Resource utilisation and latency  — **N**

If a preliminary VCU118 or HLS synthesis number exists by March, a table (LUT,
FF, BRAM, DSP, latency, per module type, standard vs low-pT) turns this from a
simulation study into a proposal. **Decide early whether this is in scope** —
without it, expect the question in every talk.

---

## 7. Tier 2 — the physics case

### F12. Event trigger efficiency per physics process ("money plot")  — **A**

- **layout** upper panel: efficiency per signal for standard L1 (pT > 2 GeV),
  low-pT L1 (pT > 0.5 GeV) and offline; lower panel: low-pT/standard gain on a
  log scale, factor printed at each point
- **denominator** common to all three series: events with ≥ k primary charged
  particles with pT > 0.5 GeV, |η| < 2.4 — physics acceptance factored out
- **groups** UPC γγ (ee, μμ, ττ, hh) · UPC γA (ρ⁰, J/ψ, ψ(2S), D⁰, dijet) ·
  hadronic MB (pPb, PbPb)
- **macro** `analyses/trigger_summary/scripts/money_plot.C`
- **needs** per-signal samples; several do not exist yet (see §9)

> **Draft caption.** Event trigger efficiency for a range of ultraperipheral
> and hadronic heavy-ion processes, for a Level-1 track trigger requiring at
> least k tracks above the threshold accessible to each configuration —
> 2 GeV for the standard algorithm, 0.5 GeV for the low-pT one — with offline
> tracking shown for reference. The denominator is common to all three: events
> containing at least k primary charged particles with pT > 0.5 GeV and
> |η| < 2.4, so that differences reflect tracking and threshold acceptance
> rather than the production cross section. The lower panel gives the gain of
> the low-pT configuration over the standard one.

### F13. Trigger rate vs multiplicity threshold  — **A**

Shows the thing is deployable inside an L1 budget, not merely efficient.

- **x** N_ch threshold (reconstructed L1 tracks, pT > 0.4 GeV, |η| < 2.4) ·
  **y** rate [Hz], log
- **normalisation** 50 kHz PbPb, 3 MHz pPb minimum-bias interaction rate
- **curves** standard vs low-pT, with a horizontal line at the allotted L1
  bandwidth
- **macro** `plot_rate_offline.C` (offline version exists; L1 version needed)

> **Draft caption.** Level-1 trigger rate as a function of the required number
> of reconstructed tracks with pT > 0.4 GeV and |η| < 2.4, assuming minimum-bias
> interaction rates of 50 kHz for PbPb and 3 MHz for pPb, for the standard and
> low-pT configurations.

---

## 8. Tier 3 — supporting, probably not in the first note

### F14. Track parameter resolution vs pT  — **A**

σ(pT)/pT, σ(z0), σ(η), σ(φ0) versus pT, extended into the 0.3–2 GeV region.
Needed by anyone who would build a selection on these tracks, and it shows the
fixed-point granularity is adequate at low pT.
**macro** `compare_res.C` → `plot_res.C`

### F15. Electron vs muon efficiency  — **R** (L1 and offline both measured)

A known-limitation figure. Bremsstrahlung breaks the single-helix assumption:
L1 efficiency at 2–2.5 GeV is 0.52 for electrons against 0.96 for muons, and
offline reaches 0.91 against 0.99. The η dependence follows the tracker
material. Strong material for the internal note and a follow-up; a distraction
in a first DP note whose job is to establish the main result.
**macros** `plot_eff_species.C`, `plot_eff_electron_quality.C`

### F16. Track quality variables  — **R**

χ²/dof, bend-χ², N_stub for standard vs low-pT. Mostly useful to justify that
the quality cuts still separate signal from combinatorics once the windows are
widened. **macro** `compare_trkquality.C` → `plot_trkquality.C`

### F17. Projection residuals vs pT  — **N**

Internal note only. The φ residual against the exact helix, per layer, versus
pT — exposes the third-order `asin` truncation in the projection arithmetic
(~2.5 cm at L4 for 0.5 GeV, against match windows of 1–3 mm). Needs one added
column in the `layerresiduals.txt` dump at `MatchProcessor.cc:574`.

---

## 9. Blockers, in priority order

1. **Truth samples with `ptMinTP = 0.3`.** The existing standard-configuration
   productions were made with 1.0 GeV, so below 1 GeV the *denominator* does
   not exist and no efficiency can be measured. Every Tier-1 figure depends on
   this. `simulations/commands_cmsDriver` already sets 0.3; the samples have to
   be regenerated from step 2. **This is the critical path, not the plotting.**
2. **The low-pT configuration itself**, validated and frozen: reachability
   guard, stub definition, match windows, nonant duplication, rinv bit width.
   Every "A" figure waits on one configuration that everyone agrees to call
   final.
3. **HYDJET statistics.** 1 800 events is enough for a cross-check, not for a
   published efficiency curve at high N_ch.
4. **Per-signal samples for F12.** STARlight covers ℓℓ, ττ, ρ⁰, J/ψ, ψ(2S);
   γγ→hh needs SuperChic; inclusive D⁰ and γA dijets need PYTHIA
   photoproduction. Decide which columns are worth generating and drop the rest
   from the figure rather than leaving placeholders.
5. **Truncation counters** for F10, and **synthesis numbers** for F11.

---

## 10. Production checklist

For each approved figure, before it goes into the note:

- [ ] both configurations run on the *same* events, file for file
- [ ] truth threshold is 0.3 GeV and stated in the caption
- [ ] primaries-only, and said so
- [ ] statistics sufficient that no bin has a Clopper–Pearson interval wider
      than the effect being claimed
- [ ] axis ranges identical across figures that will be compared
- [ ] the same colour and marker for a given configuration in every figure
- [ ] PDF and PNG, and the ROOT file with the histograms kept alongside
- [ ] the macro and its arguments recorded in this repository

---

*Figures live in `analyses/*/figures/`; the macros that make them in
`analyses/*/scripts/`. Nothing in `output/` or `figures/` is tracked by git —
regenerate rather than copy.*
