# Projection arithmetic

What the tracklet emulation approximates when it predicts where a track crosses
a layer, and why those approximations fail at low pT. No data needed — this is
the helix and the arithmetic the code actually performs.

```
cd scripts
root -l -b -q 'plot_projection_error.C()'
```

---

## 1. The exact helix

A prompt track in a uniform field **B** along z is described by four parameters:
the signed curvature `rinv = 1/R`, the azimuth at the point of closest approach
`phi0`, the longitudinal impact `z0`, and `t = pz/pT = sinh(eta)`.

```
rinv = 0.01 * c * B / pT  =  0.0114 / pT        [cm^-1, pT in GeV, B = 3.8 T]
```

Define the dimensionless

```
x  =  r * rinv / 2
```

`x` has two meanings at once, and that is the root of everything below:

- **geometric** — the track reaches radius `r` only if `x <= 1`. `x = 1` is a
  tangential crossing at the turning diameter.
- **angular** — `asin(x)` is the azimuth the track has swept between the beam
  line and radius `r`.

### Barrel: given the radius r, find phi and z

```
phi(r) =  phi0 - asin(x)
z(r)   =  z0 + (2t/rinv) * asin(x)

dphi/dr = -(rinv/2) / sqrt(1 - x^2)
dz/dr   =       t   / sqrt(1 - x^2)
```

These are implemented exactly in `TrackletCalculatorBase::exactproj`
(`src/TrackletCalculatorBase.cc:186-188`) and used only for validation.

### Disks: given z, find r and phi

```
u      =  rinv * (z - z0) / (2t)          [ u = asin(x), the swept angle ]

phi(z) =  phi0 - u
r(z)   =  (2/rinv) * sin(u)

dphi/dz = -rinv / (2t)
dr/dz   =  cos(u) / t
```

`TrackletCalculatorBase::exactprojdisk`. **Note that `phi(z)` is linear in z** —
there is no arcsine in the disk azimuth at all.

---

## 2. What the emulation computes instead

An FPGA has no `asin`, `sin` or `sqrt`, so `IMATH_TrackletCalculator.h` evaluates
truncated series.

### Barrel (lines 245-299)

```
asin(x)   ->   x * (1 + x^2/6)                       [Taylor, O(x^3)]

phi(r)    ->   phi0 - x * (1 + x^2/6)                 phiL_i   (line 260)
z(r)      ->   z0 + t*r * (1 + x^2/6)                 zL_i     (line 285)

dphi/dr   ->   -rinv/2                                der_phiL (line 266)
dz/dr     ->    t                                     der_zL   (line 298)
```

The identity `z0 + (2t/rinv)*(x + x^3/6) = z0 + t*r*(1 + x^2/6)` is why the same
`x10_i = 1 + x^2/6` factor appears in both the phi and z chains.

**The two derivatives drop `1/sqrt(1-x^2)` entirely.**

### Disks (lines 325-381)

```
sin(u)    ->   u * (1 - u^2/6)                       [Taylor, O(u^3)]

phi(z)    ->   phi0 - u                               phiD_i   -- EXACT
r(z)      ->   ((z-z0)/t) * (1 - u^2/6)               rD_i     (line 368)

dphi/dz   ->   -rinv/(2t)                             der_phiD -- EXACT
dr/dz     ->    1/t                                   der_rD   (line 380)
```

**Disk azimuth and its derivative are exact.** Only the radius is approximated,
and `dr/dz` drops the `cos(u)`.

### The match

`MatchProcessor::matchCalculator` (`src/MatchProcessor.cc:532-538`) evaluates
the projection at the layer's **nominal** radius `rmean` and then takes one
linear step to the stub's own radius:

```
dr    =  r_stub - rmean                               |dr| <= drmax = 3.75 cm

dphi  =  phi_stub - [ phi(rmean) + dr * (dphi/dr) ]
dz    =  z_stub   - [ z(rmean)   + dr * (dz/dr)   ]
```

A stub matches if `|dphi| * rmean < rphimatchcut_` and `|dz| < zmatchcut_`.

---

### What the predicted trajectory actually is

The emulation does not fit a curve to the layer. It holds one number constant:

```
dphi/dr = -rinv/2                  ->   phi(r) = phi_anchor - (rinv/2)*(r - rmean)
```

with `rinv` taken from the two seed stubs. `phi` linear in `r` is a straight
line in the coordinates the algorithm works in, but in the transverse plane it
traces a spiral, and it is less curved than the track:

| pT | true helix radius | radius of the assumed trajectory | ratio |
|---|---|---|---|
| 2.0 GeV | 175.5 cm | 182.2 cm | 1.04 |
| 0.6 GeV | 52.7 cm | **73.9 cm** | **1.40** |

The ratio is `sec(theta)` to first order -- the same factor as in A2 below. So
the assumption is not "the track is straight", it is "the track curves at the
rate it would if it were still travelling radially".

## 3. The three approximations, and their errors

Write the position error on the layer as the angular error times the radius.

### A1 — the series truncation

```
asin(x) - x*(1 + x^2/6)  =  3x^5/40 + 15x^7/336 + 105x^9/3456 + ...

eps_1  =  r * ( 3x^5/40 + 15x^7/336 + ... )
```

Independent of `dr`: it applies to every stub in the layer. The series converges
slowly near `x = 1` because `asin` has a square-root singularity there — its
derivative diverges, and no polynomial has a vertical tangent.

### A2 — the dropped slope factor

```
1/sqrt(1-x^2)  =  1 + x^2/2 + 3x^4/8 + ...

eps_2  =  r * |dr| * (rinv/2) * ( 1/sqrt(1-x^2) - 1 )
```

Grows linearly with `dr`, so it is worst for stubs at the edge of a layer. This
is the same arcsine problem seen through its derivative.

### A3 — the linear step

```
d^2phi/dr^2  =  -(rinv/2)^2 * x / (1-x^2)^(3/2)

eps_3  =  r * (dr^2/2) * (rinv/2)^2 * x / (1-x^2)^(3/2)
```

Second order in `dr`, and the smallest of the three.

### Sizes at L4 (r = 68.7 cm), |dr| = drmax, window 1.9 mm

| pT [GeV] | x | eps_1 (series) | eps_2 (slope) | eps_3 (linear) |
|---|---|---|---|---|
| 2.0 | 0.196 | 0.02 mm | 0.14 mm | 0.01 mm |
| 1.0 | 0.391 | 0.52 mm | 1.27 mm | 0.08 mm |
| 0.7 | 0.559 | 3.49 mm | 4.32 mm | 0.31 mm |
| 0.5 | 0.783 | 25.0 mm | 17.8 mm | 2.04 mm |

**A1 and A2 are comparable and both dominate A3 by an order of magnitude.**
Nicole's TODO addresses the linear step (A3 and the slope inside A2); the
series (A1) is untouched and equally large.

---

## 4. What a fix has to do

The series cannot simply be extended: for the L3 projection to stay inside its
1.0 mm window needs 5th order at 0.50 GeV, 7th at 0.45, 9th at 0.40, and more
than 15th below 0.35 — the Taylor series is the wrong tool near `x = 1`.

The practical options:

1. **Table the correction.** `D(x) = asin(x) - x - x^3/6` depends only on `x`,
   which is already computed in the chain as `x8_i`. A 256-entry interpolated
   ROM gives better than 0.1 mm at L3. `imath` already contains a LUT-valued
   variable, `VarInv` (`interface/imath.h:956`), whose `print(..., HLS)` emits
   the ROM straight into the firmware — the same pattern applies here.
2. **Range reduction.** `asin(x) = pi/2 - 2*asin(sqrt((1-x)/2))` maps `x = 0.9`
   onto an argument of 0.224, where the existing cubic is accurate to 1e-4 rad.
   Costs a square root, itself a table.

The slope needs the same treatment: `dphi/dr = -(rinv/2) * g(x)` with
`g(x) = 1/sqrt(1-x^2)` from a second small table. Note `g` grows to 2.29 at
`x = 0.9`, so the derivative word range has to be rechecked —
`nbitsphiprojderL123_` / `L456_` and `SS_phiderL_shift_`.

**The disks need neither.** Azimuth is already exact there; only `r(z)` is
approximated, and `sin` converges fast (next coefficient 1/120 against the
arcsine's 3/40, no singularity). Two more terms take the D1 radius error at
0.33 GeV from 11.7 cm to 0.32 cm at eta = 1.0, inside the 0.5 cm `rcutPS_`.

---

## 5. Figures

Three standalone figures, each PDF and PNG, all from `scripts/plot_projection_error.C`:

| file | what it shows |
|---|---|
| `projection_geometry` | the transverse plane for a 0.5 GeV track: the true helix, the locus the cubic predicts, the layer crossings of each, and the turning point at `x = 1`. A 2 GeV track is overlaid, where the two curves coincide. This is the figure that defines `x` geometrically: `x = r/2R`. |
| `projection_error_vs_pt` | the position error versus pT for L3-L6, against the real `rphimatchcut_` windows. Where a curve crosses its own dotted line, that layer stops matching. |
| `projection_error_in_layer` | two panels, two tracks at different azimuth: one crossing the inner ladder of L4, one crossing the outer ladder. Each panel is drawn in the frame of its own track, so that track is the line y = 0 with its stub on it, and the shaded band is the match window around where the algorithm looks. The ladders and sensor planes are at their measured radii. At 0.80 GeV the inner-ladder stub is 0.4 mm from the prediction and matches; the outer-ladder one is 2.6 mm away and is lost. |
| `projection_error_in_layer_0p6` | the same at 0.60 GeV, where the series offset alone is 8.3 mm. Both ladders are lost -- 4.5 mm on the inner, 11.0 mm on the outer -- so the layer contributes nothing at all. The vertical range scales with pT; the layout is identical. |
| `projection_transverse`, `projection_transverse_0p8` | the canonical transverse view, zoomed onto one patch of L4 and drawn to scale, at 0.60 and 0.80 GeV. Radius runs upward, azimuth across. The ladders alternate between two radii around phi, so at a given azimuth there is either an inner-ladder module or an outer-ladder one; track A crosses the inner, track B the outer. For each: the true helix, the straight line the algorithm predicts, the stub, and the match window around the prediction at that radius. The square marker shows where the predicted line is anchored -- at `rmean`, a radius with no module on it, already displaced there by the series truncation -- and from that point it is extended with the small-angle slope. The left panel is the same two tracks from the beam line, with the patch marked. L4 is drawn there as what it is -- flat modules alternating between two radii all the way round in phi. The two circled markers are the seed stubs on L1 and L2, which fix `rinv` and `phi0`; the squares are where that seed then projects, layer by layer, and the dashed curve is the locus of those projections, `phi0 - (x + x^3/6)`. Stubs are placed on their real ladder radii, never at `rmean`, and the projection is stepped linearly from `rmean` to the stub's own radius exactly as `MatchProcessor` does. At 0.6 GeV the resulting gap is 0.01 cm at L1, 0.04 at L2, 0.03 at L3, 0.47 at L4 and 3.10 at L5, against windows of a millimetre or two. The dotted radial line and the angle `theta` at the crossing complete the statement: `-rinv/2` is the slope of a track crossing radially, and the real track is at `theta = asin(x)`, so its slope is `sec(theta)` steeper. |

| `projection_slope` | where `1/sqrt(1-x^2)` comes from. `asin(x)` is the angle between the track and the radial direction where it crosses the layer, so the factor is `sec(theta)`. At 2 GeV the crossing angle at L4 is 11 deg and the secant is 1.02; at 0.5 GeV it is 52 deg and the secant is 1.61, so the direction the code assumes is 60% too shallow and the prediction drifts 20 mm across the layer, against a 1.9 mm window. |

```
root -l -b -q 'plot_projection_error.C()'     # all four
root -l -b -q 'plot_projection_error.C(1)'    # just the geometry
```

**The real geometry of a barrel layer**, measured from the stub and cluster
positions in an EPOS pPb ntuple (Extended2026D110), not assumed:

| layer | nominal `rmean` | ladder radii [cm] | full spread [cm] |
|---|---|---|---|
| L1 (PS) | 24.9 | 21.89, 24.27 | 21.24 - 29.01 |
| L2 (PS) | 37.2 | 34.76, 37.13 | 34.20 - 40.88 |
| L3 (PS) | 52.3 | 49.93, 52.31 | 49.40 - 55.88 |
| L4 (2S) | 68.7 | 66.99, 70.15 | 66.63 - 70.92 |
| L5 (2S) | 86.0 | 84.29, 87.45 | 83.94 - 88.19 |
| L6 (2S) | 108.3 | 106.58, 109.74 | 106.23 - 110.46 |

Resolving the clusters rather than the stubs shows the two-sensor structure of
each module: in L4 the sensor planes sit at 66.80, 67.25, 69.97 and 70.41 cm,
i.e. modules of two sensors about 4.5 mm apart on ladders about 3 cm apart.

**Note that `rmean` falls in the gap between the two ladders**, so no stub is
ever at the radius where the projection is evaluated; every stub needs the
linear step, by roughly +-1.7 cm. That step is exactly linear: the projection
value and its derivative are evaluated once at `rmean`, and
`MatchProcessor.cc:532` forms `phi(rmean) + dr * dphi/dr` with
`dr = r_stub - rmean`. Nothing is re-evaluated at the stub's own radius.

**What "stub radius" means.** A barrel layer is not a mathematical cylinder.
Its modules sit on ladders staggered between two radii so they overlap in phi
without gaps, and toward the ends of the barrel they are tilted to point at the
interaction region. A stub is formed from two clusters in the two sensors of one
module, so its radius is that module's actual radial position, which differs
from the layer's nominal `rmean` by up to `drmax = 3.75 cm`. The emulation
stores exactly that offset: the stub word carries `r` as a 7-bit signed number
with `krbarrel = 2*drmax/2^7 = 0.059 cm` per count. The projection is evaluated
once, at `rmean`, so every stub needs the linear step to its own radius -- which
is why the slope matters at all.
