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
| `projection_error_in_layer` | inside L4 at 0.7 GeV. The shaded band is the range of radii a stub can have in one layer, `\|dr\| <= drmax = 3.75 cm`; the dotted vertical line is `rmean`, the single radius at which the projection is evaluated. Two arrows mark the two errors: at `rmean` the blue curve is zero by construction, so the gap to the red one there is the series truncation; away from `rmean` both rise together, and that common tilt is `slope x dr`. Dotted horizontals are the match window. |

```
root -l -b -q 'plot_projection_error.C()'     # all three
root -l -b -q 'plot_projection_error.C(1)'    # just the geometry
```
