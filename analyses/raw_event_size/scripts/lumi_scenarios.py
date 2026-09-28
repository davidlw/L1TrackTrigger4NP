#!/usr/bin/env python3
"""Peak collision rates at CMS for Pb-Pb and p-Pb in Run 4 / Run 5, and the p-Pb
luminosity gain from running at 8.54 instead of 5.36 TeV.

Inputs are from R. Bruce, "Projected heavy-ion performance in Run 4 and Run 5"
(LPC, 2026-07-20/09-28):
  - Pb-Pb: levelled plateau ~6.4e27 at ATLAS/CMS in 2024-26 fills (p. 20);
    unlevelled start-of-fill peaks vs beta* from the CTE simulation (p. 27, IP2,
    same beta* and colliding bunches as IP1/5); 1032 colliding bunches at IP5.
  - p-Pb: CTE initial conditions (p. 28): levelling targets 1e30 / 2e30 at IP1/5,
    IP2 levelled at 5e29, beam intensities, emittances, cross sections.
The unlevelled p-Pb potential and the fill-averaged luminosity are computed
here (simple model: Pb burn-off at all IPs, no IBS/radiation damping).
Usage: python3 lumi_scenarios.py
"""
import math

F_REV = 11245.0
SIG_HAD_PBPB = 7.8e-24         # cm^2, at 5.36 TeV
PBPB_COLL = 1032

PBPB = [  # (label, peak L [cm^-2 s^-1])
    ("Run 4, levelled (2026-like, beta*=0.5 m)", 6.4e27),
    ("Run 4, no levelling (beta*=0.5 m)", 13e27),
    ("Run 5, no levelling, beta*=0.4 m", 16e27),
    ("Run 5, no levelling, beta*=0.3 m", 20e27),
    ("Run 5, no levelling, beta*=0.25 m (less likely)", 24e27),
]

PPB = {  # per energy: gamma-dependent geometric emittances (Pb, p) [m], cross sections [cm^2]
    5.36: dict(eps=(1.33e-9, 3.3e-10), had=2.084e-24, tot=2.1587e-24),
    8.54: dict(eps=(8.34e-10, 2.1e-10), had=2.13e-24, tot=2.2092e-24),
}
N_P, N_PB, BETA, HALF_XING, SIGMA_Z = 3e10, 2.31e8, 0.5, 75e-6, 0.098


def ppb_potential(E, ncoll):
    eps_pb, eps_p = PPB[E]["eps"]
    s2 = BETA * (eps_pb + eps_p)                       # Sigma_x^2 = sigma_1^2 + sigma_2^2 [m^2]
    phi = HALF_XING * SIGMA_Z / math.sqrt(s2 / 2)      # Piwinski angle
    return F_REV * ncoll * N_P * N_PB / (2 * math.pi * s2 * 1e4) / math.sqrt(1 + phi ** 2)


def ppb_fill_average(E, level, ncoll=1056, n_pb_bunches=1240, turnaround=3 * 3600, dt=60):
    """Best fill-averaged IP5 luminosity (optimal fill length), Pb burn-off only."""
    L0, N0 = ppb_potential(E, ncoll), n_pb_bunches * N_PB
    N, t, integ, best = N0, 0.0, 0.0, 0.0
    while t < 30 * 3600:
        x = N / N0
        l15, l2, l8 = min(L0 * x, level), min(L0 * x, 5e29), min(0.5 * L0 * x, 5e29)
        N -= PPB[E]["tot"] * (2 * l15 + l2 + l8) * dt
        integ += l15 * dt
        t += dt
        best = max(best, integ / (t + turnaround))
    return best


def mu(rate, ncoll):
    return rate / (ncoll * F_REV)


print("Pb-Pb at CMS (sigma_had = 7.8 b, %d colliding bunches)" % PBPB_COLL)
for label, L in PBPB:
    r = L * SIG_HAD_PBPB
    print(f"  {label:50s} L={L:.1e}  rate={r / 1e3:5.0f} kHz  mu={mu(r, PBPB_COLL):.3f}")

print("\np-Pb at CMS")
for E in (5.36, 8.54):
    for ncoll, scheme in ((1056, "25 ns p"), (955, "50 ns p")):
        pot = ppb_potential(E, ncoll)
        for label, L in (("levelled 1e30", 1e30), ("levelled 2e30", 2e30), ("no levelling", pot)):
            L = min(L, pot)
            r = L * PPB[E]["had"]
            nu = mu(r, ncoll) / (1 - math.exp(-mu(r, ncoll)))
            print(f"  {E} TeV {scheme}  {label:14s} L={L:.2e}  rate={r / 1e6:4.2f} MHz  "
                  f"mu={mu(r, ncoll):.2f}  collisions/triggered crossing={nu:.2f}")

print("\np-Pb: 8.54 vs 5.36 TeV (25 ns p, IP1/5)")
print(f"  potential peak ratio: {ppb_potential(8.54, 1056) / ppb_potential(5.36, 1056):.2f}")
for level in (1e30, 2e30, 1e40):
    a, b = ppb_fill_average(5.36, level), ppb_fill_average(8.54, level)
    tag = "no levelling" if level > 1e35 else f"levelled {level:.0e}"
    print(f"  {tag:16s} fill-averaged L: {a:.2e} -> {b:.2e}  gain x{b / a:.2f}"
          f"  ({a * 86400 / 1e33:.0f} -> {b * 86400 / 1e33:.0f} nb^-1/day at 100% efficiency)")
