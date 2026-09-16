# Trigger summary ("money plot")

Event trigger efficiency per physics process, for the default L1 tracking, the
new low-pT L1 tracking, and offline highPurity tracks, with an L1/offline ratio
panel.

```
eff = N(events with >= k tracks, pT > pT_algo, |eta| < 2.4)
      / N(events with >= k PRIMARY charged particles, pT > 0.4 GeV, |eta| < 2.4)
```

The denominator is the fiducial truth at the low-pT threshold, common to all
three series, not all generated events -- so acceptance is factored out and
processes are comparable. The numerator threshold is each trigger's own:
default L1 fires on pT > 2 GeV (its tracking does not exist below), low-pT L1
and offline on pT > 0.4 GeV. The default column therefore folds threshold
acceptance and tracking inefficiency together, which is the intended message.

Signals are grouped (`--- title` lines in the table): UPC gamma-gamma, UPC
gamma-A, hadronic MB. For exclusive signals "event has >= k particles" means
"all decay daughters are in the fiducial cut".

```
scripts/money_plot.C          drawing macro; reads the table below
scripts/money_plot_input.txt  one row per process: label Nfid Ndefault NlowpT Noffline
figures/                      money_plot.pdf/.png
```

```
cd scripts
root -l -b -q 'money_plot.C("money_plot_input.txt","../figures/money_plot",1,0.4,2.0,0.4,0.4)'
```

Arguments: k, pT_fid (denominator), then the numerator thresholds for default
L1, low-pT L1, offline. They only label the canvas; the counts must match.

The shipped input table is a placeholder (the macro draws a watermark while
the word `TEMPLATE` is present in the file). Fill the counts from the
threshold-scan histograms of `../offline_tracking_performance/scripts/mbeff_offline.C`
(offline + truth) and `../L1_tracking_performance/scripts/compare_mbeff.C` (L1),
reading the numerator/denominator bin at the chosen `pTmin` and `k`.
