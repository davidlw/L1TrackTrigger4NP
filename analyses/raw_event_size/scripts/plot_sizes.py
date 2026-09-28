#!/usr/bin/env python3
"""Plots for the raw event size estimate, from size_model_results.json.

  raw_size_breakdown.png : mean size per sample, stacked by subdetector group,
                           with the low/high variant band on the total
  raw_size_vs_activity.png : Hydjet per-event size vs IT pixel digis

Usage (any ROOT with PyROOT):  python3 plot_sizes.py results.json outdir
"""
import json, os, sys
import ROOT
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)

res = json.load(open(sys.argv[1]))
outdir = sys.argv[2] if len(sys.argv) > 2 else "."
MB = 1e6

# Subdetector groups, in fixed categorical order (colour follows the group).
GROUPS = [
    ("Tracker (IT+OT)", ["Inner Tracker", "Outer Tracker"], "#2a78d6"),
    ("HGCAL", ["HGCAL"], "#eb6834"),
    ("MTD (BTL+ETL)", ["MTD BTL", "MTD ETL"], "#1baf7a"),
    ("ECAL barrel + HCAL", ["ECAL Barrel", "HCAL (HB+HO+HF)"], "#eda100"),
    ("Muons", ["Muons (DT+CSC+GEM+RPC)"], "#e87ba4"),
    ("L1 trigger + TPGs", ["L1 trigger + TPGs"], "#008300"),
]
pu200_key = next(k for k in res["results"] if k.startswith("PU200"))
SAMPLES = [("UPC (STARlight QED #mu#mu)", "UPC (STARlight QED mumu)"),
           ("PbPb MB (Hydjet)", "Hydjet PbPb MB"),
           ("pp PU140 t#bar{t} (ref.)", "PU140 ttbar (validation)"),
           ("pp PU200 t#bar{t} (RelVal)", pu200_key)]

c = ROOT.TCanvas("c", "", 1150, 620)
c.SetLeftMargin(0.21); c.SetRightMargin(0.26); c.SetBottomMargin(0.13)
stack = ROOT.THStack("s", "")
keep = []
for name, dets, col in GROUPS:
    h = ROOT.TH1F("h" + name, "", len(SAMPLES), 0, len(SAMPLES))
    for i, (lab, key) in enumerate(SAMPLES):
        h.SetBinContent(i + 1, sum(res["results"][key]["mean"][d] for d in dets) / MB)
        h.GetXaxis().SetBinLabel(i + 1, lab)
    h.SetFillColor(ROOT.TColor.GetColor(col)); h.SetLineColor(ROOT.kWhite); h.SetLineWidth(2)
    h.SetBarWidth(0.6); h.SetBarOffset(0.2)
    stack.Add(h); keep.append((name, h))
stack.SetMaximum(10)
stack.Draw("hbar")
stack.GetXaxis().SetLabelSize(0.045); stack.GetXaxis().SetTickLength(0)
stack.GetYaxis().SetLabelSize(0.04)
xt = ROOT.TLatex(); xt.SetNDC(); xt.SetTextSize(0.042); xt.SetTextAlign(31)
xt.DrawLatex(0.74, 0.03, "mean raw event size [MB]")
# total and low/high band
band = ROOT.TGraphAsymmErrors()
tex = ROOT.TLatex(); tex.SetTextSize(0.035); tex.SetTextAlign(12)
for i, (lab, key) in enumerate(SAMPLES):
    r = res["results"][key]
    band.SetPoint(i, r["total"] / MB, i + 0.5)
    band.SetPointError(i, (r["total"] - r["total_low"]) / MB, (r["total_high"] - r["total"]) / MB, 0, 0)
    tex.DrawLatex(r["total_high"] / MB + 0.15, i + 0.5, "%.1f MB" % (r["total"] / MB))
band.SetLineWidth(2); band.SetMarkerStyle(20); band.SetMarkerSize(0.8)
band.Draw("P same")
leg = ROOT.TLegend(0.755, 0.3, 0.995, 0.9); leg.SetBorderSize(0); leg.SetTextSize(0.033)
for name, h in keep: leg.AddEntry(h, name, "f")
leg.AddEntry(band, "total, low/high variant", "lep")
leg.Draw()
t2 = ROOT.TLatex(); t2.SetNDC(); t2.SetTextSize(0.035)
t2.DrawLatex(0.21, 0.93, "CMS Phase-2 simulation, raw size model calibrated to the DAQ TDR")
c.SaveAs(os.path.join(outdir, "raw_size_breakdown.png"))
c.SaveAs(os.path.join(outdir, "raw_size_breakdown.pdf"))

# size vs activity (Hydjet)
r = res["results"]["Hydjet PbPb MB"]; u = res["results"]["UPC (STARlight QED mumu)"]
c2 = ROOT.TCanvas("c2", "", 800, 600)
c2.SetLeftMargin(0.13); c2.SetBottomMargin(0.13); c2.SetGridy()
g = ROOT.TGraph()
for i, (n, t) in enumerate(zip(r["per_event_itDigi"], r["per_event_total"])):
    g.SetPoint(i, n / 1e3, t / MB)
g.SetMarkerStyle(20); g.SetMarkerSize(0.8); g.SetMarkerColor(ROOT.TColor.GetColor("#2a78d6"))
g.SetTitle(";Inner Tracker pixel digis per event [10^{3}];modelled raw event size [MB]")
g.Draw("AP")
g.GetXaxis().SetLimits(-20, 1000); g.SetMinimum(0); g.SetMaximum(11)
for ax in (g.GetXaxis(), g.GetYaxis()): ax.SetTitleSize(0.045); ax.SetLabelSize(0.04)
lines = []
for y, txt, st in ((r["total"] / MB, "PbPb MB mean %.1f MB" % (r["total"] / MB), 2),
                   (u["total"] / MB, "UPC mean %.1f MB" % (u["total"] / MB), 7)):
    ln = ROOT.TLine(-20, y, 1000, y); ln.SetLineStyle(st); ln.SetLineWidth(2); ln.Draw(); lines.append(ln)
    tx = ROOT.TLatex(600, y + 0.2, txt); tx.SetTextSize(0.035); tx.Draw(); lines.append(tx)
t3 = ROOT.TLatex(); t3.SetNDC(); t3.SetTextSize(0.035)
t3.DrawLatex(0.13, 0.92, "Hydjet PbPb MB, %d events (one point per event)" % len(r["per_event_total"]))
c2.SaveAs(os.path.join(outdir, "raw_size_vs_activity.png"))
c2.SaveAs(os.path.join(outdir, "raw_size_vs_activity.pdf"))
