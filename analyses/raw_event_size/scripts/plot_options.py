#!/usr/bin/env python3
"""Stored event size per readout option, for pp PU200, hadronic PbPb and UPC PbPb.

One panel per event class (shared MB axis); one bar per option, stacked by
subdetector group, with the same colours as plot_sizes.py.

Usage (any ROOT with PyROOT):  python3 plot_options.py results.json outdir
"""
import json, os, sys
import ROOT
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)

res = json.load(open(sys.argv[1]))
outdir = sys.argv[2] if len(sys.argv) > 2 else "."
MB = 1e6

GROUPS = [
    ("Tracker (IT+OT)", ["Inner Tracker", "Outer Tracker"], "#2a78d6"),
    ("HGCAL", ["HGCAL"], "#eb6834"),
    ("MTD (BTL+ETL)", ["MTD BTL", "MTD ETL"], "#1baf7a"),
    ("ECAL barrel + HCAL", ["ECAL Barrel", "HCAL (HB+HO+HF)"], "#eda100"),
    ("Muons", ["Muons (DT+CSC+GEM+RPC)"], "#e87ba4"),
    ("L1 trigger + TPGs", ["L1 trigger + TPGs"], "#008300"),
]
PANELS = [("pp, <PU> = 200 (t#bar{t} RelVal)", "pp"),
          ("PbPb hadronic (Hydjet MB)", "hadronic"),
          ("PbPb UPC (STARlight QED #mu#mu)", "upc")]
opts = res["options"]
names = list(opts["pp"].keys())          # option order as defined in size_model.py
rows = list(reversed(names))             # first option at the top


def label(n):
    return n.replace(" [recommended]", "  #bf{[recommended]}")


XMAX = 10
c = ROOT.TCanvas("c", "", 1250, 1150)
top, bottom, legend_h = 0.94, 0.07, 0.07
ph = (top - bottom - legend_h) / len(PANELS)
keep = []
for ip, (title, cls) in enumerate(PANELS):
    y1 = top - (ip + 1) * ph
    pad = ROOT.TPad(f"p{ip}", "", 0.0, y1, 1.0, y1 + ph)
    pad.SetLeftMargin(0.40); pad.SetRightMargin(0.10)
    pad.SetTopMargin(0.14); pad.SetBottomMargin(0.12 if ip == len(PANELS) - 1 else 0.06)
    pad.SetGridx(); pad.Draw(); pad.cd(); keep.append(pad)
    stack = ROOT.THStack(f"s{ip}", ""); keep.append(stack)
    for g, dets, col in GROUPS:
        h = ROOT.TH1F(f"h{ip}{g}", "", len(rows), 0, len(rows))
        for i, n in enumerate(rows):
            h.SetBinContent(i + 1, sum(opts[cls][n][d] for d in dets) / MB)
            h.GetXaxis().SetBinLabel(i + 1, label(n))
        h.SetFillColor(ROOT.TColor.GetColor(col)); h.SetLineColor(ROOT.kWhite); h.SetLineWidth(2)
        h.SetBarWidth(0.7); h.SetBarOffset(0.15)
        stack.Add(h); keep.append(h)
    stack.SetMaximum(XMAX)
    stack.Draw("hbar")
    stack.GetXaxis().SetLabelSize(0.085); stack.GetXaxis().SetTickLength(0)
    stack.GetYaxis().SetLabelSize(0.075 if ip == len(PANELS) - 1 else 0)
    stack.GetYaxis().SetNdivisions(510)
    tex = ROOT.TLatex(); tex.SetTextSize(0.075); tex.SetTextAlign(12)
    for i, n in enumerate(rows):
        tot = opts[cls][n]["TOTAL"] / MB
        tex.DrawLatex(tot + 0.12, i + 0.5, "%.2f" % tot)
    t = ROOT.TLatex(); t.SetNDC(); t.SetTextSize(0.09); t.SetTextFont(62)
    t.DrawLatex(0.40, 0.9, title)
    c.cd()

# x-axis title, legend, header
c.cd()
xt = ROOT.TLatex(); xt.SetNDC(); xt.SetTextSize(0.022); xt.SetTextAlign(31)
xt.DrawLatex(0.90, bottom + legend_h - 0.035, "mean stored raw event size per event [MB]")
leg = ROOT.TLegend(0.05, 0.005, 0.98, bottom + 0.02); leg.SetNColumns(6)
leg.SetBorderSize(0); leg.SetTextSize(0.017); leg.SetFillStyle(0)
for g, dets, col in GROUPS:
    b = ROOT.TH1F("leg" + g, "", 1, 0, 1); b.SetFillColor(ROOT.TColor.GetColor(col)); b.SetLineColor(ROOT.kWhite)
    leg.AddEntry(b, g, "f"); keep.append(b)
leg.Draw()
hdr = ROOT.TLatex(); hdr.SetNDC(); hdr.SetTextSize(0.019)
hdr.DrawLatex(0.02, 0.978, "CMS Phase-2 raw event size model: channel-level ZS everywhere except HCAL")
hdr.DrawLatex(0.02, 0.957, "Risk not shown: ECAL barrel without its planned sparse readout adds +2.0 MB to every bar")
c.SaveAs(os.path.join(outdir, "raw_size_options.png"))
c.SaveAs(os.path.join(outdir, "raw_size_options.pdf"))
