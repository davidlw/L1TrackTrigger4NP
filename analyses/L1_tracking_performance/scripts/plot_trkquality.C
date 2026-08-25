// ---------------------------------------------------------------------------
// plot_trkquality.C -- overlay L1 track chi2/ndf distributions, Default vs Dummy.
//
// Each canvas is a shape comparison (unit-normalised, log y) with a dummy/default
// ratio panel underneath. Shapes rather than raw counts, because the two samples
// have slightly different track yields and the question is about the shape.
//
// Usage: root -l -b -q 'plot_trkquality.C("trkquality_qed_mumu.root")'
// ---------------------------------------------------------------------------

#include <cstdio>

namespace {

// Write every figure as both PDF and PNG.
void SaveBoth(TCanvas* c, const char* pdfPath) {
  c->SaveAs(pdfPath);
  TString png(pdfPath);
  png.ReplaceAll(".pdf", ".png");
  c->SaveAs(png);
}

void Overlay(TFile* f, const char* base, const char* xtitle, const char* title,
             const char* out, double xmax, int rebin, bool logy = true) {
  auto* hd = (TH1D*)f->Get(Form("%s_def", base));
  auto* hu = (TH1D*)f->Get(Form("%s_dum", base));
  if (!hd || !hu) {
    printf("[warn] missing %s_def / %s_dum\n", base, base);
    return;
  }
  hd = (TH1D*)hd->Clone(Form("%s_def_c", base));
  hu = (TH1D*)hu->Clone(Form("%s_dum_c", base));
  if (rebin > 1) { hd->Rebin(rebin); hu->Rebin(rebin); }

  // Include the overflow in the normalisation so the fractions are honest, then
  // report separately how much sits beyond the axis.
  double nd = hd->Integral(0, hd->GetNbinsX() + 1);
  double nu = hu->Integral(0, hu->GetNbinsX() + 1);
  double ovd = hd->GetBinContent(hd->GetNbinsX() + 1);
  double ovu = hu->GetBinContent(hu->GetNbinsX() + 1);
  if (nd > 0) hd->Scale(1.0 / nd);
  if (nu > 0) hu->Scale(1.0 / nu);

  hd->SetLineColor(kBlue + 1);
  hd->SetLineWidth(2);
  hu->SetLineColor(kRed + 1);
  hu->SetLineWidth(2);
  hu->SetLineStyle(1);

  auto* c = new TCanvas(Form("c_%s", base), "", 820, 760);
  auto* p1 = new TPad(Form("p1_%s", base), "", 0, 0.30, 1, 1);
  auto* p2 = new TPad(Form("p2_%s", base), "", 0, 0.0, 1, 0.30);
  p1->SetBottomMargin(0.02);
  p1->SetLeftMargin(0.13);
  p2->SetTopMargin(0.02);
  p2->SetBottomMargin(0.34);
  p2->SetLeftMargin(0.13);
  p2->SetGridy();
  p1->Draw();
  p2->Draw();

  p1->cd();
  if (logy) p1->SetLogy();
  hd->GetXaxis()->SetRangeUser(0, xmax);
  hd->GetXaxis()->SetLabelSize(0);
  hd->GetYaxis()->SetTitle("fraction of tracks");
  hd->GetYaxis()->SetTitleSize(0.05);
  hd->GetYaxis()->SetTitleOffset(1.25);
  hd->GetYaxis()->SetLabelSize(0.042);
  double mx = std::max(hd->GetMaximum(), hu->GetMaximum());
  hd->SetMaximum(mx * (logy ? 5 : 1.35));
  if (logy) hd->SetMinimum(std::max(1e-6, 0.2 / std::max(nd, nu)));
  hd->Draw("hist");
  hu->Draw("hist same");

  auto* leg = new TLegend(0.56, 0.70, 0.88, 0.87);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->AddEntry(hd, Form("Default Stub  (mean %.2f)", hd->GetMean()), "l");
  leg->AddEntry(hu, Form("Dummy Stub  (mean %.2f)", hu->GetMean()), "l");
  leg->Draw();

  TLatex tx;
  tx.SetNDC();
  tx.SetTextSize(0.042);
  tx.DrawLatex(0.13, 0.94, title);

  p2->cd();
  auto* r = (TH1D*)hu->Clone(Form("r_%s", base));
  r->Divide(hd);
  r->SetLineColor(kBlack);
  r->SetMarkerColor(kBlack);
  r->SetMarkerStyle(20);
  r->SetMarkerSize(0.7);
  r->GetXaxis()->SetRangeUser(0, xmax);
  r->GetYaxis()->SetRangeUser(0.0, 3.0);
  r->GetYaxis()->SetTitle("dummy / default");
  r->GetYaxis()->SetNdivisions(505);
  r->GetYaxis()->SetTitleSize(0.115);
  r->GetYaxis()->SetTitleOffset(0.52);
  r->GetYaxis()->SetLabelSize(0.10);
  r->GetXaxis()->SetTitle(xtitle);
  r->GetXaxis()->SetTitleSize(0.125);
  r->GetXaxis()->SetTitleOffset(1.15);
  r->GetXaxis()->SetLabelSize(0.10);
  r->Draw("ep");
  auto* l = new TLine(0, 1, xmax, 1);
  l->SetLineStyle(2);
  l->SetLineColor(kRed + 1);
  l->Draw();

  SaveBoth(c, out);

  printf("\n%s\n", title);
  printf("  default: N=%.0f  mean=%.4f  median=%.4f  beyond axis=%.3f%%\n", nd,
         hd->GetMean(), [&] {
           double q = 0.5, v;
           hd->GetQuantiles(1, &v, &q);
           return v;
         }(), nd > 0 ? 100 * ovd / nd : 0);
  printf("  dummy  : N=%.0f  mean=%.4f  median=%.4f  beyond axis=%.3f%%\n", nu,
         hu->GetMean(), [&] {
           double q = 0.5, v;
           hu->GetQuantiles(1, &v, &q);
           return v;
         }(), nu > 0 ? 100 * ovu / nu : 0);
  // Tail fractions -- where the two really separate.
  for (double thr : {1.0, 2.0, 5.0, 10.0}) {
    int b = hd->FindBin(thr);
    double fd = hd->Integral(b, hd->GetNbinsX() + 1);
    double fu = hu->Integral(b, hu->GetNbinsX() + 1);
    printf("    fraction > %5.1f :  default %7.4f   dummy %7.4f   ratio %.2f\n", thr, fd, fu,
           fd > 0 ? fu / fd : 0);
  }
}

// Overlay for a discrete variable (nstub, seed). Integer bins, linear y, and a
// full per-bin fraction table, which is more informative than a smooth curve.
void OverlayDiscrete(TFile* f, const char* base, const char* xtitle, const char* title,
                     const char* out, double xlo, double xhi) {
  auto* hd = (TH1D*)f->Get(Form("%s_def", base));
  auto* hu = (TH1D*)f->Get(Form("%s_dum", base));
  if (!hd || !hu) {
    printf("[warn] missing %s_def / %s_dum\n", base, base);
    return;
  }
  hd = (TH1D*)hd->Clone(Form("%s_def_d", base));
  hu = (TH1D*)hu->Clone(Form("%s_dum_d", base));
  double nd = hd->Integral(0, hd->GetNbinsX() + 1);
  double nu = hu->Integral(0, hu->GetNbinsX() + 1);
  double meand = hd->GetMean(), meanu = hu->GetMean();
  double rmsd = hd->GetRMS(), rmsu = hu->GetRMS();
  if (nd > 0) hd->Scale(1.0 / nd);
  if (nu > 0) hu->Scale(1.0 / nu);

  hd->SetLineColor(kBlue + 1);
  hd->SetLineWidth(2);
  hu->SetLineColor(kRed + 1);
  hu->SetLineWidth(2);

  auto* c = new TCanvas(Form("cd_%s", base), "", 820, 760);
  auto* p1 = new TPad(Form("pd1_%s", base), "", 0, 0.30, 1, 1);
  auto* p2 = new TPad(Form("pd2_%s", base), "", 0, 0.0, 1, 0.30);
  p1->SetBottomMargin(0.02);
  p1->SetLeftMargin(0.13);
  p2->SetTopMargin(0.02);
  p2->SetBottomMargin(0.34);
  p2->SetLeftMargin(0.13);
  p2->SetGridy();
  p1->Draw();
  p2->Draw();

  p1->cd();
  hd->GetXaxis()->SetRangeUser(xlo, xhi);
  hd->GetXaxis()->SetLabelSize(0);
  hd->GetYaxis()->SetTitle("fraction of tracks");
  hd->GetYaxis()->SetTitleSize(0.05);
  hd->GetYaxis()->SetTitleOffset(1.25);
  hd->GetYaxis()->SetLabelSize(0.042);
  hd->SetMaximum(std::max(hd->GetMaximum(), hu->GetMaximum()) * 1.35);
  hd->SetMinimum(0);
  hd->Draw("hist");
  hu->Draw("hist same");

  auto* leg = new TLegend(0.56, 0.70, 0.88, 0.87);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->AddEntry(hd, Form("Default Stub  (mean %.3f)", meand), "l");
  leg->AddEntry(hu, Form("Dummy Stub  (mean %.3f)", meanu), "l");
  leg->Draw();

  TLatex tx;
  tx.SetNDC();
  tx.SetTextSize(0.042);
  tx.DrawLatex(0.13, 0.94, title);

  p2->cd();
  auto* r = (TH1D*)hu->Clone(Form("rd_%s", base));
  r->Divide(hd);
  r->SetLineColor(kBlack);
  r->SetMarkerColor(kBlack);
  r->SetMarkerStyle(20);
  r->SetMarkerSize(0.9);
  r->GetXaxis()->SetRangeUser(xlo, xhi);
  r->GetYaxis()->SetRangeUser(0.0, 2.5);
  r->GetYaxis()->SetTitle("dummy / default");
  r->GetYaxis()->SetNdivisions(505);
  r->GetYaxis()->SetTitleSize(0.115);
  r->GetYaxis()->SetTitleOffset(0.52);
  r->GetYaxis()->SetLabelSize(0.10);
  r->GetXaxis()->SetTitle(xtitle);
  r->GetXaxis()->SetTitleSize(0.125);
  r->GetXaxis()->SetTitleOffset(1.15);
  r->GetXaxis()->SetLabelSize(0.10);
  r->Draw("ep");
  auto* l = new TLine(xlo, 1, xhi, 1);
  l->SetLineStyle(2);
  l->SetLineColor(kRed + 1);
  l->Draw();
  SaveBoth(c, out);

  printf("\n%s\n", title);
  printf("  default: N=%.0f  mean=%.4f  RMS=%.4f\n", nd, meand, rmsd);
  printf("  dummy  : N=%.0f  mean=%.4f  RMS=%.4f\n", nu, meanu, rmsu);
  printf("  %8s %12s %12s %10s\n", xtitle, "default", "dummy", "dum/def");
  for (int i = 1; i <= hd->GetNbinsX(); ++i) {
    double a = hd->GetBinContent(i), b = hu->GetBinContent(i);
    if (a <= 0 && b <= 0) continue;
    printf("  %8.0f %11.4f%% %11.4f%% %10s\n", hd->GetBinCenter(i), 100 * a, 100 * b,
           a > 0 ? Form("%.3f", b / a) : "-");
  }
}

}  // namespace

void plot_trkquality(const char* fname = "trkquality_qed_mumu.root", const char* tag = "") {
  gROOT->SetBatch(true);
  gStyle->SetOptStat(0);

  TFile* f = TFile::Open(fname);
  if (!f || f->IsZombie()) {
    printf("[error] cannot open %s\n", fname);
    return;
  }

  Overlay(f, "chi2dof", "#chi^{2}/ndf", "QED #mu#mu, L1 tracks -- total #chi^{2}/ndf",
          Form("../figures/trk_chi2dof_overlay%s.pdf", tag), 10.0, 2);
  Overlay(f, "chi2rphidof", "#chi^{2}_{r#phi}/ndf",
          "QED #mu#mu, L1 tracks -- r#phi #chi^{2}/ndf",
          Form("../figures/trk_chi2rphidof_overlay%s.pdf", tag), 10.0, 2);
  Overlay(f, "chi2rzdof", "#chi^{2}_{rz}/ndf", "QED #mu#mu, L1 tracks -- rz #chi^{2}/ndf",
          Form("../figures/trk_chi2rzdof_overlay%s.pdf", tag), 10.0, 2);

  OverlayDiscrete(f, "nstub", "trk_nstub", "QED #mu#mu, L1 tracks -- stubs per track",
                  Form("../figures/trk_nstub_overlay%s.pdf", tag), 2.5, 8.5);
  OverlayDiscrete(f, "seed", "trk_seed", "QED #mu#mu, L1 tracks -- seed type",
                  Form("../figures/trk_seed_overlay%s.pdf", tag), -0.5, 12.5);

  printf("\nfigures written to ../figures/\n");
}
