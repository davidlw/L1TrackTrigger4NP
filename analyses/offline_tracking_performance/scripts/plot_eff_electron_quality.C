// ---------------------------------------------------------------------------
// plot_eff_electron_quality.C -- offline tracking efficiency for electrons with
// and without the highPurity requirement, against muons (highPurity), from
// three offline_perf.C outputs made with pdgSel = 11 / 11 / 13.
//
// The all-generalTracks curve is the ceiling for anything built on top of the
// tracks (lowPtGsfElectrons re-seed from generalTracks): the gap between the
// two electron curves is what an electron-specific refit could recover, the
// gap to the muons is electrons with no track at all.
//
// Usage (from scripts/):
//   root -l -b -q 'plot_eff_electron_quality.C("../output/offperf_qedee_hp_prim_e.root",
//        "../output/offperf_qedee_all_prim_e.root","../output/offperf_qedmumu_hp_prim_mu.root",
//        "../figures/eff_electron_quality_offline")'
// ---------------------------------------------------------------------------

#include <cstdio>

namespace {

TGraphAsymmErrors* Eff(TFile* f, const char* var, int color, int marker, int style = 1) {
  auto* hn = (TH1D*)f->Get(Form("num_%s_off", var));
  auto* hd = (TH1D*)f->Get(Form("den_%s_off", var));
  if (!hn || !hd) {
    printf("[error] %s: missing num/den_%s_off\n", f->GetName(), var);
    return nullptr;
  }
  auto* g = new TGraphAsymmErrors(hn, hd, "cl=0.683 b(1,1) mode");
  g->SetLineColor(color);
  g->SetMarkerColor(color);
  g->SetMarkerStyle(marker);
  g->SetMarkerSize(1.2);
  g->SetLineWidth(2);
  g->SetLineStyle(style);
  return g;
}

void Panel(TFile* fEhp, TFile* fEall, TFile* fMu, const char* var, const char* xTitle, double xlo, double xhi,
           bool logx, const char* note, double legx) {
  if (logx)
    gPad->SetLogx();
  gPad->SetLeftMargin(0.12);
  gPad->SetGridy();
  auto* fr = gPad->DrawFrame(xlo, 0, xhi, 1.12);
  fr->GetXaxis()->SetTitle(xTitle);
  fr->GetYaxis()->SetTitle("offline tracking efficiency");
  fr->GetYaxis()->SetTitleOffset(1.3);
  if (logx) {
    fr->GetXaxis()->SetMoreLogLabels();
    fr->GetXaxis()->SetNoExponent();
  }
  auto* gMu = Eff(fMu, var, kAzure + 2, 20);
  auto* gEa = Eff(fEall, var, kRed + 1, 25, 2);
  auto* gEh = Eff(fEhp, var, kRed + 1, 21);
  if (!gMu || !gEa || !gEh)
    return;
  gMu->Draw("PZ same");
  gEa->Draw("PZ same");
  gEh->Draw("PZ same");
  auto* leg = new TLegend(legx, 0.16, legx + 0.42, 0.40);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetHeader("STARlight #gamma#gamma#rightarrow#it{l}^{+}#it{l}^{-}, primaries");
  leg->AddEntry(gMu, "#mu^{#pm}, highPurity", "pl");
  leg->AddEntry(gEa, "e^{#pm}, all generalTracks", "pl");
  leg->AddEntry(gEh, "e^{#pm}, highPurity", "pl");
  leg->Draw();
  TLatex tx;
  tx.SetNDC();
  tx.SetTextSize(0.038);
  tx.DrawLatex(0.14, 0.93, note);
}

}  // namespace

void plot_eff_electron_quality(const char* fileEhp = "../output/offperf_qedee_hp_prim_e.root",
                               const char* fileEall = "../output/offperf_qedee_all_prim_e.root",
                               const char* fileMu = "../output/offperf_qedmumu_hp_prim_mu.root",
                               const char* outBase = "../figures/eff_electron_quality_offline") {
  gROOT->SetBatch(true);
  gStyle->SetOptStat(0);
  auto *a = TFile::Open(fileEhp), *b = TFile::Open(fileEall), *m = TFile::Open(fileMu);
  if (!a || !b || !m)
    return;
  auto* c = new TCanvas("c", "", 1300, 560);
  c->Divide(2, 1);
  c->cd(1);
  Panel(a, b, m, "pt", "tracking particle p_{T} [GeV]", 0.28, 12, true, "|#eta| < 2.4", 0.46);
  c->cd(2);
  Panel(a, b, m, "eta1", "tracking particle #eta", -2.4, 2.4, false, "p_{T} > 0.6 GeV", 0.30);
  c->SaveAs(Form("%s.pdf", outBase));
  c->SaveAs(Form("%s.png", outBase));
}
