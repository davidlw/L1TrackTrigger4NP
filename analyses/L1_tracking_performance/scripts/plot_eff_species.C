// ---------------------------------------------------------------------------
// plot_eff_species.C -- tracking efficiency for one particle species against
// another, from two efficiency outputs made with pdgSel (e.g. 11 on the
// STARlight e+e- sample and 13 on the mu+mu- sample).
//
// Works for both producers, they write the same histogram names up to a slot:
//   L1       compare_eff.C   -> num_pt_def / den_pt_def   (slot "def"; only the
//            default half is read, so put the same directory in both slots or
//            ignore the dummy half)          eta panels: eta1 = pT>1, eta2 = pT>2
//   offline  offline_perf.C  -> num_pt_off / den_pt_off   (slot "off")
//                                            eta panels: eta0..eta3 = pT>0.3,0.6,1,2
//
// Left: vs pT (log x). Right: vs eta for the chosen threshold panel.
// Numbers are dumped as text.
//
// Usage:
//   root -l -b -q 'plot_eff_species.C("../output/eff_qed_ee_electrons_default.root","e^{#pm}",
//                                     "../output/eff_qed_mumu_muons_default.root","#mu^{#pm}",
//                                     "Default stub", "../figures/eff_electron_vs_muon_default")'
//   offline:  ... , "Offline highPurity", "../figures/eff_electron_vs_muon_offline",
//                   "off", "eta1", "p_{T} > 0.6 GeV", 0.28, 12, "offline tracking efficiency")
// ---------------------------------------------------------------------------

#include <cstdio>

namespace {

const char* gSlot = "def";

TGraphAsymmErrors* Eff(TFile* f, const char* var, int color, int marker) {
  auto* hn = (TH1D*)f->Get(Form("num_%s_%s", var, gSlot));
  auto* hd = (TH1D*)f->Get(Form("den_%s_%s", var, gSlot));
  if (!hn || !hd) {
    printf("[error] %s: missing num/den_%s_%s\n", f->GetName(), var, gSlot);
    return nullptr;
  }
  auto* g = new TGraphAsymmErrors(hn, hd, "cl=0.683 b(1,1) mode");
  g->SetLineColor(color);
  g->SetMarkerColor(color);
  g->SetMarkerStyle(marker);
  g->SetMarkerSize(1.3);
  g->SetLineWidth(2);
  return g;
}

void Dump(TFile* fa, const char* la, TFile* fb, const char* lb, const char* var) {
  auto *na = (TH1D*)fa->Get(Form("num_%s_%s", var, gSlot)), *da = (TH1D*)fa->Get(Form("den_%s_%s", var, gSlot));
  auto *nb = (TH1D*)fb->Get(Form("num_%s_%s", var, gSlot)), *db = (TH1D*)fb->Get(Form("den_%s_%s", var, gSlot));
  printf("\n%-14s %9s %8s   %9s %8s\n", var, Form("N(%s)", la), "eff", Form("N(%s)", lb), "eff");
  for (int i = 1; i <= na->GetNbinsX(); ++i) {
    double Da = da->GetBinContent(i), Db = db->GetBinContent(i);
    if (Da < 20 && Db < 20) continue;
    printf("%6.2f-%-6.2f %9.0f %8.3f   %9.0f %8.3f\n", na->GetBinLowEdge(i), na->GetBinLowEdge(i + 1), Da,
           Da > 0 ? na->GetBinContent(i) / Da : 0, Db, Db > 0 ? nb->GetBinContent(i) / Db : 0);
  }
}

}  // namespace

void plot_eff_species(const char* fileA, const char* labelA, const char* fileB, const char* labelB,
                      const char* algo = "Default stub", const char* outBase = "../figures/eff_species",
                      const char* slot = "def", const char* etaVar = "eta2",
                      const char* etaLabel = "p_{T} > 2 GeV", double ptLo = 0.9, double ptHi = 12,
                      const char* yTitle = "L1 tracking efficiency") {
  gSlot = slot;
  gROOT->SetBatch(true);
  gStyle->SetOptStat(0);
  auto *fa = TFile::Open(fileA), *fb = TFile::Open(fileB);
  if (!fa || !fb)
    return;

  const int cA = kRed + 1, cB = kAzure + 2;
  auto* c = new TCanvas("c", "", 1300, 560);
  c->Divide(2, 1);

  // ---- vs pT
  c->cd(1);
  gPad->SetLogx();
  gPad->SetLeftMargin(0.12);
  gPad->SetGridy();
  auto *gA = Eff(fa, "pt", cA, 21), *gB = Eff(fb, "pt", cB, 20);
  if (!gA || !gB)
    return;
  auto* fr = gPad->DrawFrame(ptLo, 0, ptHi, 1.12);
  fr->GetXaxis()->SetTitle("tracking particle p_{T} [GeV]");
  fr->GetYaxis()->SetTitle(yTitle);
  fr->GetXaxis()->SetMoreLogLabels();
  fr->GetXaxis()->SetNoExponent();
  fr->GetYaxis()->SetTitleOffset(1.3);
  gB->Draw("PZ same");
  gA->Draw("PZ same");
  auto* leg = new TLegend(0.50, 0.18, 0.88, 0.36);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetHeader(algo);
  leg->AddEntry(gA, labelA, "pl");
  leg->AddEntry(gB, labelB, "pl");
  leg->Draw();
  TLatex tx;
  tx.SetNDC();
  tx.SetTextSize(0.038);
  tx.DrawLatex(0.14, 0.93, "|#eta| < 2.4, N(tp_{nmatch} > 0) / N(tp)");

  // ---- vs eta, pT > 2
  c->cd(2);
  gPad->SetLeftMargin(0.12);
  gPad->SetGridy();
  auto *eA = Eff(fa, etaVar, cA, 21), *eB = Eff(fb, etaVar, cB, 20);
  if (!eA || !eB)
    return;
  auto* fr2 = gPad->DrawFrame(-2.4, 0, 2.4, 1.12);
  fr2->GetXaxis()->SetTitle("tracking particle #eta");
  fr2->GetYaxis()->SetTitle(yTitle);
  fr2->GetYaxis()->SetTitleOffset(1.3);
  eB->Draw("PZ same");
  eA->Draw("PZ same");
  auto* leg2 = new TLegend(0.35, 0.18, 0.75, 0.36);
  leg2->SetBorderSize(0);
  leg2->SetFillStyle(0);
  leg2->SetHeader(algo);
  leg2->AddEntry(eA, labelA, "pl");
  leg2->AddEntry(eB, labelB, "pl");
  leg2->Draw();
  tx.DrawLatex(0.14, 0.93, etaLabel);

  c->SaveAs(Form("%s.pdf", outBase));
  c->SaveAs(Form("%s.png", outBase));

  Dump(fa, "A", fb, "B", "pt");
  Dump(fa, "A", fb, "B", etaVar);
}
