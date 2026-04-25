#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TLine.h>
#include <TLegend.h>
#include <TPaveText.h>
#include <TText.h>
#include <vector>
#include <map>
#include <iostream>

void OTModuleOccupancy() {
    TFile *f = TFile::Open("StarLightJpsiPhase2_PrivateMC_L1TrackNtuplev2.root");
    if (!f || f->IsZombie()) return;
    TTree *tree = (TTree*)f->Get("L1TrackHitNtupleMaker/eventTree");

    std::vector<int> *c_lyr = 0, *c_isPS = 0, *c_sens = 0, *c_half = 0, *c_chip = 0;
    std::vector<unsigned int> *c_detid = 0;
    
    tree->SetBranchAddress("cluster_layer", &c_lyr);
    tree->SetBranchAddress("cluster_isPS", &c_isPS);
    tree->SetBranchAddress("cluster_sensor", &c_sens);
    tree->SetBranchAddress("cluster_detid", &c_detid);
    tree->SetBranchAddress("cluster_halfModule", &c_half);
    tree->SetBranchAddress("cluster_chipId", &c_chip);

    std::map<TString, std::map<int, std::vector<TH1F*>>> masterMap;
    struct StatBox { long nC=0; long fC=0; long nM=0; long fM=0; };
    std::map<TString, std::map<int, StatBox>> masterStats;

    Long64_t nentries = tree->GetEntries();
    std::cout << "Processing " << nentries << " events..." << std::endl;

    for (Long64_t i = 0; i < nentries; i++) {
        tree->GetEntry(i);
        if (!c_lyr || c_lyr->empty()) continue;

        std::map<uint64_t, int> chip_hits;
        std::map<uint64_t, bool> chip_isPS_map;

        for (size_t j = 0; j < c_lyr->size(); j++) {
            uint64_t key = ((uint64_t)((*c_lyr)[j]) << 48) | ((uint64_t)((*c_detid)[j]) << 16) | 
                           ((uint64_t)((*c_half)[j]) << 12) | ((uint64_t)((*c_chip)[j]) << 4) | ((*c_sens)[j]);
            chip_hits[key]++;
            chip_isPS_map[key] = (bool)(*c_isPS)[j];
        }

        std::map<uint64_t, std::pair<int, int>> mod_agg;
        for (auto const& [key, n_hit] : chip_hits) {
            int lyr = (key >> 48) & 0xFF;
            int sns = key & 0xF;
            bool isPS = chip_isPS_map[key];
            
            TString typeStr = isPS ? "PS" : "2S";
            TString regStr = (lyr <= 10) ? "Barrel" : (lyr < 110 ? "EndcapNeg" : "EndcapPos");
            TString figKey = Form("%s_%s_S%d", typeStr.Data(), regStr.Data(), sns);

            int chipLim = isPS ? (sns == 0 ? 5 : 8) : 3;

            masterStats[figKey][lyr].nC++;
            if (n_hit > chipLim) masterStats[figKey][lyr].fC++;

            if (masterMap[figKey].find(lyr) == masterMap[figKey].end()) {
                TString lbl = (lyr > 10) ? Form("Disk %d", lyr%10) : Form("Layer %d", lyr);
                masterMap[figKey][lyr].push_back(new TH1F(Form("h_%s_%d", figKey.Data(), lyr), lbl+";Clusters;Entries", 45, 0, 45));
                masterMap[figKey][lyr].push_back(new TH1F(Form("h_%s_%d_Tr", figKey.Data(), lyr), "Trunc", 45, 0, 45));
                for(auto h : masterMap[figKey][lyr]) h->SetDirectory(0);
            }
            uint64_t modKey = key & 0xFFFFFFFFFFFF000F;
            mod_agg[modKey].first += n_hit;
            mod_agg[modKey].second += std::min(n_hit, chipLim);
        }

        for (auto const& [mKey, sums] : mod_agg) {
            int lyr = (mKey >> 48) & 0xFF;
            int sns = mKey & 0xF;
            for (auto & [fKey, lyrMap] : masterMap) {
                if (fKey.EndsWith(Form("S%d", sns)) && lyrMap.count(lyr)) {
                    int modLim = fKey.Contains("PS") ? 17 : 16;
                    masterStats[fKey][lyr].nM++;
                    if (sums.first > modLim) masterStats[fKey][lyr].fM++;
                    lyrMap[lyr][0]->Fill(sums.first);
                    lyrMap[lyr][1]->Fill(sums.second);
                }
            }
        }
    }

    // --- Dynamic Canvas Scaling ---
    for (auto const& [fKey, lyrMap] : masterMap) {
        int nPads = lyrMap.size();
        int canW = 1500;
        int canH = (nPads <= 3) ? 500 : 900; // Shorter canvas for Barrel (3 layers)
        
        TCanvas *c = new TCanvas("c"+fKey, fKey, canW, canH);
        
        if (nPads <= 3) c->Divide(3, 1); // 1 row, 3 columns
        else           c->Divide(3, 2); // 2 rows, 3 columns

        int iPad = 1;
        for (auto const& [lyr, h] : lyrMap) {
            c->cd(iPad++); gPad->SetLogy();
            h[0]->SetLineColor(kBlue); h[0]->SetLineWidth(2); h[0]->Draw("HIST");
            h[1]->SetLineColor(kRed); h[1]->SetLineStyle(7); h[1]->SetLineWidth(2); h[1]->Draw("HIST SAME");

            double mLim = fKey.Contains("PS") ? 17.0 : 16.0;
            TLine *l = new TLine(mLim, 0.5, mLim, h[0]->GetMaximum()*1.5);
            l->SetLineColor(kGreen+2); l->SetLineWidth(3); l->Draw();

            TPaveText *pt = new TPaveText(0.48, 0.55, 0.88, 0.75, "NDC");
            pt->SetFillColor(0); pt->SetBorderSize(1); pt->SetTextSize(0.045);
            auto s = masterStats[fKey][lyr];
            TText *tx1 = pt->AddText(Form("Chip Fail: %.2f%%", (s.nC > 0 ? 100.0*s.fC/s.nC : 0)));
            tx1->SetTextColor(kRed);
            TText *tx2 = pt->AddText(Form("Mod Fail:  %.2f%%", (s.nM > 0 ? 100.0*s.fM/s.nM : 0)));
            tx2->SetTextColor(kGreen+2);
            pt->Draw();
        }
        c->SaveAs("Genesis_UPC_Fig_" + fKey + ".png");
    }
}
