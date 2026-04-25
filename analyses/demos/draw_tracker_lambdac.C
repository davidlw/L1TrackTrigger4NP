#include <iostream>
#include <vector>
#include <cmath>
#include "TCanvas.h"
#include "TH1F.h"
#include "TLine.h"
#include "TPolyLine.h"
#include "TEllipse.h"
#include "TLatex.h"
#include "TMarker.h"

static const double ot_radii_mm[] = {250.0, 360.0, 500.0, 680.0, 890.0, 1080.0};

// -------------------------------------------------------------------
// HELPER: Physical Bending with Stubs
// -------------------------------------------------------------------
void draw_track_and_stubs(double pT_GeV, double phi_start_deg, int color, int charge) {
    double B = 3.8;
    double R_mm = (pT_GeV * 1000.0) / (0.3 * B);
    double target_r_mm = 1100.0;
    double alpha = phi_start_deg * M_PI / 180.0;
    
    double cos_phi = 1.0 - (pow(target_r_mm, 2) / (2.0 * pow(R_mm, 2)));
    double max_phi = (cos_phi < -1.0) ? M_PI : acos(cos_phi);
    
    int n = 250;
    std::vector<double> vx(n), vy(n);
    for(int i=0; i<n; ++i) {
        double phi = (i * max_phi) / (n - 1); 
        // +1 charge = Clockwise (CW), -1 charge = Counter-Clockwise (CCW)
        vx[i] = charge * R_mm * (sin(alpha + charge * phi) - sin(alpha));
        vy[i] = charge * R_mm * (cos(alpha) - cos(alpha + charge * phi));
    }
    TPolyLine* line = new TPolyLine(n, &vx[0], &vy[0]);
    line->SetLineColor(color); 
    line->SetLineWidth(4); 
    line->Draw();

    for (int layer = 0; layer < 6; ++layer) {
        double r = ot_radii_mm[layer];
        double c_int = 1.0 - (pow(r, 2) / (2.0 * pow(R_mm, 2)));
        if (c_int < -1.0) continue;
        double p_int = acos(c_int);
        double ix = charge * R_mm * (sin(alpha + charge * p_int) - sin(alpha));
        double iy = charge * R_mm * (cos(alpha) - cos(alpha + charge * p_int));
        TMarker *m = new TMarker(ix, iy, 29);
        m->SetMarkerColor(color); m->SetMarkerSize(1.6); m->Draw();
    }
}

void draw_tracker_lambdac() {
    TCanvas *c1 = new TCanvas("c1", "CMS Tracker Lambda_c Analysis", 1000, 1000);
    TH1F* frame = c1->DrawFrame(-1350, -1350, 1350, 1350);
    frame->SetTitle("UPC #Lambda_{c}^{+} #to p K^{-} #pi^{+} (p_{T}^{#Lambda_c} = 0);x [mm];y [mm]");

    // 1. Shaded Visibility Region: N2 Window (20 to 100 deg)
    TEllipse *shade = new TEllipse(0, 0, 1150, 1150, 20, 100);
    shade->SetFillColorAlpha(kYellow-9, 0.4); 
    shade->SetLineColor(kYellow+2);
    shade->SetLineStyle(2); 
    shade->Draw();

    // 2. Full Geometry (Small module segments)
    for (int layer = 0; layer < 6; ++layer) {
        double r = ot_radii_mm[layer];
        int n_mod = (int)(2 * M_PI * r / 35.0);
        double d_phi = (2.0 * M_PI / n_mod) * 0.45;
        for (int m = 0; m < n_mod; m++) {
            double phi = (m * 2.0 * M_PI) / n_mod;
            double cur_r = (m % 2 == 0) ? (r - 8) : (r + 8);
            TLine *mod = new TLine(cur_r*cos(phi-d_phi), cur_r*sin(phi-d_phi), cur_r*cos(phi+d_phi), cur_r*sin(phi+d_phi));
            mod->SetLineColor((layer < 3) ? kBlue-10 : kRed-10); 
            mod->Draw();
        }
    }

    // 3. Sector Boundaries
    for (int i = 0; i < 9; ++i) {
        double angle = i * (40.0 * M_PI / 180.0);
        TLine *b = new TLine(0, 0, 1280 * cos(angle), 1280 * sin(angle));
        b->SetLineColor(kGray+1); b->Draw(); 
        TLatex *l = new TLatex(1200 * cos(angle+0.1), 1200 * sin(angle+0.1), Form("N%d", i+1));
        l->SetTextSize(0.018); l->Draw();
    }

    // --- 4. THE Lambda_c DECAY (3 tracks) ---
    TMarker *v = new TMarker(0, 0, 20); v->Draw();

    // Proton (Orange, +1): 0.9 GeV, starts at 60 deg (N2 center), bends CW
    draw_track_and_stubs(0.90, 60.0, kOrange+7, 1);   
    
    // Kaon- (Blue, -1): 0.7 GeV, starts at 180 deg (N5), bends CCW
    draw_track_and_stubs(0.70, 180.0, kAzure+7, -1);  

    // Pion+ (Red, +1): 0.5 GeV, starts at 300 deg (N8), bends CW
    draw_track_and_stubs(0.50, 300.0, kRed+1, 1);
/*
    // 5. Annotations
    TLatex *t = new TLatex(); t->SetNDC(); t->SetTextSize(0.028);
    t->DrawLatex(0.13, 0.85, "#Lambda_{c}^{+} #to p (Orange) + K^{-} (Blue) + #pi^{+} (Red)");
    t->DrawLatex(0.13, 0.81, "3-Body decay results in 3 isolated orphans across different Nonants.");
    t->DrawLatex(0.13, 0.77, "GTT cross-sector correlation is mandatory for mass reconstruction.");
*/
    c1->Update();
}
