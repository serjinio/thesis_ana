#include <iostream>

#include "TChain.h"


void theta_phi_comparison() {
    TChain chain{"scattree"};

    for (int i = 272; i < 274; i++){
        TString name;
        name.Form("scattree_rootfiles/scattree_run%i.root", i);
        std::cout << "Adding file to chain: " << name << std::endl;
        chain.Add(name);
    }

    chain.Print();

    TFile f_hist_out("out/theta_phi_comparison.root", "RECREATE");
    TCanvas c1("c1");
    gStyle->SetOptStat(1111111);

    c1.Divide(1);

    c1.cd(1);
    //gPad->SetLogz();
    chain.Draw("s1dc_theta*57.3:fdc0_theta*57.3 >>s1dc_fdc0_theta(400,-20,20,400,-20,20)",
               "triggers[5]==1", "colz");
    c1.Print("out/theta_phi_comparison.pdf(", "pdf");
    //f_hist_out.Write();

    c1.Clear();
    c1.Divide(1);

    c1.cd(1);
    chain.Draw("s1dc_phi*57.3:fdc0_phi*57.3 >>s1dc_fdc0_phi(400,-20,20,400,-20,20)",
               "triggers[5]==1", "colz");
    c1.Print("out/theta_phi_comparison.pdf", "pdf");

    c1.Clear();
    c1.Divide(1);

    c1.cd(1);
    chain.Draw("s1dc_theta*57.3:fdc0_theta*57.3 >>s1dc_fdc0_theta2(400,-100,100,400,-100,100)",
               "triggers[5]==1", "colz");
    c1.Print("out/theta_phi_comparison.pdf", "pdf");

    c1.Clear();
    c1.Divide(1,2);

    c1.cd(1);
    chain.Draw("s1dc_theta*57.3:fdc0_theta*57.3 >>s1dc_fdc_theta2(400,-20,20,400,-20,20)",
                   "triggers[5]==1","colz");
    TH1* s1_fdc0_good_projx = static_cast<TH2*>(gPad->GetPrimitive("s1dc_fdc_theta2"))->ProjectionX("s1_fdc0_good_projx", 1, 400);

    c1.cd(2);
    chain.Draw("s1dc_theta*57.3:fdc0_theta*57.3 >>s1dc_fdc_theta3(400,-20,20,400,85,95)",
               "triggers[5]==1","colz");
    TH1* s1_bad_projx = static_cast<TH2*>(gPad->GetPrimitive("s1dc_fdc_theta3"))->ProjectionX("s1_bad_projx", 1, 400);
    c1.Print("out/theta_phi_comparison.pdf", "pdf");

    c1.Clear();
    c1.Divide(1,3);

    c1.cd(1);
    s1_fdc0_good_projx->SetTitle("S1 - good tracks distribution");
    s1_fdc0_good_projx->GetYaxis()->SetTitle("Events number");
    s1_fdc0_good_projx->GetXaxis()->SetTitle("#theta [lab. deg.]");
    s1_fdc0_good_projx->Draw();

    c1.cd(2);
    s1_bad_projx->SetTitle("S1 - bad tracks distribution");
    s1_bad_projx->GetYaxis()->SetTitle("Events number");
    s1_bad_projx->GetXaxis()->SetTitle("#theta [lab. deg.]");
    s1_bad_projx->Draw();

    c1.cd(3);
    auto s1_total_tracks = (TH1*)s1_fdc0_good_projx->Clone("s1_total");
    s1_total_tracks->Add(s1_bad_projx);
    s1_total_tracks->SetTitle("S1 - total tracks distribution");
    s1_total_tracks->GetYaxis()->SetTitle("Events number");
    s1_total_tracks->GetXaxis()->SetTitle("#theta [lab. deg.]");
    s1_total_tracks->Draw();
    c1.Print("out/theta_phi_comparison.pdf", "pdf");

    c1.Clear();
    c1.Divide(1);
    c1.cd(1);
    auto s1_eff = (TH1*)s1_fdc0_good_projx->Clone("s1_eff");
    s1_eff->Divide(s1_total_tracks);
    s1_eff->SetTitle("S1 efficiency");
    s1_eff->GetYaxis()->SetTitle("Tracking efficiency [%]");
    s1_eff->GetXaxis()->SetTitle("#theta [lab. deg.]");
    s1_eff->Draw();
    c1.Print("out/theta_phi_comparison.pdf", "pdf");

    c1.Clear();
    c1.Divide(1,2);

    c1.cd(1);
    chain.Draw("s1dc_theta*57.3:fdc0_theta*57.3 >>s1dc_fdc_theta4(400,-20,20,400,-20,20)",
               "triggers[5]==1","colz");
    TH1* s1_fdc0_good_projy = static_cast<TH2*>(gPad->GetPrimitive("s1dc_fdc_theta4"))->ProjectionY("s1_fdc0_good_projy", 1, 400);

    c1.cd(2);
    chain.Draw("s1dc_theta*57.3:fdc0_theta*57.3 >>s1dc_fdc_theta5(400,85,95,400,-20,20)",
               "triggers[5]==1","colz");
    TH1* fdc0_bad_projy = static_cast<TH2*>(gPad->GetPrimitive("s1dc_fdc_theta5"))->ProjectionY("fdc0_bad_projx", 1, 400);
    c1.Print("out/theta_phi_comparison.pdf", "pdf");

    c1.Clear();
    c1.Divide(1,3);

    c1.cd(1);
    s1_fdc0_good_projy->SetTitle("FDC0 - good tracks distribution");
    s1_fdc0_good_projy->GetYaxis()->SetTitle("Events number");
    s1_fdc0_good_projy->GetXaxis()->SetTitle("#theta [lab. deg.]");
    s1_fdc0_good_projy->Draw();

    c1.cd(2);
    fdc0_bad_projy->SetTitle("FDC0 - bad tracks distribution");
    fdc0_bad_projy->GetYaxis()->SetTitle("Events number");
    fdc0_bad_projy->GetXaxis()->SetTitle("#theta [lab. deg.]");
    fdc0_bad_projy->Draw();

    c1.cd(3);
    auto fdc0_total_tracks = (TH1*)s1_fdc0_good_projy->Clone("fdc0_total");
    fdc0_total_tracks->Add(fdc0_bad_projy);
    fdc0_total_tracks->SetTitle("FDC0 - total tracks distribution");
    fdc0_total_tracks->GetYaxis()->SetTitle("Events number");
    fdc0_total_tracks->GetXaxis()->SetTitle("#theta [lab. deg.]");
    fdc0_total_tracks->Draw();
    c1.Print("out/theta_phi_comparison.pdf", "pdf");

    c1.Clear();
    c1.Divide(1);
    c1.cd(1);
    auto fdc0_eff = (TH1*)s1_fdc0_good_projy->Clone("fdc0_eff");
    fdc0_eff->Divide(fdc0_total_tracks);
    fdc0_eff->SetTitle("FDC0 efficiency");
    fdc0_eff->GetYaxis()->SetTitle("Tracking efficiency [%]");
    fdc0_eff->GetXaxis()->SetTitle("#theta [lab. deg.]");
    fdc0_eff->Draw();
    c1.Print("out/theta_phi_comparison.pdf)", "pdf");


    f_hist_out.Write();
    f_hist_out.Close();
}
