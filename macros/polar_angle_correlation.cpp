#include <iostream>

#include "TChain.h"


void polar_angle_correlation() {
    TChain chain{"scattree"};

    //for (int i = 272; i < 283; i++){
    for (int i = 272; i < 283; i++){
        TString name;
        name.Form("scattree_rootfiles/scattree_run%i.root", i);
        std::cout << "Adding file to chain: " << name << std::endl;
        chain.Add(name);
    }

    TCutG* esr_nai2_gcut = (TCutG*)
        (new TFile("cuts/esr_nai2_cut.root"))->Get("CUTG");
    esr_nai2_gcut->SetName("esr_nai2_gcut");
    TCutG* esl_nai2_gcut = (TCutG*)
        (new TFile("cuts/esl_nai2_cut.root"))->Get("CUTG");
    esl_nai2_gcut->SetName("esl_nai2_gcut");

    TCutG* esl_xpos_aang_gcut = (TCutG*)
        (new TFile("cuts/esl_rdc_xpos_aang_cut.root"))->Get("CUTG");
    esl_xpos_aang_gcut->SetName("esl_rdc_xpos_aang_gcut");

    TFile hist_out("out/polar_angle_correlation.root", "RECREATE");
    TCanvas c1("c1");
    gStyle->SetOptStat(1111111);

    c1.Divide(2);

    c1.cd(1);
    chain.Draw(
        "p_theta*57.3:fdc0_theta*57.3 >>h1(200,0,20,200,50,75)",
        "triggers[5]==1 && esl_rdc_xpos_aang_gcut",
        "colz");
    c1.cd(2);
    chain.Draw("p_theta*57.3:fdc0_theta*57.3 >>h2(200,0,20,200,50,75)",
               "triggers[5]==1", "colz");
    c1.Print("out/polar_angle_correlation.pdf", "pdf");


    hist_out.Write();
    hist_out.Close();
}
