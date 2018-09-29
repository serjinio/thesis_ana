#include <iostream>

#include "TChain.h"


void vertex_pos_from_rdc() {

    TChain chain{"scattree"};

    //for (int i = 272; i < 283; i++){
    for (int i = 272; i < 283; i++){
        TString name;
        name.Form("scattree_rootfiles/scattree_run%i.root", i);
        std::cout << "Adding file to chain: " << name << std::endl;
        chain.Add(name);
    }

    TCutG* esl_xpos_aang_gcut = (TCutG*)
        (new TFile("cuts/esl_rdc_xpos_aang_cut.root"))->Get("CUTG");
    esl_xpos_aang_gcut->SetName("esl_rdc_xpos_aang_gcut");

    TFile hist_out("out/vertex_pos_from_rdc.root", "RECREATE");

    TCanvas c1("c1");
    gStyle->SetOptStat(1111111);

    c1.Divide(2);

    c1.cd(1);
    chain.Draw("esl_aang >>esl_aang(200,-2000,2000)", "triggers[5]==1", "colz");
    c1.cd(2);
    chain.Draw("esr_aang >>esr_aang(200,-2000,2000)", "triggers[5]==1", "colz");
    c1.Print("out/vertex_pos_from_rdc.pdf(", "pdf");

    c1.Clear();
    c1.Divide(2);

    c1.cd(1);
    chain.Draw("esl_aang:esl_xpos >>esl_aang_xpos(400,-220,220,400,-1000,1000)",
               "triggers[5]==1", "colz");
    c1.cd(2);
    TString draw_cmd = "esl_naie[2]:esl_de >>nail2(400,0,2000,400,0,4000)";
    TCut nai_min_e_cut = TString::Format("esl_naie[2]>100").Data();
    chain.Draw(draw_cmd, "triggers[5]==1" && nai_min_e_cut,
               "colz");
    c1.Print("out/vertex_pos_from_rdc.pdf", "pdf");

    c1.Clear();
    c1.Divide(2);

    c1.cd(1);
    chain.Draw("esl_aang:esl_xpos >>esl_aang_xpos_cut1(400,-220,220,400,-1000,1000)",
               "triggers[5]==1 && esl_rdc_xpos_aang_gcut", "colz");
    c1.cd(2);
    draw_cmd = "esl_naie[2]:esl_de >>nail3_cut1(400,0,2000,400,0,4000)";
    chain.Draw(draw_cmd,
               "triggers[5]==1 && esl_rdc_xpos_aang_gcut" && nai_min_e_cut,
               "colz");
    c1.Print("out/vertex_pos_from_rdc.pdf)", "pdf");


    hist_out.Write();
    hist_out.Close();
}
