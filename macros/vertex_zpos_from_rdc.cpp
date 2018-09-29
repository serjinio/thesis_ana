#include <iostream>
#include <string>

#include "TChain.h"


void hdraw(TTree& tree, TString name, TString draw_cmd,
           TString binning, TString cuts = "", TString title = "",
           TString xaxis_title = "", TString yaxis_title = "",
           TString draw_opts = "colz") {
    TString hstr = TString::Format("%s >>%s%s", draw_cmd.Data(), name.Data(),
                                   binning.Data());
    tree.Draw(hstr, cuts, draw_opts);
    TH1 *hist = static_cast<TH1*>(gPad->GetPrimitive(name));
    hist->GetYaxis()->SetTitle(yaxis_title);
    hist->GetXaxis()->SetTitle(xaxis_title);
    hist->SetTitle(title);
}


void vertex_zpos_from_rdc() {
    TChain chain{"scattree"};

    //for (int i = 272; i < 283; i++){
    for (int i = 272; i < 274; i++){
        TString name;
        name.Form("scattree_rootfiles/scattree_run%i.root", i);
        std::cout << "Adding file to chain: " << name << std::endl;
        chain.Add(name);
    }


    TCut target_cut{"tgt_up_xpos*tgt_up_xpos + tgt_up_ypos*tgt_up_ypos < 64"};
    TCut vertex_zpos_cut{"es_vertex_zpos > -40 && es_vertex_zpos < 40"};
    TCut esl_vertex_zpos_cut{"esl_vertex_zpos > -7-13 && esl_vertex_zpos < -7+13"};
    TCut esr_vertex_zpos_cut{"esr_vertex_zpos > -9-18 && esr_vertex_zpos < -9+18"};
    TCutG* esr_nai2_gcut = (TCutG*)
        (new TFile("cuts/esr_nai2_cut.root"))->Get("CUTG");
    esr_nai2_gcut->SetName("esr_nai2_gcut");
    TCutG* esl_nai2_gcut = (TCutG*)
        (new TFile("cuts/esl_nai2_cut.root"))->Get("CUTG");
    esl_nai2_gcut->SetName("esl_nai2_gcut");


    TFile hist_out("out/vertex_zpos_from_rdc.root", "RECREATE");
    TCanvas c1("c1");
    TH1 *hist;
    //gStyle->SetOptTitle(0);
    gStyle->SetOptStat(1111111);

    c1.Divide(2);

    c1.cd(1);
    hdraw(chain, "esl_vert", "esl_vertex_zpos", "(200,-200,200)",
          "triggers[5]==1",
          "ESPRI left vertex position [mm]",
          "Vertex Z position [mm]", "Event count");
    c1.cd(2);
    hdraw(chain, "esr_vert", "esr_vertex_zpos", "(200,-200,200)",
          "triggers[5]==1",
          "ESPRI right vertex position [mm]",
          "Vertex Z position [mm]", "Event count");
    c1.Print("out/vertex_zpos_from_rdc.pdf(", "pdf");

    c1.Clear();
    c1.Divide(1);

    c1.cd(1);
    chain.Draw(
        "esl_p_theta*57.3:s1dc_theta*57.3 >>polar_esl(400,0,20,400,50,75)",
        "triggers[5]==1" && target_cut && vertex_zpos_cut,
        "colz");
    hist = static_cast<TH1*>(gPad->GetPrimitive("polar_esl"));
    hist->GetYaxis()->SetTitle("Proton angle [lab. deg]");
    hist->GetXaxis()->SetTitle("Fragment angle [lab. deg]");
    hist->SetTitle("ESPRI left - Polar angle correlation");

    c1.Print("out/vertex_zpos_from_rdc.pdf", "pdf");

    c1.Clear();
    c1.Divide(1);

    c1.cd(1);
    chain.Draw(
        "esr_p_theta*57.3:s1dc_theta*57.3 >>polar_esr(400,0,20,400,50,75)",
        "triggers[5]==1" && target_cut && vertex_zpos_cut,
        "colz");
    hist = static_cast<TH1*>(gPad->GetPrimitive("polar_esr"));
    hist->GetYaxis()->SetTitle("Proton angle [lab. deg]");
    hist->GetXaxis()->SetTitle("Fragment angle [lab. deg]");
    hist->SetTitle("ESPRI right - Polar angle correlation");

    c1.Print("out/vertex_zpos_from_rdc.pdf", "pdf");

    c1.Clear();
    c1.Divide(1);

    c1.cd(1);
    chain.Draw(
        "p_theta*57.3:s1dc_theta*57.3 >>polar_total_s1dc(400,0,20,400,50,75)",
        "triggers[5]==1" && target_cut && vertex_zpos_cut,
        "colz");
    hist = static_cast<TH1*>(gPad->GetPrimitive("polar_total_s1dc"));
    hist->GetYaxis()->SetTitle("Proton angle [lab. deg]");
    hist->GetXaxis()->SetTitle("Fragment angle [lab. deg]");
    hist->SetTitle("Total - Polar angle correlation (S1DC)");

    c1.Print("out/vertex_zpos_from_rdc.pdf", "pdf");

    c1.cd(1);
    chain.Draw(
        "p_theta*57.3:fdc0_theta*57.3 >>polar_total_fdc0(400,0,20,400,50,75)",
        "triggers[5]==1" && target_cut && vertex_zpos_cut,
        "colz");
    hist = static_cast<TH1*>(gPad->GetPrimitive("polar_total_fdc0"));
    hist->GetYaxis()->SetTitle("Proton angle [lab. deg]");
    hist->GetXaxis()->SetTitle("Fragment angle [lab. deg]");
    hist->SetTitle("Total - Polar angle correlation (FDC0)");

    c1.Print("out/vertex_zpos_from_rdc.pdf)", "pdf");


    hist_out.Write();
    hist_out.Close();
}
