#include <iostream>
#include <string>

#include "TChain.h"


void hdraw(TTree& tree, TString name, TString draw_cmd,
           TString binning, TCut cuts = "", TString title = "",
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


void Ep_vs_theta() {
    TChain chain{"scattree"};

    //for (int i = 272; i < 283; i++){
    for (int i = 272; i < 283; i++){
        TString name;
        name.Form("scattree_rootfiles/scattree_run%i.root", i);
        std::cout << "Adding file to chain: " << name << std::endl;
        chain.Add(name);
    }


    TCut esl_phi_corr_cut = "abs(-1.3 * s1dc_phi*57.3 - 25 - esl_p_phi*57.3) < 15";
    TCut esr_phi_corr_cut = "abs(-1.35 * s1dc_phi*57.3 + 33 - esr_p_phi*57.3) < 15";

    TCut target_cut{"tgt_up_xpos*tgt_up_xpos + tgt_up_ypos*tgt_up_ypos < 64"};
    TCut vertex_zpos_cut{"es_vertex_zpos > -40 && es_vertex_zpos < 40"};

    TCut esl_vertex_zpos_cut{"esl_vertex_zpos > -7-6 && esl_vertex_zpos < -7+6"};
    TCut esr_vertex_zpos_cut{"esr_vertex_zpos > -9-9 && esr_vertex_zpos < -9+9"};

    TCutG* esr_nai2_gcut = (TCutG*)
        (new TFile("cuts/esr_nai2_cut.root"))->Get("CUTG");
    esr_nai2_gcut->SetName("esr_nai2_gcut");
    TCutG* esl_nai2_gcut = (TCutG*)
        (new TFile("cuts/esl_nai2_cut.root"))->Get("CUTG");
    esl_nai2_gcut->SetName("esl_nai2_gcut");

    TCut polar_angle_4he_corr_cut = "abs(s1dc_theta*57.3 + p_theta*57.3 * 0.25 - 25) < 1";

    TFile hist_out("out/Ep_vs_theta.root", "RECREATE");
    TCanvas c1("c1");
    TH1 *hist;
    gStyle->SetOptStat(1111111);

    c1.Divide(2,4);

    for (int i = 0; i < 7; i++) {
        TString name = TString::Format("naiL%i", i);
        TString draw_cmd = TString::Format("esl_naie[%i]:p_theta*57.3", i);
        c1.cd(i+1);
        hdraw(chain, name, draw_cmd, "(100,50,75,60,100,4000)",
              "triggers[5]==1" && target_cut && polar_angle_4he_corr_cut &&
              esl_vertex_zpos_cut && esl_phi_corr_cut,
              "ESPRI left - Proton E vs. #theta [mm]",
              "Proton #theta [lab. deg.]", "Proton E [ch]");
        gPad->SetLogz();
    }
    c1.Print("out/Ep_vs_theta.pdf(", "pdf");

    for (int i = 0; i < 7; i++) {
        TString name = TString::Format("naiR_%i", i);
        TString draw_cmd = TString::Format("esr_naie[%i]:p_theta*57.3", i);
        c1.cd(i+1);
        hdraw(chain, name, draw_cmd, "(100,50,75,60,100,4000)",
              "triggers[5]==1" && target_cut && polar_angle_4he_corr_cut &&
              esr_vertex_zpos_cut && esr_phi_corr_cut,
              "ESPRI right - Proton E vs. #theta [mm]",
              "Proton #theta [lab. deg.]", "Proton E [ch]");
        gPad->SetLogz();
    }
    c1.Print("out/Ep_vs_theta.pdf", "pdf");

    c1.Clear();
    c1.Divide(1);

    c1.cd(1);
    chain.Draw(
        "p_theta*57.3:s1dc_theta*57.3 >>polar(400,0,15,400,55,70)",
        "triggers[5]==1" && target_cut &&
        (esr_vertex_zpos_cut || esl_vertex_zpos_cut),
        "colz");
    hist = static_cast<TH1*>(gPad->GetPrimitive("polar"));
    hist->GetYaxis()->SetTitle("Proton angle [lab. deg]");
    hist->GetXaxis()->SetTitle("Fragment angle [lab. deg]");
    hist->SetTitle("Polar angle correlation");
    c1.Print("out/Ep_vs_theta.pdf)", "pdf");


    hist_out.Write();
    hist_out.Close();
}
