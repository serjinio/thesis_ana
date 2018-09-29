#include <iostream>

#include "TChain.h"


void hdraw(TTree& tree, TString name, TString draw_cmd,
           TString binning, TCut cuts = "", TString title = "",
           TString xaxis_title = "", TString yaxis_title = "",
           TString draw_opts = "colz") {
    TString hstr = TString::Format("%s >>%s%s", draw_cmd.Data(), name.Data(),
                                   binning.Data());
    tree.Draw(hstr, cuts, draw_opts);
    TH1 *hist = static_cast<TH1*>(gPad->GetPrimitive(name));
    if (yaxis_title != "") {
        hist->GetYaxis()->SetTitle(yaxis_title);
    }
    if (xaxis_title != "") {
        hist->GetXaxis()->SetTitle(xaxis_title);
    }
    if (title != "") {
        hist->SetTitle(title);
    }
}


void proton_e_de() {
    TChain chain{"scattree"};

    //for (int i = 272; i < 283; i++){
    for (int i = 272; i < 297; i++){
        TString name;
        name.Form("scattree_rootfiles/scattree_run%i.root", i);
        std::cout << "Adding file to chain: " << name << std::endl;
        chain.Add(name);
    }


    TCut esl_phi_corr_cut = "abs(-1.3 * s1dc_phi*57.3 - 25 - esl_p_phi*57.3) < 15";
    TCut esr_phi_corr_cut = "abs(-1.35 * s1dc_phi*57.3 + 33 - esr_p_phi*57.3) < 15";

    TCut target_cut{"tgt_up_xpos*tgt_up_xpos + tgt_up_ypos*tgt_up_ypos < 64"};
    TCut vertex_zpos_cut{"es_vertex_zpos > -40 && es_vertex_zpos < 40"};

    TCut rdc_cut = "proton_theta*57.3<70 && proton_theta*57.3>55";
    TCut esl_rdc_cut = "esl_xpos > -220 && esl_xpos < 160";
    TCut esr_rdc_cut = "esr_xpos > -160 && esr_xpos < 220";

    TFile hist_out("out/proton_pid.root", "RECREATE");
    TCanvas c1("c1");
    gStyle->SetOptStat(1111111);

    c1.Divide(3,3);

    for (int i = 0 ; i < 7; i++) {
        c1.cd(i+1);
        TCut nai_min_e_cut = TString::Format("esl_naie[%i]>100", i).Data();
        auto draw_cmd = TString::Format("esl_naie[%i]:esl_de", i);
        hdraw(chain, TString::Format("hnail%i", i), draw_cmd,
              "(100,0,2000,100,0,4000)",
              "triggers[5]==1" && nai_min_e_cut && target_cut &&
              vertex_zpos_cut && esl_phi_corr_cut,
              TString::Format("NaI #%i - dE-E", i), "dE [ch]", "E [ch]", "colz");
        gPad->SetLogz();
    }
    c1.Print("out/proton_pid.pdf(", "pdf");

    c1.Clear();
    c1.Divide(3,3);

    for (int i = 0 ; i < 7; i++) {
        c1.cd(i+1);
        TCut nai_min_e_cut = TString::Format("esr_naie[%i]>100", i).Data();
        auto draw_cmd = TString::Format("esr_naie[%i]:esr_de", i);
        hdraw(chain, TString::Format("hnail%i", i), draw_cmd,
              "(100,0,2000,100,0,4000)",
              "triggers[5]==1" && nai_min_e_cut && target_cut &&
              vertex_zpos_cut && esr_phi_corr_cut,
              TString::Format("NaI #%i - dE-E", i), "dE [ch]", "E [ch]", "colz");
        gPad->SetLogz();
    }
    c1.Print("out/proton_pid.pdf)", "pdf");

    // c1.Clear();
    // c1.Divide(2);
    // c1.cd(1);
    // chain.Draw("esl_ypos:esl_xpos >>h_rdc1(400,-250,250,400,-250,250)",
    //            "triggers[5]==1", "colz");
    // c1.cd(2);
    // chain.Draw("esr_ypos:esr_xpos >>h_rdc2(400,-250,250,400,-250,250)",
    //            "triggers[5]==1", "colz");
    // c1.Print("out/proton_pid.pdf)", "pdf");

    hist_out.Write();
    hist_out.Close();
}
