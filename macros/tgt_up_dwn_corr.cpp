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


void tgt_up_dwn_corr() {
    TChain chain{"scattree"};

    //for (int i = 272; i < 283; i++){
    for (int i = 268; i < 269; i++) {
        TString name;
        name.Form("scattree_rootfiles/scattree_run%i.root", i);
        std::cout << "Adding file to chain: " << name << std::endl;
        chain.Add(name);
    }

    TCutG* tgt_up_dwn_xpos_gcut = (TCutG*)
        (new TFile("cuts/tgt_up_dwn_xpos_gcut.root"))->Get("CUTG");
    tgt_up_dwn_xpos_gcut->SetName("tgt_up_dwn_xpos_gcut");
    TCutG* tgt_up_dwn_ypos_gcut = (TCutG*)
        (new TFile("cuts/tgt_up_dwn_ypos_gcut.root"))->Get("CUTG");
    tgt_up_dwn_ypos_gcut->SetName("tgt_up_dwn_ypos_gcut");

    TCut tgt_up_dwn_xpos_cut = "abs((tgt_up_xpos - tgt_dwn_xpos) - -3.32) < 7";
    TCut tgt_up_dwn_ypos_cut = "abs((tgt_up_ypos - tgt_dwn_ypos) - 1.13e-2) < 7";

    TCut phi_corr_cut_1d = "abs(abs(p_phi*57.3-s1dc_phi*57.3) - 180) < 15";

    TCut target_cut{"tgt_up_xpos*tgt_up_xpos + tgt_up_ypos*tgt_up_ypos < 100"};
    TCut vertex_zpos_cut{"es_vertex_zpos > -45 && es_vertex_zpos < 45"};

    TFile hist_out("out/tgt_up_dwn_corr.root", "RECREATE");
    TCanvas c1("c1");
    gStyle->SetOptStat(1111111);


    // beam profiles on DCs
    c1.Clear();
    c1.Divide(2,2);
    c1.cd(1);
    hdraw(chain, "bdc1_pos", "bdc1_ypos:bdc1_xpos",
          "(400,-40,40,400,-40,40)",
          "triggers[1]==1",
          "BDC1 X - Y positions",
          "Xpos [mm]", "Ypos [mm]");
    gPad->SetLogz();
    c1.cd(2);
    hdraw(chain, "bdc2_pos", "bdc2_ypos:bdc2_xpos",
          "(400,-40,40,400,-40,40)",
          "triggers[1]==1",
          "BDC2 X - Y positions",
          "Xpos [mm]", "Ypos [mm]");
    gPad->SetLogz();
    c1.cd(3);
    hdraw(chain, "fdc0_pos", "fdc0_ypos:fdc0_xpos",
          "(400,-40,40,400,-40,40)",
          "triggers[1]==1",
          "FDC0 X - Y positions",
          "Xpos [mm]", "Ypos [mm]");
    gPad->SetLogz();
    c1.cd(4);
    hdraw(chain, "s1dc_pos", "s1dc_ypos:s1dc_xpos",
          "(400,-40,40,400,-40,40)",
          "triggers[1]==1",
          "S1DC X - Y positions",
          "Xpos [mm]", "Ypos [mm]");
    gPad->SetLogz();
    c1.Print("out/tgt_up_dwn_corr.pdf(", "pdf");

    // beam profiles on DCs - extrapolations from up & down
    c1.Clear();
    c1.Divide(2);
    c1.cd(1);
    hdraw(chain, "tgt_up_pos", "tgt_up_ypos:tgt_up_xpos",
          "(400,-40,40,400,-40,40)",
          "triggers[1]==1",
          "Target up X - Y positions extrapolation",
          "Xpos [mm]", "Ypos [mm]");
    gPad->SetLogz();
    c1.cd(2);
    hdraw(chain, "tgt_dwn_pos", "tgt_dwn_ypos:tgt_dwn_xpos",
          "(400,-40,40,400,-40,40)",
          "triggers[1]==1",
          "Target down X - Y positions extrapolation",
          "Xpos [mm]", "Ypos [mm]");
    gPad->SetLogz();
    c1.Print("out/tgt_up_dwn_corr.pdf", "pdf");

    // BDC1&2 correlations
    c1.Clear();
    c1.Divide(2);
    c1.cd(1);
    hdraw(chain, "bdc1_bdc2_xpos", "bdc1_xpos:bdc2_xpos",
          "(400,-40,40,400,-40,40)",
          "triggers[1]==1",
          "BDC1 - BDC2 Xpos correlation",
          "BDC2 xpos [mm]", "BDC1 xpos [mm]");
    gPad->SetLogz();
    c1.cd(2);
    hdraw(chain, "bdc1_bdc2_ypos", "bdc1_ypos:bdc2_ypos",
          "(400,-40,40,400,-40,40)",
          "triggers[1]==1",
          "BDC1 - BDC2 Ypos correlation",
          "BDC2 ypos [mm]", "BDC1 ypos [mm]");
    gPad->SetLogz();
    c1.Print("out/tgt_up_dwn_corr.pdf", "pdf");

    // on target correlations
    c1.Clear();
    c1.Divide(2,2);
    c1.cd(1);
    hdraw(chain, "tgt_up_dwn_xpos", "tgt_up_xpos:tgt_dwn_xpos",
          "(400,-40,40,400,-40,40)",
          "triggers[1]==1",
          "Target up - down xpos correlation",
          "Target down xpos [mm]", "Target up xpos [mm]");
    gPad->SetLogz();
    c1.cd(2);
    hdraw(chain, "tgt_up_dwn_ypos", "tgt_up_ypos:tgt_dwn_ypos",
          "(400,-40,40,400,-40,40)",
          "triggers[1]==1",
          "Target up - down ypos correlation",
          "Target down ypos [mm]", "Target up ypos [mm]");
    gPad->SetLogz();
    c1.cd(3);
    hdraw(chain, "tgt_up_dwn_aang", "tgt_up_aang*1e3:tgt_dwn_aang*1e3",
          "(400,-40,40,400,-40,40)",
          "triggers[1]==1",
          "Target up - down aang correlation",
          "Target down aang [mrad]", "Target up aang [mrad]");
    gPad->SetLogz();
    c1.cd(4);
    hdraw(chain, "tgt_up_dwn_bang", "tgt_up_bang*1e3:tgt_dwn_bang*1e3",
          "(400,-40,40,400,-40,40)",
          "triggers[1]==1",
          "Target up - down bang correlation",
          "Target down bang [mrad]", "Target up (BDCs) bang [mrad]");
    gPad->SetLogz();
    c1.Print("out/tgt_up_dwn_corr.pdf", "pdf");

    // on target correlations - 1d
    c1.Clear();
    c1.Divide(2,2);
    c1.cd(1);
    hdraw(chain, "tgt_up_dwn_xpos_diff", "tgt_up_xpos - tgt_dwn_xpos",
          "(400,-40,40)",
          "triggers[1]==1 && abs(tgt_up_xpos - tgt_dwn_xpos)>0.01",
          "Target up - down xpos difference",
          "xpos difference [mm]", "Counts");
    gPad->SetLogy();
    c1.cd(2);
    hdraw(chain, "tgt_up_dwn_ypos_diff", "tgt_up_ypos - tgt_dwn_ypos",
          "(400,-40,40)",
          "triggers[1]==1 && abs(tgt_up_ypos - tgt_dwn_ypos)>0.01",
          "Target up - down ypos difference",
          "ypos difference [mm]", "Counts");
    gPad->SetLogy();
    c1.cd(3);
    hdraw(chain, "tgt_up_dwn_aang_diff", "tgt_up_aang*1e3 - tgt_dwn_aang*1e3",
          "(400,-40,40)",
          "triggers[1]==1 && abs(tgt_up_aang*1e3 - tgt_dwn_aang*1e3)>0.01",
          "Target up - down aang difference",
          "aang difference [mrad]", "Counts");
    gPad->SetLogy();
    c1.cd(4);
    hdraw(chain, "tgt_up_dwn_bang_diff", "tgt_up_bang*1e3 - tgt_dwn_bang*1e3",
          "(400,-40,40)",
          "triggers[1]==1 && abs(tgt_up_bang*1e3 - tgt_dwn_bang*1e3)>0.01",
          "Target up - down bang difference",
          "bang difference [mrad]", "Counts");
    gPad->SetLogy();
    c1.Print("out/tgt_up_dwn_corr.pdf", "pdf");

    // on target correlations - 1d (with cuts)
    c1.Clear();
    c1.Divide(2,2);
    c1.cd(1);
    hdraw(chain, "tgt_up_dwn_xpos_diff", "tgt_up_xpos - tgt_dwn_xpos",
          "(400,-40,40)",
          "triggers[1]==1 && abs(tgt_up_xpos - tgt_dwn_xpos)>0.01" &&
          tgt_up_dwn_xpos_cut && tgt_up_dwn_ypos_cut,
          "Target up - down xpos difference with cut",
          "xpos difference [mm]", "Counts");
    gPad->SetLogy();
    c1.cd(2);
    hdraw(chain, "tgt_up_dwn_ypos_diff", "tgt_up_ypos - tgt_dwn_ypos",
          "(400,-40,40)",
          "triggers[1]==1 && abs(tgt_up_ypos - tgt_dwn_ypos)>0.01" &&
          tgt_up_dwn_xpos_cut && tgt_up_dwn_ypos_cut,
          "Target up - down ypos difference with cut",
          "ypos difference [mm]", "Counts");
    gPad->SetLogy();
    c1.cd(3);
    hdraw(chain, "tgt_up_dwn_aang_diff", "tgt_up_aang*1e3 - tgt_dwn_aang*1e3",
          "(400,-40,40)",
          "triggers[1]==1 && abs(tgt_up_aang*1e3 - tgt_dwn_aang*1e3)>0.01" &&
          tgt_up_dwn_xpos_cut && tgt_up_dwn_ypos_cut,
          "Target up - down aang difference with cut",
          "aang difference [mrad]", "Counts");
    gPad->SetLogy();
    c1.cd(4);
    hdraw(chain, "tgt_up_dwn_bang_diff", "tgt_up_bang*1e3 - tgt_dwn_bang*1e3",
          "(400,-40,40)",
          "triggers[1]==1 && abs(tgt_up_bang*1e3 - tgt_dwn_bang*1e3)>0.01" &&
          tgt_up_dwn_xpos_cut && tgt_up_dwn_ypos_cut,
          "Target up - down bang difference with cut",
          "bang difference [mrad]", "Counts");
    gPad->SetLogy();
    c1.Print("out/tgt_up_dwn_corr.pdf", "pdf");

    // BDC2 - FDC0 correlations
    c1.Clear();
    c1.Divide(2,2);
    c1.cd(1);
    hdraw(chain, "bdc2_fdc_xpos", "bdc2_xpos:fdc0_xpos",
          "(400,-40,40,400,-40,40)",
          "triggers[1]==1",
          "BDC2 - FDC0 X pos correlation",
          "FDC0 xpos [mm]", "BDC2 xpos [mm]");
    gPad->SetLogz();
    c1.cd(2);
    hdraw(chain, "bdc2_fdc_ypos", "bdc2_ypos:fdc0_ypos",
          "(400,-40,40,400,-40,40)",
          "triggers[1]==1",
          "BDC2 - FDC0 Ypos correlation",
          "FDC0 ypos [mm]", "BDC2 ypos [mm]");
    gPad->SetLogz();
    c1.cd(3);
    hdraw(chain, "bdc2_fdc_aang", "tgt_up_aang*1e3:fdc0_aang*1e3",
          "(400,-80,80,400,-80,80)",
          "triggers[1]==1",
          "BDCs - FDC0 aang correlation",
          "FDC0 aang [mrad]", "Target up (BDCs) aang [mrad]");
    gPad->SetLogz();
    c1.cd(4);
    hdraw(chain, "bdc2_fdc_bang", "tgt_up_bang*1e3:fdc0_bang*1e3",
          "(400,-80,80,400,-80,80)",
          "triggers[1]==1",
          "BDCs - FDC0 bang correlation",
          "FDC0 bang [mrad]", "Target up (BDCs) bang [mrad]");
    gPad->SetLogz();
    c1.Print("out/tgt_up_dwn_corr.pdf", "pdf");

    // FDC0 - S1DC correlations
    c1.Clear();
    c1.Divide(2,2);
    c1.cd(1);
    hdraw(chain, "fdc0_s1dc_xpos", "fdc0_xpos:s1dc_xpos",
          "(400,-40,40,400,-40,40)",
          "triggers[1]==1",
          "FDC0 - S1DC X pos correlation",
          "S1DC xpos [mm]", "FDC0 xpos [mm]");
    gPad->SetLogz();
    c1.cd(2);
    hdraw(chain, "fdc0_s1dc_ypos", "fdc0_ypos:s1dc_ypos",
          "(400,-40,40,400,-40,40)",
          "triggers[1]==1",
          "FDC0 - S1DC Ypos correlation",
          "S1DC ypos [mm]", "FDC0 ypos [mm]");
    gPad->SetLogz();
    c1.cd(3);
    hdraw(chain, "fdc0_s1dc_aang", "fdc0_aang*1e3:s1dc_aang*1e3",
          "(400,-80,80,400,-80,80)",
          "triggers[1]==1",
          "FDC0 - S1DC aang correlation",
          "S1DC aang [mrad]", "FDC0 aang [mrad]");
    gPad->SetLogz();
    c1.cd(4);
    hdraw(chain, "fdc0_s1dc_bang", "fdc0_bang*1e3:s1dc_bang*1e3",
          "(400,-80,80,400,-80,80)",
          "triggers[1]==1",
          "FDC0 - S1DC bang correlation",
          "S1DC bang [mrad]", "FDC0 bang [mrad]");
    gPad->SetLogz();
    c1.Print("out/tgt_up_dwn_corr.pdf", "pdf");

    // BDC2 - S1DC correlations
    c1.Clear();
    c1.Divide(2,2);
    c1.cd(1);
    hdraw(chain, "bdc2_fdc_xpos", "bdc2_xpos:s1dc_xpos",
          "(400,-40,40,400,-40,40)",
          "triggers[1]==1",
          "BDC2 - S1DC Xpos correlation",
          "S1DC xpos [mm]", "BDC2 xpos [mm]");
    gPad->SetLogz();
    c1.cd(2);
    hdraw(chain, "bdc2_fdc_ypos", "bdc2_ypos:s1dc_ypos",
          "(400,-40,40,400,-40,40)",
          "triggers[1]==1",
          "BDC2 - S1DC Ypos correlation",
          "S1DC ypos [mm]", "BDC2 ypos [mm]");
    gPad->SetLogz();
    c1.cd(3);
    hdraw(chain, "bdc2_fdc_aang", "tgt_up_aang*1e3:s1dc_aang*1e3",
          "(400,-80,80,400,-80,80)",
          "triggers[1]==1",
          "BDCs - S1DC aang correlation",
          "S1DC aang [mrad]", "Target up (BDCs) aang [mrad]");
    gPad->SetLogz();
    c1.cd(4);
    hdraw(chain, "bdc2_fdc_bang", "tgt_up_bang*1e3:s1dc_bang*1e3",
          "(400,-80,80,400,-80,80)",
          "triggers[1]==1",
          "BDCs - S1DC bang correlation",
          "S1DC bang [mrad]", "Target up (BDCs) bang [mrad]");
    gPad->SetLogz();
    c1.Print("out/tgt_up_dwn_corr.pdf", "pdf");

    // FDC0 - FDC2 correlations
    c1.Clear();
    c1.Divide(2,2);
    c1.cd(1);
    hdraw(chain, "fdc0_fdc2_xpos", "fdc0_xpos:fdc2_xpos",
          "(400,0,600,400,-80,80)",
          "triggers[1]==1",
          "FDC0 - FDC2 X pos correlation",
          "FDC2 xpos [mm]", "FDC0 xpos [mm]");
    gPad->SetLogz();
    c1.cd(2);
    hdraw(chain, "fdc0_fdc2_ypos", "fdc0_ypos:fdc2_ypos",
          "(400,-600,600,400,-80,80)",
          "triggers[1]==1",
          "FDC0 - FDC2 Ypos correlation",
          "FDC2 ypos [mm]", "FDC0 ypos [mm]");
    gPad->SetLogz();
    c1.cd(3);
    hdraw(chain, "fdc0_fdc2_aang", "fdc0_aang*1e3:fdc2_aang*1e3",
          "(400,-800,800,400,-80,80)",
          "triggers[1]==1",
          "FDC0 - FDC2 aang correlation",
          "FDC2 aang [mrad]", "FDC0 aang [mrad]");
    gPad->SetLogz();
    c1.cd(4);
    hdraw(chain, "fdc0_fdc2_bang", "fdc0_bang*1e3:fdc2_bang*1e3",
          "(200,-200,200,200,-80,80)",
          "triggers[1]==1",
          "FDC0 - FDC2 bang correlation",
          "FDC2 bang [mrad]", "FDC0 bang [mrad]");
    gPad->SetLogz();
    c1.Print("out/tgt_up_dwn_corr.pdf", "pdf");

    // S1DC - FDC2 correlations
    c1.Clear();
    c1.Divide(2,2);
    c1.cd(1);
    hdraw(chain, "s1dc_fdc2_xpos", "s1dc_xpos:fdc2_xpos",
          "(400,0,600,400,-80,80)",
          "triggers[1]==1",
          "S1DC - FDC2 X pos correlation",
          "FDC2 xpos [mm]", "S1DC xpos [mm]");
    gPad->SetLogz();
    c1.cd(2);
    hdraw(chain, "s1dc_fdc2_ypos", "s1dc_ypos:fdc2_ypos",
          "(400,-600,600,400,-80,80)",
          "triggers[1]==1",
          "S1DC - FDC2 Ypos correlation",
          "FDC2 ypos [mm]", "S1DC ypos [mm]");
    gPad->SetLogz();
    c1.cd(3);
    hdraw(chain, "s1dc_fdc2_aang", "s1dc_aang*1e3:fdc2_aang*1e3",
          "(400,-800,800,400,-80,80)",
          "triggers[1]==1",
          "S1DC - FDC2 aang correlation",
          "FDC2 aang [mrad]", "S1DC aang [mrad]");
    gPad->SetLogz();
    c1.cd(4);
    hdraw(chain, "s1dc_fdc2_bang", "s1dc_bang*1e3:fdc2_bang*1e3",
          "(200,-200,200,200,-80,80)",
          "triggers[1]==1",
          "S1DC - FDC2 bang correlation",
          "FDC2 bang [mrad]", "S1DC bang [mrad]");
    gPad->SetLogz();
    c1.Print("out/tgt_up_dwn_corr.pdf)", "pdf");


    hist_out.Write();
    hist_out.Close();
}
