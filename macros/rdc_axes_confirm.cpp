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


void rdc_axes_confirm() {

  TChain chain{"scattree"};

  //for (int i = 272; i < 283; i++){
  for (int i = 272; i < 280; i++) {
    TString name;
    name.Form("scattree_rootfiles/scattree_run%i.root", i);
    std::cout << "Adding file to chain: " << name << std::endl;
    chain.Add(name);
  }

  TCut phi_corr_cut_1d = "abs(abs(p_phi*57.3-s1dc_phi*57.3) - 180) < 10";
  TCut target_cut{"tgt_up_xpos*tgt_up_xpos + tgt_up_ypos*tgt_up_ypos < 36"};
  TCut vertex_zpos_cut{"es_vertex_zpos > -45 && es_vertex_zpos < 30"};

  TFile hist_out("out/tgt_up_dwn_corr.root", "RECREATE");
  TCanvas c1("c1");
  gStyle->SetOptStat(1111);
  
  c1.Clear();
  c1.Divide(2,2);
  c1.cd(1);
  hdraw(chain, "rdcl_yaxis_up", "esl_ypos:esl_xpos", "(200,-220,220,200,-220,220)",
        "triggers[5]==1 && esl_naie[0]>100", "RDC left and NaI # 1",
        "RDC pos X [mm]", "RDC pos Y [mm]");
  gPad->SetLogz();
  c1.cd(2);
  hdraw(chain, "rdcr_yaxis_up", "esr_ypos:esr_xpos", "(200,-220,220,200,-220,220)",
        "triggers[5]==1 && esr_naie[0]>100", "RDC right and NaI # 1",
        "RDC pos X [mm]", "RDC pos Y [mm]");
  gPad->SetLogz();
  c1.cd(3);
  hdraw(chain, "rdcl_yaxis_down", "esl_ypos:esl_xpos", "(200,-220,220,200,-220,220)",
        "triggers[5]==1 && esl_naie[6]>100", "RDC left and NaI # 7",
        "RDC pos X [mm]", "RDC pos Y [mm]");
  gPad->SetLogz();
  c1.cd(4);
  hdraw(chain, "rdcr_yaxis_down", "esr_ypos:esr_xpos", "(200,-220,220,200,-220,220)",
        "triggers[5]==1 && esr_naie[6]>100", "RDC right and NaI # 7",
        "RDC pos X [mm]", "RDC pos Y [mm]");
  gPad->SetLogz();
  
  c1.Print("out/rdc_axes_confirm.pdf(", "pdf");

  // S1DC - RDC
  c1.Clear();
  c1.Divide(2,2);
  c1.cd(1);
  hdraw(chain, "rdcl_xaxis_corr", "esl_xpos:s1dc_xpos", "(200,-250,250,200,-220,220)",
        "triggers[5]==1" && target_cut && vertex_zpos_cut && phi_corr_cut_1d, 
        "RDC left X - S1DC X correlation", "S1DC X [mm]", "RDC X [mm]");
  c1.cd(3);
  hdraw(chain, "rdcr_xaxis_corr", "esr_xpos:s1dc_xpos", "(200,-250,250,200,-220,220)",
        "triggers[5]==1" && target_cut && vertex_zpos_cut && phi_corr_cut_1d, 
        "RDC right X - S1DC X correlation", "S1DC X [mm]", "RDC X [mm]");
  c1.cd(2);
  hdraw(chain, "rdcl_yaxis_corr", "esl_ypos:s1dc_ypos", "(200,-65,65,200,-220,220)",
        "triggers[5]==1" && target_cut && vertex_zpos_cut &&
        phi_corr_cut_1d,
        "RDC left Y - S1DC Y correlation", "S1DC Y [mm]", "RDC Y [mm]");
  c1.cd(4);
  hdraw(chain, "rdcr_yaxis_corr", "esr_ypos:s1dc_ypos", "(200,-65,65,200,-220,220)",
        "triggers[5]==1" && target_cut && vertex_zpos_cut && 
        phi_corr_cut_1d,
        "RDC right Y - S1DC Y correlation", "S1DC Y [mm]", "RDC Y [mm]");
  c1.Print("out/rdc_axes_confirm.pdf", "pdf");

  // FDC0 - RDC
  c1.Clear();
  c1.Divide(2,2);
  c1.cd(1);
  hdraw(chain, "rdcl_xaxis_corr", "esl_xpos:fdc0_xpos", "(200,-80,80,200,-220,220)",
        "triggers[5]==1" && target_cut && vertex_zpos_cut && phi_corr_cut_1d,
        "RDC left X - FDC0 X correlation", "FDC0 X [mm]", "RDC X [mm]");
  c1.cd(3);
  hdraw(chain, "rdcr_xaxis_corr", "esr_xpos:fdc0_xpos", "(200,-80,80,200,-220,220)",
        "triggers[5]==1" && target_cut && vertex_zpos_cut && phi_corr_cut_1d,
        "RDC right X - FDC0 X correlation", "FDC0 X [mm]", "RDC X [mm]");
  c1.cd(2);
  hdraw(chain, "rdcl_yaxis_corr", "esl_ypos:fdc0_ypos", "(200,-40,40,200,-220,220)",
        "triggers[5]==1" && target_cut && vertex_zpos_cut &&
        phi_corr_cut_1d, 
        "RDC left Y - FDC0 Y correlation", "FDC0 Y [mm]", "RDC Y [mm]");
  c1.cd(4);
  hdraw(chain, "rdcr_yaxis_corr", "esr_ypos:fdc0_ypos", "(200,-40,40,200,-220,220)",
        "triggers[5]==1" && target_cut && vertex_zpos_cut && 
        phi_corr_cut_1d,
        "RDC right Y - FDC0 Y correlation", "FDC0 Y [mm]", "RDC Y [mm]");
  c1.Print("out/rdc_axes_confirm.pdf)", "pdf");

  hist_out.Write();
  hist_out.Close();
}
