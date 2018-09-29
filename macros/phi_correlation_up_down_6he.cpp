#include <iostream>

#include "init_6he_ds.hpp"
#include "elacuts.hpp"
#include "TChain.h"

#include <string.h>

#define __FILENAME__ (strrchr(__FILE__, '/') ? strrchr(__FILE__, '/') + 1 : __FILE__)

using namespace s13::ana;

void inline print_pdf_first(TCanvas& canv) {
  TString output_filename = TString::Format("out/%s.pdf(", __FILENAME__);
  canv.Print(output_filename, "pdf");
}

void inline print_pdf(TCanvas& canv) {
  TString output_filename = TString::Format("out/%s.pdf", __FILENAME__);
  canv.Print(output_filename, "pdf");
}

void inline print_pdf_last(TCanvas& canv) {
  TString output_filename = TString::Format("out/%s.pdf)", __FILENAME__);
  canv.Print(output_filename, "pdf");
}

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


void draw_theta_corr(TCanvas& canv, TChain& chain, ElasticScatteringCuts& cuts,
                     TString title, bool draw2d = true) {
  TCut p_angle_acceptance_cut = "p_theta*57.3>55 && p_theta*57.3<70";
  
  if (draw2d) {
    hdraw(chain, cuts.GetTotalCutName(),
          "p_theta*57.3:s1dc_theta*57.3 - he_theta_theor*57.3",
          "(200,-8,8,200,50,75)",
          cuts.GetTotalCut(),
          title + ", " + cuts.GetTotalCutTitle(),
          "Fragment: #theta_{exp} - #theta_{theor} [lab. deg.]",
          "Proton #theta angle [lab. deg.]", "colz");
  } else {
    hdraw(chain, cuts.GetTotalCutName()+"_1d",
          "s1dc_theta*57.3 - he_theta_theor*57.3",
          "(200,-8,8)",
          cuts.GetTotalCut() && p_angle_acceptance_cut,
          title + " 1D, " + cuts.GetTotalCutTitle(),
          "Fragment: #theta_{exp} - #theta_{theor} [lab. deg.]",
          "Counts", "colz");
  }
}

void phi_correlation_up_down_6he() {
  init_dataset();
  init_carbon_dataset();

  TString output_root_filename = TString::Format("out/%s.root", __FILENAME__);
  TFile hist_out(output_root_filename, "RECREATE");
  TCanvas c1("c1");
  gStyle->SetOptStat(1111111);

  ElasticScatteringCuts cuts;
  cuts.apply_vertexXY_cut = true;
  cuts.apply_theta_cut = false;
  cuts.apply_phi_cut = true;
  cuts.apply_hodf_he6_cut = true;
  cuts.phi_width = 16;
  cuts.vertexXY_radius = 9;
  cuts.espri_selector = Espri::both;

  // c1.Clear();
  // c1.Divide(2,2);
  // c1.cd(1);
  // draw_theta_corr(c1, Espri::right, g_chain_up,
  //                 "#theta_{exp}-#theta_{theor} pol. up");
  // c1.cd(3);
  // draw_theta_corr(c1, Espri::right, g_chain_up,
  //                 "#theta_{exp}-#theta_{theor} pol. up", false);
  // c1.cd(2);
  // draw_theta_corr(c1, Espri::left, g_chain_up,
  //                 "#theta_{exp}-#theta_{theor} pol. up");
  // c1.cd(4);
  // draw_theta_corr(c1, Espri::left, g_chain_up,
  //                 "#theta_{exp}-#theta_{theor} pol. up", false);
  // c1.Print("out/phi_correlation_up_down_6he.pdf(", "pdf");

  // c1.Clear();
  // c1.Divide(2,2);
  // c1.cd(1);
  // draw_theta_corr(c1, Espri::right, g_chain_down,
  //                 "#theta_{exp}-#theta_{theor} pol. down");
  // c1.cd(3);
  // draw_theta_corr(c1, Espri::right, g_chain_down,
  //                 "#theta_{exp}-#theta_{theor} pol. down", false);
  // c1.cd(2);
  // draw_theta_corr(c1, Espri::left, g_chain_down,
  //                 "#theta_{exp}-#theta_{theor} pol. down");
  // c1.cd(4);
  // draw_theta_corr(c1, Espri::left, g_chain_down,
  //                 "#theta_{exp}-#theta_{theor} pol. down", false);
  // c1.Print("out/phi_correlation_up_down_6he.pdf", "pdf");


  // c1.Clear();
  // c1.Divide(1,1);
  // c1.cd(1);
  // draw_theta_corr(c1, Espri::both, g_chain_total,
  //                 "#theta_{exp}-#theta_{theor} pol. up & down");
  // c1.Print("out/phi_correlation_up_down_6he.pdf", "pdf");

  // c1.Clear();
  // c1.Divide(1,1);
  // c1.cd(1);
  // draw_theta_corr(c1, Espri::both, g_chain_total,
  //                 "#theta_{exp}-#theta_{theor}  up & down", false);
  // c1.Print("out/phi_correlation_up_down_6he.pdf)", "pdf");


  c1.Clear();
  c1.Divide(4,2);
  c1.cd(1);
  cuts.apply_phi_cut = false;
  cuts.apply_hodf_he6_cut = false;
  draw_theta_corr(c1, g_c_chain, cuts,
                  "#theta_{exp}-#theta_{theor} pol. up & down");
  c1.cd(5);
  draw_theta_corr(c1, g_c_chain, cuts,
                  "#theta_{exp}-#theta_{theor} pol. up & down", false);
  cuts.apply_phi_cut = true;
  c1.cd(2);
  draw_theta_corr(c1, g_c_chain, cuts,
                  "#theta_{exp}-#theta_{theor} pol. up & down");
  c1.cd(6);
  draw_theta_corr(c1, g_c_chain, cuts,
                  "#theta_{exp}-#theta_{theor} pol. up & down", false);
  cuts.apply_phi_cut = false;
  cuts.apply_hodf_he6_cut = true;
  c1.cd(3);
  draw_theta_corr(c1, g_c_chain, cuts,
                  "#theta_{exp}-#theta_{theor} pol. up & down");
  c1.cd(7);
  draw_theta_corr(c1, g_c_chain, cuts,
                  "#theta_{exp}-#theta_{theor} pol. up & down", false);
  cuts.apply_phi_cut = true;
  c1.cd(4);
  draw_theta_corr(c1, g_c_chain, cuts,
                  "#theta_{exp}-#theta_{theor} pol. up & down");
  c1.cd(8);
  draw_theta_corr(c1, g_c_chain, cuts,
                  "#theta_{exp}-#theta_{theor} pol. up & down", false);
  print_pdf(c1);


  hist_out.Write();
  hist_out.Close();
}
