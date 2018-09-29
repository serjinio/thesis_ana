#include <iostream>

#include "common.hpp"
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

void on_target_beam_profile() {
  init_dataset();

  TString output_root_filename = TString::Format("out/%s.root", __FILENAME__);
  TFile hist_out(output_root_filename, "RECREATE");
  TCanvas c1("c1");
  gROOT->SetStyle("Video");
  gStyle->SetOptStat(nullptr);

  ElasticScatteringCuts cuts;
  cuts.apply_vertexXY_cut = false;
  cuts.apply_theta_cut = false;
  cuts.apply_phi_cut = false;
  cuts.apply_vertexZ_cut = false;
  cuts.espri_selector = Espri::both;

  c1.Clear(); c1.Divide(1,1);
  c1.cd(1);
  hdraw(g_chain_total, "tgt_up_profile", "tgt_up_ypos:tgt_up_xpos",
        "(400,-40,40,400,-40,40)",
        "triggers[1]==1",
        "Beam on target",
        "X [mm]", "Y [mm]", "col");
  gPad->SetLogz();
  draw_marker_circle(0,0,7,5);
  print_pdf_first(c1);

  c1.Clear(); c1.Divide(1,1);
  gStyle->SetOptStat(1111);
  c1.cd(1);
  hdraw(g_chain_total, "tgt_up_xpos", "tgt_up_xpos",
        "(80,-40,40)",
        "triggers[1]==1 && tgt_up_ypos > -40 && tgt_up_ypos < 40",
        "Beam on target - X projection",
        "X [mm]", "Y [mm]");
  gPad->SetLogz();
  draw_marker_line(-7, 0, -7, 390e3, 3);
  draw_marker_line(7, 0, 7, 390e3, 3);
  print_pdf(c1);

  c1.Clear(); c1.Divide(1,1);
  c1.cd(1);
  hdraw(g_chain_total, "tgt_up_ypos", "tgt_up_ypos",
        "(80,-40,40)",
        "triggers[1]==1 && tgt_up_xpos > -40 && tgt_up_xpos < 40",
        "Beam on target - Y projection",
        "X [mm]", "Y [mm]");
  gPad->SetLogz();
  draw_marker_line(-7, 0, -7, 300e3, 3);
  draw_marker_line(7, 0, 7, 300e3, 3);
  print_pdf_last(c1);


  hist_out.Write();
  hist_out.Close();
}
