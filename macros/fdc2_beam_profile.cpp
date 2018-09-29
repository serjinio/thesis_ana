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

void fdc2_beam_profile() {
  init_dataset();

  TString output_root_filename = TString::Format("out/%s.root", __FILENAME__);
  TFile hist_out(output_root_filename, "RECREATE");
  TCanvas c1("c1");
  gROOT->SetStyle("Video");
  gStyle->SetOptStat(nullptr);

  ElasticScatteringCuts cuts;
  cuts.apply_vertexXY_cut = true;
  cuts.apply_theta_cut = false;
  cuts.apply_phi_cut = false;
  cuts.apply_vertexZ_cut = false;
  cuts.espri_selector = Espri::both;
  cuts.apply_hodf_he6_cut=  false;

  c1.Clear(); c1.Divide(1,1);
  c1.cd(1);
  hdraw(g_chain_total, "fdc2_profile", "fdc2_ypos:-fdc2_xpos",
        "(200,-1200,1200,200,-400,400)",
        cuts.GetTotalCut(),
        "Fragments disctribution on FDC2",
        "X [mm]", "Y [mm]", "col");
  gPad->SetLogz();
  print_pdf_first(c1);

  c1.Clear(); c1.Divide(1,1);
  c1.cd(1);
  hdraw(g_chain_total, "fdc2_profile", "-fdc2_xpos",
        "(200,-1200,1200,200,-400,400)",
        cuts.GetTotalCut() && "fdc2_ypos > -400 && fdc2_ypos < 400",
        "Fragments disctribution on FDC2 - X projection",
        "X [mm]", "Counts", "col");
  gPad->SetLogy();
  print_pdf(c1);

  c1.Clear(); c1.Divide(1,1);
  //gStyle->SetOptStat(1111);
  c1.cd(1);
  hdraw(g_chain_total, "hodf_profile", "hodf_q:Iteration$",
        "(24,0,23,400,0,40)",
        cuts.GetTotalCut(),
        "Fragmens distribution on HODF",
        "Pla. ID", "dE [MeV]", "col");
  gPad->SetLogz();
  // draw_marker_line(-12, 0, -12, 140e3, 3);
  // draw_marker_line(12, 0, 12, 140e3, 3);
  print_pdf_last(c1);


  hist_out.Write();
  hist_out.Close();
}
