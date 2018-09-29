#include <iostream>

#include "init_ds.hpp"
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

/*
  Function returns a graph of 6He-p kinematics
*/
TGraph* BuildHe6KinematicsGraph() {
  static const double p_lab_theta[] =
    {50, 52, 53, 55., 57., 59., 61., 64, 67, 70, 72, 73, 75};
  static const double he_lab_theta[] =
    {8.5, 8.2, 8.05, 7.72, 7.37, 7.0, 6.62, 6.01, 5.37, 4.72, 4.25, 4.03, 3.57};
  static const int num_points = 13;

  TGraph* graph = new TGraph(num_points, he_lab_theta, p_lab_theta);
  graph->SetLineColor(kRed);
  graph->SetLineWidth(5);
  graph->SetMarkerColor(5);
  graph->SetMarkerStyle(22);
  graph->SetTitle("6he-p kinematics");
  graph->GetXaxis()->SetTitle("He #theta [deg]");
  graph->GetYaxis()->SetTitle("Proton #theta [deg]");
  graph->SetMinimum(50);
  graph->SetMaximum(75);

  return graph;
};

void draw_theta_corr(TCanvas& canv, TChain& chain, ElasticScatteringCuts& cuts,
                     TString title, TString draw_opts = "colz", bool draw2d = true,
                     int bin_num = 200, TString name = "") {
  TCut p_angle_acceptance_cut = "p_theta*57.3>55 && p_theta*57.3<70";
  TString hist_name = name == "" ? cuts.GetTotalCutName() : name;

  if (draw2d) {
    TString limits = TString::Format("(%i,0,10,%i,50,75)", bin_num, bin_num);
    hdraw(chain, hist_name,
          "p_theta*57.3:s1dc_theta*57.3",
          limits,
          cuts.GetTotalCut(),
          //title + ", " + cuts.GetTotalCutTitle(),
          title,
          "#theta_{6He} [lab. deg.]",
          "#theta_{p} [lab. deg.]", draw_opts);
  } else {
    if (name == "") {
      hist_name += "_1d";
    }
    TString limits = TString::Format("(%i,-8,8)", bin_num);
    hdraw(chain, hist_name,
          "s1dc_theta*57.3 - he_theta_theor*57.3",
          limits,
          cuts.GetTotalCut() && p_angle_acceptance_cut,
          //title + " 1D, " + cuts.GetTotalCutTitle(),
          title + " - 1D",
          "Fragment: #theta_{exp} - #theta_{theor} [lab. deg.]",
          "Counts", draw_opts);
  }
}

TH1* get_scaled_carbon_theta_corr_1d(ElasticScatteringCuts& cuts) {
  const double k_carbon_data_scaling_factor = 11.89;
  TCanvas* canv = new TCanvas("c_tmp_carbon");

  canv->Clear();
  canv->Divide(1,1);
  canv->cd(1);
  draw_theta_corr(*canv, g_chain_carbon, cuts,
                  "#theta correlation (6He - C)", "",
                  false, 90, "tmp_carbon_1d");
  TH1* hist = static_cast<TH1*>(gPad->GetPrimitive("tmp_carbon_1d"));
  std::cout << "got carbon hist: " << hist << std::endl << std::flush;
  hist->Scale(k_carbon_data_scaling_factor);
  return hist;
}

TH1* get_polt_theta_corr_1d(ElasticScatteringCuts& cuts) {
  TCanvas* canv = new TCanvas("c_tmp_polt");

  canv->Clear();
  canv->Divide(1,1);
  canv->cd(1);
  draw_theta_corr(*canv, g_chain_total, cuts,
                  "#theta correlation (6He - p), pol. up & down", "",
                  false, 90, "tmp_polt_1d");
  TH1* hist = static_cast<TH1*>(gPad->GetPrimitive("tmp_polt_1d"));
  return hist;
}

void theta_corr_6he_carbon() {
  init_6he_dataset();
  init_carbon_dataset();

  TString output_root_filename = TString::Format("out/%s.root", __FILENAME__);
  TFile hist_out(output_root_filename, "RECREATE");
  TCanvas c1("c1");
  gROOT->SetStyle("Video");
  gStyle->SetLineWidth(2);
  gStyle->SetOptStat(00000000);

  ElasticScatteringCuts cuts;
  cuts.apply_vertexXY_cut = true;
  cuts.apply_theta_cut = false;
  cuts.apply_phi_cut = true;
  cuts.apply_vertexZ_cut = false;
  cuts.apply_espri_aang_cut = true;
  cuts.apply_espri_e_cut = false;
  cuts.apply_hodf_he6_cut = true;

  cuts.phi_width = 20;
  cuts.espri_aang_width = 180;
  cuts.espri_e_width = 0.9;
  cuts.vertexXY_radius = 10;
  cuts.espri_selector = Espri::both;

  c1.Clear(); c1.Divide(1,1);
  c1.cd(1);
  draw_theta_corr(c1, g_chain_total, cuts,
                  "Kinematical correlation (6He-p)", "col",
                  true, 180);
  auto kin_graph = BuildHe6KinematicsGraph();
  c1.cd(1);
  kin_graph->Draw("C same");
  print_pdf(c1);

  // c1.Clear(); c1.Divide(1,1);
  // c1.cd(1);
  // draw_theta_corr(c1, g_chain_total, cuts,
  //                 "Kinematical correlation (6He-p)", "col",
  //                 false, 90);
  // print_pdf(c1);

  // c1.Clear(); c1.Divide(1,1);
  // c1.cd(1);
  // draw_theta_corr(c1, g_chain_carbon, cuts,
  //                 "Kinematical correlation (6He-C)", "col",
  //                 true, 90);
  // print_pdf(c1);

  // c1.Clear(); c1.Divide(1,1);
  // c1.cd(1);
  // draw_theta_corr(c1, g_chain_carbon, cuts,
  //                 "Kinematical correlation (6He-C)", "col",
  //                 false, 90);
  // print_pdf(c1);

  // TH1* carbon_hist = get_scaled_carbon_theta_corr_1d(cuts);
  // TH1* polt_hist = get_polt_theta_corr_1d(cuts);
  // c1.Clear(); c1.Divide(1,1);
  // c1.cd(1);
  // polt_hist->SetTitle("6He-p, 6He-C (scaled)");
  // polt_hist->GetYaxis()->SetRangeUser(0, 1800);
  // polt_hist->Draw();
  // carbon_hist->SetLineColor(kRed);
  // carbon_hist->Draw("same");
  // print_pdf_first(c1);

  // polt_hist->Add(carbon_hist, -1);
  // polt_hist->GetYaxis()->SetRangeUser(0, 1800);
  // c1.Clear(); c1.Divide(1,1);
  // c1.cd(1);
  // polt_hist->SetTitle("6He-p (QFS carbon BG subtracted)");
  // polt_hist->Draw();
  // print_pdf_last(c1);


  hist_out.Write();
  hist_out.Close();
}
