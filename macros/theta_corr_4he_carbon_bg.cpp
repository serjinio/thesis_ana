#include <iostream>

#include "common.hpp"
#include "init_4he_ds.hpp"
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

/*
  Function returns a graph of 6He-p kinematics 
*/
TGraph* BuildHe4KinematicsGraph() {
  static const double p_lab_theta[] =
    {50, 52, 53, 55., 57., 59., 61., 64, 67, 70, 72, 73, 75};
  static const double he_lab_theta[] = 
    {12.42, 11.95, 11.71, 11.19, 10.65, 10.1, 9.82, 9.52, 
     8.6, 7.69, 6.71, 6.07, 5.75, 5.08};
  static const int num_points = 13;

  TGraph* graph = new TGraph(num_points, he_lab_theta, p_lab_theta);
  graph->SetLineColor(kRed);
  graph->SetLineWidth(5);
  graph->SetMarkerColor(5);
  graph->SetMarkerStyle(22);
  graph->SetTitle("4he-p kinematics");
  graph->GetXaxis()->SetTitle("He #theta [deg]");
  graph->GetYaxis()->SetTitle("Proton #theta [deg]");
  graph->SetMinimum(50);
  graph->SetMaximum(75);

  return graph;
};

void draw_theta_corr(TChain& chain, TCut cuts,
                     TString title, TString draw_opts = "colz", bool draw2d = true,
                     int bin_num = 100, TString name = "") {
  static int hist_num = 1;

  TCut p_angle_acceptance_cut = "p_theta*57.3>55 && p_theta*57.3<70";
  TString hist_name = name == "" ? TString::Format("thcorr_%i", hist_num) : name;
  
  if (draw2d) {
    TString limits = TString::Format("(%i,0,15,%i,50,75)", bin_num, bin_num);
    hdraw(chain, hist_name,
          "p_theta*57.3:s1dc_theta*57.3",
          limits,
          cuts,
          title,
          "Fragment #theta [lab. deg.]",
          "Proton #theta [lab. deg.]", draw_opts);
  } else {
    if (name == "") {
      hist_name += "_1d";
    }
    TString limits = TString::Format("(%i,-8,8)", bin_num);
    hdraw(chain, hist_name,
          "s1dc_theta*57.3 - he_theta_theor*57.3",
          limits,
          cuts && p_angle_acceptance_cut,
          title + " 1D",
          "Fragment: #theta_{exp} - #theta_{theor} [lab. deg.]",
          "Counts", draw_opts);
  }

  hist_num += 1;
}

TH1* get_thcorr_hist(TCut cuts, TString title) {
  static int hist_num = 101;
  TCanvas canv("c_tmp_hists");
  TString hist_name = TString::Format("thcorr_1d_%i", hist_num);

  canv.Clear(); canv.Divide(1,1);
  canv.cd(1);
  draw_theta_corr(g_chain_total, cuts,
                  title, "", false, 100, hist_name);
  TH1* hist = static_cast<TH1*>(gPad->GetPrimitive(hist_name));
  hist = (TH1*)hist->Clone();

  hist_num += 1;
  return hist;
}

void draw_phi_corr(TChain& chain,
                   TString draw_opts = "colz",
                   float phi_width = 0, float phi_width2 = 0) {
  static int hist_num = 1;
  ElasticScatteringCuts cuts;
  cuts.apply_vertexXY_cut = true;
  cuts.apply_vertexZ_cut = false;
  cuts.apply_phi_cut = false;
  cuts.apply_theta_cut = false;

  hdraw(chain, cuts.GetTotalCutName() + "_" + std::to_string(hist_num),
        "abs(p_phi*57.3-s1dc_phi*57.3)",
        "(300,120,240)",
        cuts.GetTotalCut(),
        "abs(p_{#phi} - He_{#phi}) (" + cuts.GetTotalCutTitle() + ")",
        "p_{#phi} - He_{#phi} [deg]",
        "Counts");
  if (phi_width != 0) {
    draw_marker_line(177.5 - -phi_width/2., 0,
                     177.5 - -phi_width/2., 400);
    draw_marker_line(177.5 - phi_width/2., 0,
                     177.5 - phi_width/2., 400);
  }
  if (phi_width2 != 0) {
    draw_marker_line(177.5 - -phi_width2/2., 0,
                     177.5 - -phi_width2/2., 400);
    draw_marker_line(177.5 - phi_width2/2., 0,
                     177.5 - phi_width2/2., 400);
  }
  hist_num += 1;
}

void draw_phi_cut(TCanvas& canv, ElasticScatteringCuts& cuts,
                  TChain& chain, float signal_phi_thresh,
                  float noise_phi_min_thresh, float noise_phi_max_thresh) {
  TCut signal_phi_cut =
    TString::Format("abs(abs(p_phi*57.3-s1dc_phi*57.3) - 177.5) < %.2f",
                    signal_phi_thresh).Data();
  TCut noise_phi_min_cut =
    TString::Format("abs(abs(p_phi*57.3-s1dc_phi*57.3) - 177.5) > %.2f",
                    noise_phi_min_thresh).Data();
  TCut noise_phi_max_cut =
    TString::Format("abs(abs(p_phi*57.3-s1dc_phi*57.3) - 177.5) < %.2f",
                    noise_phi_max_thresh).Data();
  TCut noise_phi_cut = noise_phi_min_cut && noise_phi_max_cut;

  canv.Clear(); canv.Divide(3,2);
  int canv_pad_counter = 1;

  TString phi_cut_title = TString::Format("4He-p: #phi cut (|#phi|<%.1f)", signal_phi_thresh);
  TString not_phi_cut_title = TString::Format("4He-p: #phi cut (%.1f<|#phi|<%.1f)",
                                              noise_phi_min_thresh, noise_phi_max_thresh);

  canv.cd(1);
  draw_phi_corr(g_chain_total, "", signal_phi_thresh);
  canv.cd(2);
  draw_theta_corr(g_chain_total, cuts.GetTotalCut() && signal_phi_cut,
                  phi_cut_title, "colz", true);
  canv.cd(3);
  draw_theta_corr(g_chain_total, cuts.GetTotalCut() && signal_phi_cut,
                  phi_cut_title + " 1D", "", false);

  canv.cd(4);
  draw_phi_corr(g_chain_total, "", noise_phi_min_thresh, noise_phi_max_thresh);
  canv.cd(5);
  draw_theta_corr(g_chain_total, cuts.GetTotalCut() && noise_phi_cut,
                  not_phi_cut_title, "colz", true);
  canv.cd(6);
  draw_theta_corr(g_chain_total, cuts.GetTotalCut() && noise_phi_cut,
                  not_phi_cut_title + " 1D", "", false);
}

void theta_corr_4he_carbon_bg() {
  init_dataset();

  TString output_root_filename = TString::Format("out/%s.root", __FILENAME__);
  TFile hist_out(output_root_filename, "RECREATE");o
  TCanvas c1("c1");
  gStyle->SetOptStat(1111111);

  ElasticScatteringCuts cuts;
  cuts.apply_vertexXY_cut = true;
  cuts.apply_theta_cut = false;
  cuts.apply_phi_cut = false;
  cuts.apply_vertexZ_cut = false;
  cuts.espri_selector = Espri::both;

  float signal_phi_thresh = 16;
  float noise_phi_min_thresh = 24;
  float noise_phi_max_thresh = 48;
  TCut signal_phi_cut =
    TString::Format("abs(abs(p_phi*57.3-s1dc_phi*57.3) - 177.5) < %.2f",
                    signal_phi_thresh).Data();
  TCut noise_phi_min_cut =
    TString::Format("abs(abs(p_phi*57.3-s1dc_phi*57.3) - 177.5) > %.2f",
                    noise_phi_min_thresh).Data();
  TCut noise_phi_max_cut =
    TString::Format("abs(abs(p_phi*57.3-s1dc_phi*57.3) - 177.5) < %.2f",
                    noise_phi_max_thresh).Data();
  TCut noise_phi_cut = noise_phi_min_cut && noise_phi_max_cut;

  draw_phi_cut(c1, cuts, g_chain_total, signal_phi_thresh,
               noise_phi_min_thresh, noise_phi_max_thresh);
  print_pdf_first(c1);

  c1.Clear(); c1.Divide(1,1);
  c1.cd(1);
  TH1* hist_signal = get_thcorr_hist(cuts.GetTotalCut() && signal_phi_cut,
                                     "4He signal events");
  TH1* hist_noise = get_thcorr_hist(cuts.GetTotalCut() && noise_phi_cut,
                                    "4He noise events");
  hist_signal->Draw();
  hist_noise->SetLineColor(kRed);
  hist_noise->Draw("same");
  print_pdf(c1);


  signal_phi_thresh = 16;
  noise_phi_min_thresh = 24;
  noise_phi_max_thresh = 72;
  signal_phi_cut =
    TString::Format("abs(abs(p_phi*57.3-s1dc_phi*57.3) - 177.5) < %.2f",
                    signal_phi_thresh).Data();
  noise_phi_min_cut =
    TString::Format("abs(abs(p_phi*57.3-s1dc_phi*57.3) - 177.5) > %.2f",
                    noise_phi_min_thresh).Data();
  noise_phi_max_cut =
    TString::Format("abs(abs(p_phi*57.3-s1dc_phi*57.3) - 177.5) < %.2f",
                    noise_phi_max_thresh).Data();
  noise_phi_cut = noise_phi_min_cut && noise_phi_max_cut;

  draw_phi_cut(c1, cuts, g_chain_total, signal_phi_thresh,
               noise_phi_min_thresh, noise_phi_max_thresh);
  print_pdf_first(c1);

  c1.Clear(); c1.Divide(1,1);
  c1.cd(1);
  hist_signal = get_thcorr_hist(cuts.GetTotalCut() && signal_phi_cut,
                                "4He signal events");
  hist_noise = get_thcorr_hist(cuts.GetTotalCut() && noise_phi_cut,
                               "4He noise events");
  hist_signal->Draw();
  hist_noise->SetLineColor(kRed);
  hist_noise->Draw("same");
  print_pdf(c1);


  signal_phi_thresh = 16;
  noise_phi_min_thresh = 16;
  noise_phi_max_thresh = 64;
  signal_phi_cut =
    TString::Format("abs(abs(p_phi*57.3-s1dc_phi*57.3) - 177.5) < %.2f",
                    signal_phi_thresh).Data();
  noise_phi_min_cut =
    TString::Format("abs(abs(p_phi*57.3-s1dc_phi*57.3) - 177.5) > %.2f",
                    noise_phi_min_thresh).Data();
  noise_phi_max_cut =
    TString::Format("abs(abs(p_phi*57.3-s1dc_phi*57.3) - 177.5) < %.2f",
                    noise_phi_max_thresh).Data();
  noise_phi_cut = noise_phi_min_cut && noise_phi_max_cut;

  draw_phi_cut(c1, cuts, g_chain_total, signal_phi_thresh,
               noise_phi_min_thresh, noise_phi_max_thresh);
  print_pdf_first(c1);

  c1.Clear(); c1.Divide(1,1);
  c1.cd(1);
  hist_signal = get_thcorr_hist(cuts.GetTotalCut() && signal_phi_cut,
                                "4He signal events");
  hist_noise = get_thcorr_hist(cuts.GetTotalCut() && noise_phi_cut,
                               "4He noise events");
  hist_signal->Draw();
  hist_noise->SetLineColor(kRed);
  hist_noise->Draw("same");
  print_pdf_last(c1);

}
