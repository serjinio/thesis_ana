
#include <iostream>
#include <fstream>

#include "TChain.h"
#include "TGraph.h"

#include "init_4he_ds.hpp"
#include "macroincludes.hpp"


using namespace s13::ana;


const Double_t gk_phi_cut_width = 16;
const Double_t gk_phi_corr_center = 177.5;
const Double_t gk_yaxis_limit = 100;

/*
  Function computes nomalization factor to estimate contribution of background
  to the elastic scattering peak. R is defined as follows:

  R = Y_bg_tot / Y_bg_out

  During assymmetry computation to get yield of background events when
  R_tot_in is known, just:

  Y_bg_in = Y_bg_out * (R - 1)

  Then elastic yield will be given by:

  Y_ela_in = Y_tot_in - Y_bg_in
*/
Double_t compute_bg_norm_factor(FitParams& fit_params_bg, FitParams& fit_params_signal,
                                Double_t phi_cut_width, Double_t phi_corr_center) {
  TF1 bg_fit_fn("bg_fit_fn", "gaus(0)", 120, 240);
  bg_fit_fn.SetParameters(&fit_params_bg.params[0]);
  TF1 signal_fit_fn("signal_fit_fn", "gaus(0)", 120, 240);
  signal_fit_fn.SetParameters(&fit_params_signal.params[0]);

  Double_t a_lim = phi_corr_center - phi_cut_width / 2.;
  Double_t b_lim = phi_corr_center + phi_cut_width / 2.;
  Double_t bg_total_int = bg_fit_fn.Integral(120, 240);
  Double_t bg_in_int = bg_fit_fn.Integral(a_lim, b_lim);
  Double_t bg_out_int = bg_total_int - bg_in_int;

  return bg_total_int / bg_out_int;
}


void compute_total_to_bg_norm_factor() {
  //init_dataset();
  RootScript script("total_to_bg_norm_factor");

  std::ifstream bg_fit_params_str, signal_fit_params_str;
  bg_fit_params_str.open("data/phi_corr_bg_fit_params.csv");
  signal_fit_params_str.open("data/phi_corr_signal_fit_params.csv");
  std::ofstream r_norm_str;
  r_norm_str.open("data/total_to_bg_norm_factor_2.csv");

  script.NewPage(4,2);
  int pad_num = 1;
  while(bg_fit_params_str) {
    std::string bg_params_line, signal_params_line;
    std::getline(bg_fit_params_str, bg_params_line);
    std::getline(signal_fit_params_str, signal_params_line);

    std::stringstream bg_line_str(bg_params_line);
    std::stringstream signal_line_str(signal_params_line);
    std::string bg_cell, signal_cell;

    std::getline(bg_line_str, bg_cell, ',');
    if (bg_cell == "") {
      std::cout << "WARN: Empty line in input, skipping!" << std::endl;
      continue;
    }

    Double_t p_theta = std::stod(bg_cell);
    // ignore same value from another file
    std::getline(signal_line_str, signal_cell, ',');

    std::cout << "Drawing fit functions for mean proton theta: " << p_theta <<
      "..." << std::endl;

    FitParams bg_fit_params, signal_fit_params;
    bg_line_str >> bg_fit_params;
    signal_line_str >> signal_fit_params;
    cout << "BG fit params: " << bg_fit_params;
    cout << "Signal fit params: " << signal_fit_params;

    auto bg_fit_fn = new TF1("bg_fit_fn", "gaus(0)", 120, 240);
    bg_fit_fn->SetParameters(&bg_fit_params.params[0]);
    auto signal_fit_fn = new TF1("signal_fit_fn", "gaus(0)", 120, 240);
    signal_fit_fn->SetParameters(&signal_fit_params.params[0]);

    script.cd(pad_num);
    signal_fit_fn->SetFillColor(kGreen);
    signal_fit_fn->SetLineColor(kGreen);
    TString title = TString::Format("p_{#theta} = %.1f;p_{#phi} - He_{#phi};Counts",
                                    p_theta);
    signal_fit_fn->SetTitle(title);

    bg_fit_fn->SetFillColor(kRed);
    bg_fit_fn->SetLineColor(kRed);
    bg_fit_fn->SetTitle(title);

    signal_fit_fn->Draw();
    bg_fit_fn->Draw("CF SAME");

    // compute and write normalization factor
    Double_t r_norm = compute_bg_norm_factor(bg_fit_params, signal_fit_params,
                                             gk_phi_cut_width, gk_phi_corr_center);
    r_norm_str << p_theta << ",\t" << r_norm << std::endl;
    auto norm_factor_text = new TText(195,90,
                                      TString::Format("R  = %.2f", r_norm));
    norm_factor_text->SetTextSize(0.05);
    norm_factor_text->Draw();

    bg_fit_fn->GetYaxis()->SetRangeUser(0, gk_yaxis_limit);
    signal_fit_fn->GetYaxis()->SetRangeUser(0, gk_yaxis_limit);
    gPad->Update();

    draw_marker_line(gk_phi_corr_center - gk_phi_cut_width/2, 0,
                     gk_phi_corr_center - gk_phi_cut_width/2, gk_yaxis_limit);
    draw_marker_line(gk_phi_corr_center + gk_phi_cut_width/2, 0,
                     gk_phi_corr_center + gk_phi_cut_width/2, gk_yaxis_limit);

    pad_num++;
  }
}
