
#include <iostream>
#include <fstream>

#include "TChain.h"
#include "TGraph.h"

#include "init_ds.hpp"
#include "macroincludes.hpp"


using namespace s13::ana;


TChain *g_local_chain = init_6he_partial_dataset(70);
//TChain *g_local_chain = &g_chain_total;


// phi cut defined as 3 deviations of phi signal gaussian
const Double_t gk_phi_cut_width = 37.8;
// center of the phi_corr: p_theta - He_theta [deg]
// ideally should be 180 [deg]
const Double_t gk_phi_corr_center = 177.5;


TH1* draw_phi_corr(Double_t min_p_theta, Double_t max_p_theta,
                   Espri espri_selector,
                   bool do_cut = false, bool in_cut = true) {
  static int hist_num = 1;
  const int k_yaxis_limit = 250;

  ElasticScatteringCuts cuts;
  cuts.espri_selector = espri_selector;
  cuts.phi_width = gk_phi_cut_width;
  cuts.apply_vertexXY_cut = true;
  cuts.apply_espri_e_cut = false;
  cuts.apply_espri_aang_cut = true;
  cuts.apply_user_cut = true;
  // proton theta range cut
  cuts.user_cut = TString::Format("p_theta*57.3>%.2f && p_theta*57.3<%.2f",
                                  min_p_theta, max_p_theta).Data();
  cuts.apply_theta_cut = false;

  TString title = TString::Format("%.0f<p_{#theta}<%.0f, ",
                                  min_p_theta, max_p_theta);

  if (!do_cut) {
    cuts.apply_phi_cut = false;
  } else {
    cuts.apply_phi_cut = true;
    if (in_cut) {
      title = title + " 'IN' #phi";
    } else {
      title = title + " 'OUT' #phi";
    }
  }

  TCut tcuts = in_cut ? cuts.GetTotalCut() : cuts.GetTotalPhiBgCut();
  TString hist_name = TString::Format("phi_corr_%i", hist_num);
  hdraw(*g_local_chain, hist_name,
        "abs(p_phi*57.3-s1dc_phi*57.3)",
        "(300,120,240)",
        tcuts,
        title + " (" + cuts.GetTotalCutTitle() + ")",
        "p_{#phi} - He_{#phi} [deg]",
        "Counts");

  if (!do_cut) { // if not cutting, then show lines where to cut
    draw_marker_line(gk_phi_corr_center - cuts.phi_width/2, 0,
                     gk_phi_corr_center - cuts.phi_width/2, k_yaxis_limit);
    draw_marker_line(gk_phi_corr_center + cuts.phi_width/2, 0,
                     gk_phi_corr_center + cuts.phi_width/2, k_yaxis_limit);
  }

  TH1* hist = static_cast<TH1*>(gPad->GetPrimitive(hist_name));
  hist->GetYaxis()->SetRangeUser(0, k_yaxis_limit);

  hist_num += 1;
  return hist;
}

FitParams fit_phi_corr(TH1* hist) {
  FitParams fit_params{3};
  const Double_t *fit_params_raw;
  const Double_t *fit_errors_raw;

  if (hist->GetEntries() == 0) {
    std::cout << "WARN: Empty histogram passed for fitting!" << std::endl;
    const Double_t zero_params[] = {0,0,0};
    fit_params.SetParams(zero_params);
    fit_params.SetErrors(zero_params);
    return fit_params;
  }

  hist->Fit("gaus", "W");
  fit_params_raw = hist->GetFunction("gaus")->GetParameters();
  fit_errors_raw = hist->GetFunction("gaus")->GetParErrors();

  fit_params.SetParams(fit_params_raw);
  fit_params.SetErrors(fit_errors_raw);
  return fit_params;
}


/*
  Function computes normalization factor to estimate contribution of background
  to the elastic scattering peak. R is defined as follows:

  R = Y_bg_tot / Y_bg_out

  During asymmetry computation to get yield of background events when
  R_tot_in is known, just:

  Y_bg_in = Y_bg_out * (R - 1)

  Then elastic yield will be given by:

  Y_ela_in = Y_tot_in - Y_bg_in
*/
Double_t compute_bg_norm_factor(FitParams& fit_params_bg, FitParams& fit_params_signal,
                                Double_t phi_cut_width, Double_t phi_corr_center) {
  Double_t phi_max = ElasticScatteringCuts::k_phi_bg_max_angle;
  Double_t phi_min = ElasticScatteringCuts::k_phi_bg_min_angle;

  TF1 bg_fit_fn("bg_fit_fn", "gaus(0)", phi_min, phi_max);
  bg_fit_fn.SetParameters(&fit_params_bg.params[0]);
  TF1 signal_fit_fn("signal_fit_fn", "gaus(0)", phi_min, phi_max);
  signal_fit_fn.SetParameters(&fit_params_signal.params[0]);

  Double_t a_lim = phi_corr_center - phi_cut_width / 2.;
  Double_t b_lim = phi_corr_center + phi_cut_width / 2.;
  Double_t bg_total_int = bg_fit_fn.Integral(phi_min, phi_max);
  Double_t bg_in_int = bg_fit_fn.Integral(a_lim, b_lim);
  Double_t bg_out_int = bg_total_int - bg_in_int;

  return bg_total_int / bg_out_int;
}

void compute_sigmas_range(RootScript& script, Double_t theta_p_start,
                          Espri espri_selector,
                          std::ofstream& bg_fit_params_file,
                          std::ofstream& signal_fit_params_file,
                          std::ofstream& r_norm_file,
                          Double_t theta_p_step = 1.) {
  int pad_num = 1;

  script.NewPage(5,3);
  for (Double_t p_theta = theta_p_start;
       p_theta <= theta_p_start + theta_p_step * 4;
       p_theta += theta_p_step) {
    TH1* hist;
    TH1* hist_overall;
    Double_t p_theta_mean = (p_theta + theta_p_step + p_theta) / 2.0;
    bg_fit_params_file << p_theta_mean << ",";
    signal_fit_params_file << p_theta_mean << ",";
    r_norm_file << p_theta_mean << ",";

    script.cd(pad_num);
    hist_overall = draw_phi_corr(p_theta, p_theta + theta_p_step,
                                 espri_selector, false);

    script.cd(pad_num + 5);
    hist = draw_phi_corr(p_theta, p_theta + theta_p_step,
                         espri_selector, true, true);
    FitParams fit_params_signal = fit_phi_corr(hist);
    signal_fit_params_file << fit_params_signal;

    script.cd(pad_num + 10);
    hist = draw_phi_corr(p_theta, p_theta + theta_p_step,
                         espri_selector, true, false);
    FitParams fit_params_bg = fit_phi_corr(hist);
    bg_fit_params_file << fit_params_bg;

    Double_t fit_params_overall[6];
    std::copy(fit_params_bg.params.begin(), fit_params_bg.params.end(),
              fit_params_overall);
    std::copy(fit_params_signal.params.begin(), fit_params_signal.params.end(),
              fit_params_overall + 3);
    auto total_fit_fn = new TF1("mstotal", "gaus(0)+gaus(3)",
                                ElasticScatteringCuts::k_phi_bg_min_angle,
                                ElasticScatteringCuts::k_phi_bg_max_angle);
    total_fit_fn->SetParameters(fit_params_overall);
    script.cd(pad_num);
    hist_overall->Fit(total_fit_fn, "R");

    // compute and write normalization factor
    Double_t r_norm = compute_bg_norm_factor(fit_params_bg, fit_params_signal,
                                             gk_phi_cut_width, gk_phi_corr_center);
    r_norm_file << r_norm << std::endl;

    pad_num += 1;
    if (p_theta >= 72.) {
      break;
    }
  }
}

void he4_signal_bg_sigmas() {
  RootScript script("he6_signal_bg_sigmas");
  std::ofstream bg_fit_params_file, signal_fit_params_file, r_norm_file;

  bg_fit_params_file.open("data/he6_phi_corr_bg_fit_params.csv");
  signal_fit_params_file.open("data/he6_phi_corr_signal_fit_params.csv");
  r_norm_file.open("data/he6_total_to_bg_norm_factor.csv");

  compute_sigmas_range(script, 53, Espri::both,
                       bg_fit_params_file, signal_fit_params_file, r_norm_file);
  compute_sigmas_range(script, 58, Espri::both,
                       bg_fit_params_file, signal_fit_params_file, r_norm_file);
  compute_sigmas_range(script, 63, Espri::both,
                       bg_fit_params_file, signal_fit_params_file, r_norm_file);
  compute_sigmas_range(script, 68, Espri::both,
                       bg_fit_params_file, signal_fit_params_file, r_norm_file);

  bg_fit_params_file.close();
  signal_fit_params_file.close();
  r_norm_file.close();

//  bg_fit_params_file.open("data/esl_phi_corr_bg_fit_params.csv");
//  signal_fit_params_file.open("data/esl_phi_corr_signal_fit_params.csv");
//  r_norm_file.open("data/esl_total_to_bg_norm_factor.csv");
//
//  compute_sigmas_range(script, 53, Espri::left,
//                       bg_fit_params_file, signal_fit_params_file, r_norm_file);
//  compute_sigmas_range(script, 63, Espri::left,
//                       bg_fit_params_file, signal_fit_params_file, r_norm_file);
//
//  bg_fit_params_file.close();
//  signal_fit_params_file.close();
//  r_norm_file.close();
//
//  bg_fit_params_file.open("data/esr_phi_corr_bg_fit_params.csv");
//  signal_fit_params_file.open("data/esr_phi_corr_signal_fit_params.csv");
//  r_norm_file.open("data/esr_total_to_bg_norm_factor.csv");
//
//  compute_sigmas_range(script, 53, Espri::right,
//                       bg_fit_params_file, signal_fit_params_file, r_norm_file);
//  compute_sigmas_range(script, 63, Espri::right,
//                       bg_fit_params_file, signal_fit_params_file, r_norm_file);
//
//  bg_fit_params_file.close();
//  signal_fit_params_file.close();
//  r_norm_file.close();
}
