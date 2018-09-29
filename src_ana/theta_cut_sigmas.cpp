/**
 *   \file theta_cut_sigmas.cpp
 *   \brief Draws CS given yields as input.
 *
 */

#include <numeric>

#include "TROOT.h"
#include "TLegend.h"

#include "cli.hpp"
#include "init_ds.hpp"
#include "consts.hpp"
#include "csalg.hpp"
#include "beam.hpp"
#include "drawing.hpp"
#include "cuts_conf.hpp"
#include "rootscript.hpp"


using Script = s13::ana::RootScript;
using Yield = s13::ana::ScatteringYield;
using Cuts = s13::ana::ElasticScatteringCuts;
using EspriSel = s13::ana::Espri;
using DsType = s13::ana::DatasetType;
using PolDir = s13::ana::PolarizationDirection;
using Dsdcs = s13::ana::Dsdcs;
using FitParams = s13::ana::FitParams;

static s13::misc::CliOptions cli_opts;


std::pair<std::shared_ptr<s13::ana::ElasticScatteringYieldsAlg>,
          std::shared_ptr<s13::ana::SignalBgMonitorAlg> >
load_yields(DsType ds_type) {
  s13::misc::MessageLogger logger("load_yields()");
  Cuts cuts =
    s13::io::parse_cuts_config(std::fstream(cli_opts.cuts_config()));
  std::string espri_sel = s13::ana::espri_selector_to_str(cuts.espri_selector);
  std::string ds_type_str = "";
  if (ds_type == DsType::he4) {
    ds_type_str = "he4";
  } else if (ds_type == DsType::he6) {
    ds_type_str = "he6";
  } else if (ds_type == DsType::carbon) {
    ds_type_str = "carbon";
  }
  std::string pol_dir_str = "";
  if (cli_opts.pol_direction() == s13::ana::PolarizationDirection::up) {
    pol_dir_str = "_pol_up";
  } else if (cli_opts.pol_direction() == s13::ana::PolarizationDirection::down) {
    pol_dir_str = "_pol_down";
  }

  std::string filename = "cs_data/yields";
  filename += "_" + ds_type_str + "_" + espri_sel + pol_dir_str + ".csv";
  logger.debug("Trying to load ScatteringYield object from file: %s...",
               filename.c_str());
  std::ifstream ifs(filename);
  if (!ifs) {
    throw std::logic_error(s13::ana::tstrfmt("Cannot find specified file: %s",
                                             filename.c_str()).Data());
  }

  logger.debug("Loading scattering yield results...");
  auto p_yields = s13::ana::make_elastic_scattering_yield_alg_from_file(ifs);
  logger.debug("Loading bin monitor results...");
  auto p_sig_bg = s13::ana::make_signal_bg_monitor_alg_from_file(ifs);

  return std::pair<decltype(p_yields), decltype(p_sig_bg)>(p_yields, p_sig_bg);
}

Yield
scale_carbon_bg(Yield carbon_yield) {
  auto scaled_yields = carbon_yield.multiply(s13::gk_carbon_bg_norm_factor);
  return scaled_yields;
}

void
draw_total_kin_corr(Script& script,
                    const std::vector<Yield>& yields_he,
                    const std::vector<Yield>& yields_carbon,
                    TString cut_conditions_str = "",
                    bool scale_carbon = false) {
  using SY = Yield;
  SY total_yield_he = std::accumulate(yields_he.begin(), yields_he.end(),
                                      yields_he[0].Empty());
  SY total_yield_carbon = std::accumulate(yields_carbon.begin(), yields_carbon.end(),
                                          yields_carbon[0].Empty());
  if (scale_carbon) {
    total_yield_carbon = scale_carbon_bg(total_yield_carbon);
  }

  TH1* hist_he =
    total_yield_he.AsHist(-4,4, "Total", "#theta_{{}^{6}He} - #theta_{{}^{6}He-theor}");
  hist_he->SetMaximum(hist_he->GetMaximum() * 1.2);
  TH1* hist_carbon =
    total_yield_carbon.AsHist(-4,4, "Carbon", "#theta_{exp} - #theta_{theor}");
  TH1* hist_he_wo_carbon = static_cast<TH1*>(hist_he->Clone());
  hist_he_wo_carbon->Add(hist_carbon, -1);
  hist_he_wo_carbon->SetTitle("Elastic");

  script.NewPage(1,1).cd();
  s13::ana::draw_thstack({hist_he_wo_carbon, hist_he, hist_carbon},
                         "He & Carbon yields", "NOSTACK HIST");
  gPad->BuildLegend(0.6, 0.65, 0.85, 0.85);
  s13::ana::draw_text(cut_conditions_str, 0.1, 0.8, 0.6, 0.9);
  gPad->Update();
}

TH1* subtract_bin_bg(TH1* signal, TH1* bg, double lab_theta) {
  TH1* signal_wo_bg = static_cast<TH1*>(signal->Clone());
  signal_wo_bg->Add(bg, -1);
  // signal_wo_bg->SetMinimum(0);
  signal_wo_bg->SetName(s13::ana::tstrfmt("%s-nobg", signal->GetName()).Data());
  signal_wo_bg->SetLineColor(kGreen);
  signal_wo_bg->SetMarkerColor(kGreen);
  return signal_wo_bg;
}

FitParams fit_theta_corr(TH1* hist, double min_theta, double max_theta) {
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

  TF1* gaus_fit_fn = new TF1("g1", "gaus", min_theta, max_theta);
  hist->Fit(gaus_fit_fn, "R");
  fit_params_raw = gaus_fit_fn->GetParameters();
  fit_errors_raw = gaus_fit_fn->GetParErrors();

  fit_params.SetParams(fit_params_raw);
  fit_params.SetErrors(fit_errors_raw);
  return fit_params;
}

std::string get_output_filename(const Cuts& cuts) {
  std::string fname = "data/test_theta_sigmas";
  if (cli_opts.dataset_type() == s13::ana::DatasetType::he4) {
    fname += "_he4";
  } else {
    fname += "_he6";
  }
  if (cuts.dsdc_selector == s13::ana::Dsdcs::s1dc) {
    fname += "_s1dc";
  } else if (cuts.dsdc_selector == s13::ana::Dsdcs::fdc0) {
    fname += "_fdc0";
  }
  if (cuts.espri_selector == s13::ana::Espri::left) {
    fname += "_esl";
  } else if (cuts.espri_selector == s13::ana::Espri::right) {
    fname += "_esr";
  }

  fname += ".csv";
  return fname;
}

std::vector<TH1*>
draw_bin_signal_bg_levels(Script& script,
                          DsType ds_type,
                          const Cuts& cuts,
                          const std::vector<Yield>& signals,
                          const std::vector<Yield>& bgs,
                          const Yield& yields) {
  s13::misc::MessageLogger logger("draw_bin_signal_bg_levels()");
  std::vector<TH1*> signals_wo_bg;
  auto bg_norm_fact_interp =
    s13::ana::make_bg_norm_factor_interp2(ds_type,
                                         cuts.dsdc_selector);
  std::ofstream theta_sigmas_ofs;
  theta_sigmas_ofs.open(get_output_filename(cuts));
  theta_sigmas_ofs << "# Gaussian fiting of theta "
                   << "dists. (theta, gaus_const, gaus_mean, "
                   << "gaus_sigma, errs...)" << std::endl;

  logger.debug("Drawing signal & BGs hists...");
  logger.debug("Yields object bins num: %d", yields.GetBinsNumber());
  for (int i = 0; i < yields.GetBinsNumber(); ++i) {
    logger.debug("Drawing for theta %.1f...", yields.GetBinArg(i));
    theta_sigmas_ofs << yields.GetBinArg(i) << ",";
    if (i % 4 == 0) {
      script.NewPage(4,2).cd(s13::ana::PadSeq::column);
    }

    TString title = s13::ana::tstrfmt("Signal/BG (#theta: %.1f)", yields.GetBinArg(i));
    TH1* hist_signal = signals.at(i).AsHist(title);
    TH1* hist_bg = bgs.at(i).AsHist("BG");
    double theta = yields.GetBinArg(i);
    double phiW = cuts.GetPhiCutWidth(theta);
    hist_bg->Scale(bg_norm_fact_interp.eval(theta, phiW) - 1);

    script.cd();
    hist_signal->SetMarkerColor(kBlue);
    hist_signal->SetLineColor(kBlue);
    hist_signal->Draw("HIST");
    hist_bg->SetLineColor(kRed);
    hist_bg->SetMarkerColor(kRed);
    hist_bg->Draw("SAME HIST");

    script.cd();
    TH1* signal_wo_bg = subtract_bin_bg(hist_signal, hist_bg, yields.GetBinArg(i));
    signal_wo_bg->SetTitle("Signal - BG");
    FitParams fps = fit_theta_corr(signal_wo_bg, -cuts.theta_width / 2, cuts.theta_width / 2);
    theta_sigmas_ofs << fps;
    signal_wo_bg->Draw();

    s13::ana::draw_marker_line_ysym(cuts.theta_width / 2,
                                    signal_wo_bg->GetMinimum(), signal_wo_bg->GetMaximum());
    signals_wo_bg.push_back(signal_wo_bg);
  }

  theta_sigmas_ofs.close();
  return signals_wo_bg;
}

void
draw_bin_signal_bg_levels_carbon(Script& script,
                                 const std::vector<Yield>& signals,
                                 const std::vector<Yield>& bgs,
                                 const Yield& yields,
                                 const Cuts& cuts) {
  s13::misc::MessageLogger logger("draw_bin_signal_bg_levels_carbon()");
  std::ofstream theta_sigmas_ofs;
  theta_sigmas_ofs.open(get_output_filename(cuts));
  theta_sigmas_ofs << "# Gaussian fiting of theta "
                   << "dists. (theta, gaus_const, gaus_mean, "
                   << "gaus_sigma, errs...)" << std::endl;

  for (int i = 0; i < yields.GetBinsNumber(); ++i) {
    logger.debug("Drawing for theta %.1f...", yields.GetBinArg(i));
    theta_sigmas_ofs << yields.GetBinArg(i) << ",";
    if (i % 4 == 0) {
      script.NewPage(4,2).cd(s13::ana::PadSeq::column);
    }

    TString title = s13::ana::tstrfmt("Signal/BG (#theta: %.1f)", yields.GetBinArg(i));
    TH1* hist_signal = signals.at(i).AsHist(title, "#theta_{exp} - #theta_{theor}");
    TH1* hist_bg = scale_carbon_bg(bgs.at(i)).AsHist("BG");
    // TH1* hist_bg = bgs.at(i).AsHist("BG");

    script.cd();
    hist_signal->Draw("HIST");
    hist_bg->SetLineColor(kRed);
    hist_bg->Draw("SAME HIST");

    script.cd();
    TH1* signal_wo_bg = static_cast<TH1*>(hist_signal->Clone());
    signal_wo_bg->SetName(TString::Format("%s_nobg", signal_wo_bg->GetName()));
    signal_wo_bg->Add(hist_bg, -1);
    // signal_wo_bg->SetMinimum(0);
    signal_wo_bg->SetTitle("Signal - BG");

    signal_wo_bg->SetLineColor(kGreen);
    FitParams fps = fit_theta_corr(signal_wo_bg, -cuts.theta_width / 2, cuts.theta_width / 2);
    theta_sigmas_ofs << fps;

    signal_wo_bg->Draw();

  }
}

void compute_theta_sigmas_he4(Script& script, const Cuts& cuts) {
  auto algs = load_yields(cli_opts.dataset_type());
  auto signal_yield = algs.first->total_yields();
  auto bg_yield = algs.first->bg_yields();
  draw_bin_signal_bg_levels(script,
                            DsType::he4,
                            cuts,
                            algs.second->signal_yields(),
                            algs.second->bgs_yields(),
                            algs.second->total_yields());
}

void compute_theta_sigmas_he6(Script& script, const Cuts& cuts) {
  auto algs = load_yields(cli_opts.dataset_type());
  auto signal_yield = algs.first->total_yields();
  auto bg_yield = algs.first->bg_yields();
  auto algs_carbon = load_yields(DsType::carbon);
  auto carbon_yield = algs_carbon.first->total_yields();
  auto carbon_yield_scaled = scale_carbon_bg(carbon_yield);

  auto he_yield = signal_yield;
  auto he_yield_wo_bg = signal_yield - carbon_yield_scaled;
  std::cout << "Yields of pol. target: " << he_yield;
  std::cout << "Yields of pol. target wo BG: " << he_yield_wo_bg;
  draw_bin_signal_bg_levels_carbon(script,
                                   algs.second->signal_yields(),
                                   algs_carbon.second->signal_yields(),
                                   algs_carbon.second->total_yields(),
                                   cuts);
  draw_total_kin_corr(script,
                      algs.second->signal_yields(),
                      algs_carbon.second->signal_yields(),
                      cuts.GetTotalCutTitle(),
                      true);
}

void theta_cut_sigmas() {
  Script script("theta_cut_sigmas");
  gStyle->SetOptFit();
  Cuts cuts =
    s13::io::parse_cuts_config(std::fstream(cli_opts.cuts_config()));
  cuts.PrepareVarThetaCut(cli_opts.dataset_type());

  if (cli_opts.dataset_type() == DsType::he4) {
    compute_theta_sigmas_he4(script, cuts);
  } else if (cli_opts.dataset_type() == DsType::he6) {
    compute_theta_sigmas_he6(script, cuts);
  } else {
    throw std::invalid_argument("Invalid argument for CS drawing");
  }
}


int main(int argc, char** argv) {
  s13::misc::MessageLogger logger("main()");
  // s13::ana::SetRootStyle(s13::ana::RootStyle::publication);

  try {
    cli_opts.parse_options(argc, argv);
  } catch (std::exception& ex) {
    logger.error("Invalid argument provided: %s", ex.what());
    return -1;
  }

  // script.SetName(cli_opts.output_file_basename("cs_"));
  theta_cut_sigmas();
}
