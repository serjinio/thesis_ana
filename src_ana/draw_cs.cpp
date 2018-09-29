/**
 *   \file draw_cs.cpp
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
using Espri = s13::ana::Espri;
using DsType = s13::ana::DatasetType;
using PolDir = s13::ana::PolarizationDirection;
using Dsdcs = s13::ana::Dsdcs;

static s13::misc::CliOptions cli_opts;


std::pair<std::shared_ptr<s13::ana::ElasticScatteringYieldsAlg>,
          std::shared_ptr<s13::ana::SignalBgMonitorAlg> >
load_yields(DsType ds_type, const Cuts& cuts) {
  s13::misc::MessageLogger logger("load_yields()");
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
  filename += "_" + ds_type_str + "_" + pol_dir_str + "_" + espri_sel + ".csv";
  logger.debug("Trying to load ScatteringYield object from file: %s...",
               filename.c_str());
  std::ifstream ifs(filename);
  if (!ifs) {
    throw std::logic_error(s13::ana::tstrfmt("Cannot find specified file: %s",
                                             filename.c_str()).Data());
  }

  auto p_yields = s13::ana::make_elastic_scattering_yield_alg_from_file(ifs);
  logger.debug("Loaded scattering yield:");
  std::cout << p_yields->total_yields();
  logger.debug("Loaded scattering yield of BG:");
  std::cout << p_yields->bg_yields();
  logger.debug("Loading bin monitor results...");
  auto p_sig_bg = s13::ana::make_signal_bg_monitor_alg_from_file(ifs);

  return std::pair<decltype(p_yields), decltype(p_sig_bg)>(p_yields, p_sig_bg);
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

double draw_bin_sn(double angle, TH1* bin_hist, const Cuts& cuts) {
  s13::misc::MessageLogger logger("draw_bin_sn()");
  const double bg_min_angle = -4;
  const double bg_max_angle = -cuts.theta_width / 2;

  double bg_int =
    bin_hist->Integral(bin_hist->GetXaxis()->FindBin(bg_min_angle),
                       bin_hist->GetXaxis()->FindBin(bg_max_angle));
  double signal_int =
    bin_hist->Integral(bin_hist->GetXaxis()->FindBin(-cuts.theta_width/2),
                       bin_hist->GetXaxis()->FindBin(cuts.theta_width/2));
  double sn = signal_int / bg_int;
  s13::ana::draw_text(s13::ana::tstrfmt("S/N: %.2f", sn));
  logger.debug("S/N for bin at %.2f deg.: %.2f (S: %.2f; N: %.2f)",
               angle, sn, signal_int, bg_int);

  double hist_ymin = bin_hist->GetMinimum();
  double hist_ymax = bin_hist->GetMaximum() * 0.3;
  s13::ana::draw_marker_line_ysym(cuts.theta_width / 2, hist_ymin, hist_ymax);
  s13::ana::draw_marker_line(bg_min_angle, hist_ymin, bg_min_angle, hist_ymax);
  s13::ana::draw_marker_line(bg_max_angle, hist_ymin, bg_max_angle, hist_ymax);
  return sn;
}

std::vector<double>
draw_bin_signal_bg_levels(Script& script,
                          DsType ds_type,
                          const Cuts& cuts,
                          const std::vector<Yield>& signals,
                          const std::vector<Yield>& bgs,
                          const Yield& yields) {
  s13::misc::MessageLogger logger("draw_bin_signal_bg_levels()");
  std::vector<TH1*> signals_wo_bg;
  std::vector<double> sn_ratios;
  auto bg_norm_fact_interp =
    s13::ana::make_bg_norm_factor_interp(ds_type,
                                         cuts.dsdc_selector,
                                         cuts.use_var_phi_cut);

  logger.debug("Drawing signal & BGs hists...");
  logger.debug("Yields object bins num: %d", yields.GetBinsNumber());
  for (int i = 0; i < yields.GetBinsNumber(); ++i) {
    logger.debug("Drawing for theta %.1f...", yields.GetBinArg(i));
    if (i % 4 == 0) {
      script.NewPage(4,2).cd(s13::ana::PadSeq::column);
    }

    TString title = s13::ana::tstrfmt("Signal/BG (#theta: %.1f)", yields.GetBinArg(i));
    TH1* hist_signal = signals.at(i).AsHist(title);
    TH1* hist_bg = bgs.at(i).AsHist("BG");
    hist_bg->Scale(bg_norm_fact_interp->Eval(yields.GetBinArg(i)) - 1);

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
    signal_wo_bg->Draw("HIST");
    double sn = draw_bin_sn(yields.GetBinArg(i), signal_wo_bg, cuts);
    sn_ratios.push_back(sn);
    // s13::ana::draw_marker_line_ysym(cuts.theta_width / 2,
    //                                 signal_wo_bg->GetMinimum(), signal_wo_bg->GetMaximum());
    signals_wo_bg.push_back(signal_wo_bg);
  }

  return sn_ratios;
}

Yield
scale_carbon_bg(Yield carbon_yield) {
  auto scaled_yields = carbon_yield.multiply(s13::gk_carbon_bg_norm_factor);
  return scaled_yields;
}

std::vector<double>
draw_bin_signal_bg_levels_carbon(Script& script,
                                 const Cuts& cuts,
                                 const std::vector<Yield>& signals,
                                 const std::vector<Yield>& bgs,
                                 const Yield& yields) {
  s13::misc::MessageLogger logger("draw_bin_signal_bg_levels_carbon()");
  std::vector<TH1*> signals_wo_bg;
  std::vector<double> sn_ratios;
  logger.debug("Drawing he signal / carbon bg histograms...");
  for (int i = 0; i < yields.GetBinsNumber(); ++i) {
    logger.debug("Drawing for theta %.1f...", yields.GetBinArg(i));
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
    signal_wo_bg->Draw("HIST");
    signals_wo_bg.push_back(signal_wo_bg);
    double sn = draw_bin_sn(yields.GetBinArg(i), signal_wo_bg, cuts);
    sn_ratios.push_back(sn);
  }

  return sn_ratios;
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

/*
  Removes background contribution using out phi method.

  Returns: ScatteringYield without BG.
*/
Yield
extract_out_phi_bg(DsType ds_type, const Cuts& cuts, Yield total_yields, Yield bg_yields) {
  auto bg_extractor = s13::ana::make_out_phi_bgextractor(ds_type, cuts);
  Yield elastic_yields =
    bg_extractor->subtract_bg(total_yields, bg_yields);
  return elastic_yields;
}

double norm_dist_p_norm_fact(double dist_sigma, double cut_width) {
  double hw = cut_width / 2;
  double p_in = TMath::Erf(hw * (TMath::Sqrt(1/2.) / dist_sigma));
  return 1 / p_in;
}

double
vertz_cut_width_norm_fact(double vertz_cut_width) {
  return norm_dist_p_norm_fact(s13::gk_vertz_distribution_sigma, vertz_cut_width);
}

double
phi_cut_width_norm_fact(double phi_cut_width, double sigma) {
  return norm_dist_p_norm_fact(sigma, phi_cut_width);
}

/*
  Converts elastic scattering signal_yield into abosolute units:
  counts -> mbarns.
*/
Yield
normalize_ela_yields(DsType ds_type, PolDir pol_dir, Yield elastic_yields, const Cuts& cuts) {
  s13::misc::MessageLogger logger("normalize_ela_yields()");
  double num_incident = s13::ana::get_num_beam(ds_type, pol_dir, cuts.vertexXY_radius);
  logger.info("Total number of incident particles: %.9e", num_incident);

  s13::ana::RdcSolidAngle2
    solid_angle(s13::ana::CoordinateFrameType::lab,
                cuts.espri_selector,
                cuts.espri_y_width * 2,
                s13::gk_cs_bin_number, s13::gk_cs_range_start, s13::gk_cs_range_end);

  double target_density = 0;
  if (ds_type == DsType::he6 || ds_type == DsType::he4) {
    target_density = s13::gk_pol_tgt_areal_number_density;
  } else {
    target_density = s13::gk_carbon_tgt_areal_number_density;
  }
  logger.info("Target areal density: %.9e", target_density);

  auto yields = s13::ana::
    cs_from_yields(elastic_yields, solid_angle,
                   num_incident, target_density);

  // P-norm for different cuts
  if (cuts.apply_phi_cut && !cuts.use_var_phi_cut) {
    double sigma = cuts.dsdc_selector == Dsdcs::s1dc ?
      s13::gk_phi_distribution_sigma_s1dc :
      s13::gk_phi_distribution_sigma_fdc0;
    yields = yields.multiply(phi_cut_width_norm_fact(cuts.phi_width, sigma));
  }
  if (cuts.apply_espri_vertz_cut) {
    yields = yields.multiply(vertz_cut_width_norm_fact(cuts.espri_vertz_width));
  }

  return yields;
}

Double_t compute_edeg_thickness(Double_t proton_theta) {
  static const Double_t edeg_dist_from_target = 1123;
  static const Double_t edeg_base_angle = 53 * s13::gk_d2r;
  // static const Double_t edeg_finish_angle = 65 * D2R;
  static const Double_t edeg_thickness = 29;
  static const Double_t edeg_length = 255;
  static const Double_t edeg_min_thickness = 2;
  static const Double_t edeg_angle_tan = edeg_thickness / edeg_length;

  if (proton_theta < edeg_base_angle) {
    return 0;
  }

  Double_t dist_from_base = 2 * edeg_dist_from_target *
    TMath::Tan((proton_theta - edeg_base_angle) / 2);
  Double_t length2 = edeg_length - dist_from_base;
  Double_t thickness = length2 * edeg_angle_tan;
  if (thickness < edeg_min_thickness) {
    thickness = 0;
  }
  // std::cout << "Proton angle: " << proton_theta / D2R <<
  //   "; edeg thickness: " << thickness << std::endl;

  return thickness;
}

Yield correct_edeg_losses(Yield cs) {
  s13::misc::MessageLogger logger("compute_edeg_losses()");
  Yield corr_facts(cs.GetBinsNumber(), cs.RangeStart(), cs.RangeEnd());
  for(int i = 0; i < cs.GetBinsNumber(); ++i) {
    Double_t p_theta = cs.GetBinArg(i) * s13::gk_d2r;
    Double_t p_theta_lab = (TMath::Pi() - p_theta) / 2;
    Double_t edeg_z = compute_edeg_thickness(p_theta_lab);
    Double_t att_factor = 1 + (0.896 - 1) / 20 * edeg_z;
    logger.debug("For proton theta CM: %.2f (theta lab: %.2f) [deg.] "
                 "edeg_z: %.2f [mm] att. factor: %.2f",
                 p_theta/s13::gk_d2r, p_theta_lab/s13::gk_d2r,
                 edeg_z, att_factor);
    corr_facts.SetBinValue(cs.GetBinArg(i), 1 / att_factor);
    corr_facts.SetBinError(cs.GetBinArg(i), 0);
  }
  logger.info("Computed correction factors for E-deg reaction loss:");
  std::cout << corr_facts;
  return cs * corr_facts;
}

Yield
correct_det_efficiencies(Yield cs, Espri espri_selector,
                         const Cuts& cuts) {
  assert(espri_selector != Espri::both);
  s13::misc::MessageLogger logger("correct_det_efficiencies()");
  bool do_correction = true;

  if (!do_correction) {
    logger.info("Detectors correction for efficiencies is turned off, "
                "skipping this step!");
    return cs;
  }

  static const double rdc_left_opt_eff = 0.736;
  static const double rdc_right_opt_eff = 0.911;
  static const double rdc_left_eff_all_planes = 0.375;
  static const double rdc_right_eff_all_planes = 0.623;
  static const double fdc0_eff = 0.986;
  static const double s1dc_eff = 0.889;
  auto corr_cs = cs;
  if (espri_selector == Espri::left) {
    double eff = rdc_left_opt_eff;
    if ((cuts.apply_rdc_xnhitplanes_cut && cuts.rdc_xnhitplanes == 4)
        && (cuts.apply_rdc_ynhitplanes_cut && cuts.rdc_ynhitplanes == 3)) {
      eff = rdc_left_eff_all_planes;
    }
    logger.info("RDC left corrected using efficiency value: %.2f%%", eff * 100);
    corr_cs = cs.multiply(1./eff);
  } else if (espri_selector == Espri::right) {
    double eff = rdc_right_opt_eff;
    if ((cuts.apply_rdc_xnhitplanes_cut && cuts.rdc_xnhitplanes == 4)
        && (cuts.apply_rdc_ynhitplanes_cut && cuts.rdc_ynhitplanes == 3)) {
      eff = rdc_right_eff_all_planes;
    }
    logger.info("RDC right corrected using efficiency value: %.2f%%", eff * 100);
    corr_cs = cs.multiply(1./eff);
  }

  if (cuts.dsdc_selector == Dsdcs::fdc0) {
    logger.info("Correcting scattering signal yields for FDC0 efficiency: %.2f%%",
                fdc0_eff * 100);
    corr_cs = corr_cs.multiply(1./fdc0_eff);
  } else if (cuts.dsdc_selector == Dsdcs::s1dc) {
    logger.info("Correcting scattering signal yields for S1DC efficiency: %.2f%%",
                s1dc_eff * 100);
    corr_cs = corr_cs.multiply(1./s1dc_eff);
  } else {
    throw std::invalid_argument("Invalid DSDC selector!");
  }

  // NOTE: custom factor! USE ONLY FOR TESTING!!!
  // corr_cs = corr_cs.multiply(0.7); // rough solid angle correction for NaIs

  return corr_cs;
}

Yield
yields_to_cm_cs(Yield ela_yields,
                const Cuts& cuts,
                DsType ds_type) {
  auto lab_cs = normalize_ela_yields(ds_type, cli_opts.pol_direction(),
                                     ela_yields, cuts);
  auto cm_cs = s13::ana::lab_cs_to_cm(lab_cs, ds_type);
  cm_cs = correct_det_efficiencies(cm_cs, cuts.espri_selector, cuts);
  cm_cs = correct_edeg_losses(cm_cs);
  std::cout << "CM CS: " << cm_cs << std::endl;
  return cm_cs;
}

Yield he6_cs_yield(const Cuts& cuts, bool subtract_carbon_bg = true) {
  auto algs = load_yields(DsType::he6, cuts);
  auto signal_yield = algs.first->total_yields();
  if (subtract_carbon_bg) {
    auto algs_carbon = load_yields(DsType::carbon, cuts);
    auto carbon_yield = algs_carbon.first->total_yields();
    auto carbon_yield_scaled = scale_carbon_bg(carbon_yield);
    signal_yield = signal_yield - carbon_yield_scaled;
  }

  // TODO: fix this method now returns just yields
  // auto he_cm_cs = yields_to_cm_cs(signal_yield,
  //                                 cli_opts.espri_selector(),
  //                                 DsType::he6);
  return signal_yield;
}

TGraph* draw_sn_ratios_graph(Script& script,
                             std::vector<double> sn_ratios_data,
                             DsType ds_type) {
  Yield sn_ratios(s13::gk_cs_bin_number, s13::gk_cs_range_start, s13::gk_cs_range_end);
  Yield inv_sn_ratios(s13::gk_cs_bin_number, s13::gk_cs_range_start, s13::gk_cs_range_end);
  for (int i = 0; i < inv_sn_ratios.GetBinsNumber(); i++) {
    auto arg = inv_sn_ratios.GetBinArg(i);
    sn_ratios.SetBinValue(arg, sn_ratios_data[i]);
    sn_ratios.SetBinError(arg, 0);
    inv_sn_ratios.SetBinValue(arg, 1 / sn_ratios_data[i]);
    inv_sn_ratios.SetBinError(arg, 0);
  }
  sn_ratios = s13::ana::lab_yield_to_cm(sn_ratios, ds_type);
  inv_sn_ratios = s13::ana::lab_yield_to_cm(inv_sn_ratios, ds_type);

  auto title = s13::ana::tstrfmt("Bin data: (S/N)");
  auto sn_graph = s13::ana::
    BuildTGraphErrors(sn_ratios, title, "CM angle [deg]", "(S/N)");
  script.NewPage(1,1).cd();
  sn_graph->Draw("ACP");
  // sn_graph->SetMinimum(0);

  title = s13::ana::tstrfmt("Bin data: 1 / (S/N)");
  auto inv_sn_graph = s13::ana::
    BuildTGraphErrors(inv_sn_ratios, title, "CM angle [deg]", "1 / (S/N)");
  script.NewPage(1,1).cd();
  inv_sn_graph->Draw("ACP");
  inv_sn_graph->SetMinimum(-0.5);
  inv_sn_graph->SetMaximum(0.5);

  return sn_graph;
}

TGraph* make_cs_graph(TString title, Yield cs) {
  auto cs_graph = s13::ana::
    BuildTGraphErrors(cs, title, "CM angle [deg]", "DCS [mb/sr]");
  return cs_graph;
}

TMultiGraph* make_cs_mgraph(TString title, std::vector<TGraph*> graphs) {
  TMultiGraph* mg = new TMultiGraph();
  mg->SetTitle(TString::Format("%s;%s;%s", title.Data(),
                               "CM angle [deg.]", "DCS [mb/sr]"));
  for (auto* el : graphs) {
    mg->Add(el);
  }
  return mg;
}

void draw_cs_mgraph(Script& script, TMultiGraph* mg) {
  script.NewPage(1,1).cd(1);
  mg->Draw("ACP");
  gPad->SetLogy();
}

TMultiGraph* draw_he4_cs(Script& script, Yield cs, TString title = "") {
  auto tit = s13::ana::tstrfmt("CS p-4He@200MeV (%s)", title.Data());
  auto* cs_graph = make_cs_graph(tit, cs);
  auto* mg = make_cs_mgraph("He4 CS", {s13::ana::
        BuildMossCmsCsGraph("data/he4_200mev_moss_cs.csv")});
  draw_cs_mgraph(script, mg);
  cs_graph->Draw("P SAME");
  // gPad->BuildLegend(0.15, 0.15, 0.45, 0.35, "", "PL");
  gPad->BuildLegend(0.15, 0.15, 0.45, 0.35, "");

  return mg;
}

TMultiGraph* draw_he6_cs(Script& script, Yield cs, TString title = "") {
  s13::misc::MessageLogger logger("draw_he6_cs()");
  auto tit = s13::ana::tstrfmt("CS p-6He@200MeV (%s)", title.Data());
  auto* cs_graph = make_cs_graph(tit, cs);
  auto* mg = make_cs_mgraph("He6 CS", {});
  s13::ana::AddHe6TheorCses(mg);
  draw_cs_mgraph(script, mg);
  cs_graph->Draw("P SAME");
  mg->SetMaximum(0.7);
  // gPad->BuildLegend(0.6, 0.6, 0.90, 0.9, "", "PL");
  gPad->BuildLegend(0.6, 0.6, 0.90, 0.9, "");
  return mg;
}

TMultiGraph* draw_qfs_cs(Script& script, Yield cs, const Cuts& cuts) {
  s13::misc::MessageLogger logger("draw_he6_cs()");
  auto* cs_graph = make_cs_graph( "QFS", cs);
  cs_graph->SetLineColor(kRed);
  auto* he6_cs_graph = make_cs_graph("p-6He (pol. target)", he6_cs_yield(cuts, false));
  auto* mg = make_cs_mgraph("6He QFS CS", {cs_graph, he6_cs_graph});
  draw_cs_mgraph(script, mg);
  return mg;
}

void write_cs_to_file(std::string filename, Yield cs) {
  std::ofstream ofs(filename);
  cs.Serialize(ofs);
}

Yield do_compute_he4_cs(Script& script, const Cuts& cuts) {
  assert(cuts.espri_selector != Espri::both);
  auto algs = load_yields(cli_opts.dataset_type(), cuts);
  auto signal_yield = algs.first->total_yields();
  auto bg_yield = algs.first->bg_yields();
  auto sn_ratios =
    draw_bin_signal_bg_levels(script,
                              DsType::he4,
                              cuts,
                              algs.second->signal_yields(),
                              algs.second->bgs_yields(),
                              algs.second->total_yields());
  // draw_total_kin_corr(script,
  //                     algs.second->signal_yields(),
  //                     algs.second->bgs_yields(),
  //                     cuts.GetTotalCutTitle(),
  //                     false);
  auto ela_yield = extract_out_phi_bg(DsType::he4, cuts,
                                      signal_yield, bg_yield);
  auto he_cm_cs = yields_to_cm_cs(ela_yield,
                                  cuts, DsType::he4);
  draw_sn_ratios_graph(script, sn_ratios, DsType::he4);
  draw_he4_cs(script, he_cm_cs,
              cuts.espri_selector == Espri::left ? "ESL" : "ESR");
  return he_cm_cs;
}

void compute_he4_cs(Script& script, Cuts& cuts) {
  cuts.espri_selector = Espri::left;
  Yield left_cs = do_compute_he4_cs(script, cuts);
  write_cs_to_file("cs_data/cs_left.csv", left_cs);
  cuts.espri_selector = Espri::right;
  Yield right_cs = do_compute_he4_cs(script, cuts);
  write_cs_to_file("cs_data/cs_right.csv", right_cs);

  Yield total_cs = (left_cs + right_cs).multiply(1/2.);
  write_cs_to_file("cs_data/cs.csv", total_cs);
  draw_he4_cs(script, total_cs, "both ESPRIs");
}

Yield do_compute_he6_cs(Script& script, const Cuts& cuts) {
  assert(cuts.espri_selector != Espri::both);
  auto algs = load_yields(cli_opts.dataset_type(), cuts);
  auto signal_yield = algs.first->total_yields();
  auto bg_yield = algs.first->bg_yields();
  auto algs_carbon = load_yields(DsType::carbon, cuts);
  auto carbon_yield = algs_carbon.first->total_yields();
  auto carbon_yield_scaled = scale_carbon_bg(carbon_yield);

  auto he_yield = signal_yield;
  auto he_yield_wo_bg = signal_yield - carbon_yield_scaled;
  std::cout << "Yields of pol. target: " << he_yield;
  std::cout << "Yields of pol. target wo BG: " << he_yield_wo_bg;
  auto sn_ratios =
    draw_bin_signal_bg_levels_carbon(script,
                                     cuts,
                                     algs.second->signal_yields(),
                                     algs_carbon.second->signal_yields(),
                                     algs_carbon.second->total_yields());
  draw_total_kin_corr(script,
                      algs.second->signal_yields(),
                      algs_carbon.second->signal_yields(),
                      cuts.GetTotalCutTitle(),
                      true);
  draw_sn_ratios_graph(script, sn_ratios, DsType::he6);

  // auto he_yield_wo_bg = extract_out_phi_bg(DsType::he6, cuts.dsdc_selector,
  //                                          signal_yield, bg_yield);
  // draw_bin_signal_bg_levels(script,
  //                           DsType::he6,
  //                           algs.second->signal_yields(),
  //                           algs.second->bgs_yields(),
  //                           algs.second->total_yields());

  auto he_cm_cs = yields_to_cm_cs(he_yield_wo_bg,
                                  cuts, cli_opts.dataset_type());
  draw_he6_cs(script, he_cm_cs,
              cuts.espri_selector == Espri::left ? "ESL" : "ESR");
  return he_cm_cs;
}

void compute_he6_cs(Script& script, Cuts& cuts) {
  cuts.espri_selector = Espri::left;
  Yield left_cs = do_compute_he6_cs(script, cuts);
  write_cs_to_file("cs_data/cs_left.csv", left_cs);
  cuts.espri_selector = Espri::right;
  Yield right_cs = do_compute_he6_cs(script, cuts);
  write_cs_to_file("cs_data/cs_right.csv", right_cs);

  Yield total_cs = (left_cs + right_cs).multiply(1/2.);
  write_cs_to_file("cs_data/cs.csv", total_cs);
  draw_he6_cs(script, total_cs, "both ESPRIs");
}

void compute_qfs_cs(Script& script, const Cuts& cuts) {
  auto algs_carbon = load_yields(DsType::carbon, cuts);
  auto carbon_yield = algs_carbon.first->total_yields();
  //auto carbon_yield_scaled = scale_carbon_bg(carbon_yield);
  auto cm_cs = yields_to_cm_cs(carbon_yield,
                               cuts, cli_opts.dataset_type());
  draw_bin_signal_bg_levels_carbon(script,
                                   cuts,
                                   algs_carbon.second->signal_yields(),
                                   algs_carbon.second->bgs_yields(),
                                   algs_carbon.second->total_yields());
  draw_qfs_cs(script, cm_cs, cuts);
  write_cs_to_file("cs_data/qfs_cm_cs.csv", cm_cs);
}

void draw_cs() {
  Script script("draw_cs", "cs_data");
  Cuts cuts =
    s13::io::parse_cuts_config(std::fstream(cli_opts.cuts_config()));
  cuts.PrepareVarThetaCut(cli_opts.dataset_type());

  if (cli_opts.dataset_type() == DsType::he4) {
    compute_he4_cs(script, cuts);
  } else if (cli_opts.dataset_type() == DsType::he6) {
    compute_he6_cs(script, cuts);
  } else if (cli_opts.dataset_type() == DsType::carbon) {
    compute_qfs_cs(script, cuts);
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
  draw_cs();
}
