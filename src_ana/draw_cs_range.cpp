/**
 *   \file draw_cs_range.cpp
 *   \brief Draws CS given yields as input.
 *
 */

#include <numeric>
#include <type_traits>

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


std::string
get_out_fname(std::string data_name, DsType ds_type,
          PolDir pol_dir, const Cuts& cuts) {
  auto ds_type_str = ds_type == DsType::carbon ? "carbon" :
    cli_opts.dataset_type_str();
  auto pol_dir_str = ds_type != DsType::carbon ?
    cli_opts.pol_direction_str() : "pol_both";
  TString fname = s13::ana::tstrfmt("cs_data/%s_%s_%s_%s.csv",
                                    ds_type_str.data(),
                                    pol_dir_str.data(),
                                    cuts.GetTotalCutName().Data(),
                                    data_name.data());
  return fname.Data();
}

std::shared_ptr<s13::ana::SignalBgMonitorAlg>
load_snmon_data(const Cuts& cuts, DsType ds_type) {
  s13::misc::MessageLogger logger("load_snmon_data()");
  auto ds_type_str = ds_type == DsType::carbon ? "carbon" :
    cli_opts.dataset_type_str();
  auto pol_dir_str = ds_type != DsType::carbon ?
    cli_opts.pol_direction_str() : "pol_both";
  TString snmon_fname = s13::ana::tstrfmt("cs_data/%s_%s_%s_snmon.csv",
                                           ds_type_str.data(),
                                           pol_dir_str.data(),
                                           cuts.GetTotalCutName().Data());
  logger.debug("Trying to load S/N monitor data from file: %s...",
               snmon_fname.Data());
  std::ifstream snmon_ifs(snmon_fname);
  if (!snmon_ifs) {
    throw std::logic_error(s13::ana::tstrfmt("Cannot find specified file: %s",
                                             snmon_fname.Data()).Data());
  }

  auto p_sig_bg = s13::ana::make_signal_bg_monitor_alg_from_file(snmon_ifs);
  return p_sig_bg;
}

std::shared_ptr<s13::ana::ElasticScatteringYieldsAlg>
load_yields_data(const Cuts& cuts, DsType ds_type) {
  s13::misc::MessageLogger logger("load_yields_data()");
  auto ds_type_str = ds_type == DsType::carbon ? "carbon" :
    cli_opts.dataset_type_str();
  auto pol_dir_str = ds_type != DsType::carbon ?
    cli_opts.pol_direction_str() : "pol_both";
  TString yields_fname = s13::ana::tstrfmt("cs_data/%s_%s_%s_yield.csv",
                                           ds_type_str.data(),
                                           pol_dir_str.data(),
                                           cuts.GetTotalCutName().Data());
  logger.debug("Trying to load ScatteringYield object from file: %s...",
               yields_fname.Data());
  std::ifstream yields_ifs(yields_fname);
  if (!yields_ifs) {
    throw std::logic_error(s13::ana::tstrfmt("Cannot find specified file: %s",
                                             yields_fname.Data()).Data());
  }

  auto p_yields = s13::ana::make_elastic_scattering_yield_alg_from_file(yields_ifs);
  logger.debug("Loaded scattering yield:");
  std::cout << p_yields->total_yields();
  logger.debug("Loaded scattering yield of BG:");
  std::cout << p_yields->bg_yields();
  logger.debug("Loading bin monitor results...");

  return p_yields;
}

TH1* bin_hist_subtract_bg(TH1* signal, TH1* bg, double lab_theta) {
  TH1* signal_wo_bg = static_cast<TH1*>(signal->Clone());
  signal_wo_bg->Add(bg, -1);
  // signal_wo_bg->SetMinimum(0);
  signal_wo_bg->SetName(s13::ana::tstrfmt("%s-nobg", signal->GetName()).Data());
  signal_wo_bg->SetLineColor(kGreen);
  signal_wo_bg->SetMarkerColor(kGreen);
  return signal_wo_bg;
}

double draw_bin_sn(double angle, TH1* bin_hist, const Cuts& cuts) {
  return s13::ana::draw_kin_corr_hist_sn(angle, bin_hist,
                                         cuts.dsdc_selector,
                                         cuts.theta_width,
                                         s13::ana::UseAbsVal(true));
}

Yield
bins_sn_vector_as_yield(std::vector<double> lab_bins_sn_data) {
  Yield sn_ratios(s13::gk_cs_bin_number, s13::gk_cs_range_start, s13::gk_cs_range_end);
  for (int i = 0; i < sn_ratios.GetBinsNumber(); i++) {
    auto arg = sn_ratios.GetBinArg(i);
    sn_ratios.SetBinValue(arg, lab_bins_sn_data[i]);
    sn_ratios.SetBinError(arg, 0);
  }
  return sn_ratios;
}

Yield
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
    s13::ana::make_bg_norm_factor_interp2(ds_type, cuts.dsdc_selector);

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
    double r_fact =
      bg_norm_fact_interp.eval(yields.GetBinArg(i),
                               cuts.GetPhiCutWidth(yields.GetBinArg(i)));
    hist_bg->Scale(r_fact - 1);

    script.cd();
    hist_signal->SetMarkerColor(kBlue);
    hist_signal->SetLineColor(kBlue);
    hist_signal->Draw("HIST");
    hist_bg->SetLineColor(kRed);
    hist_bg->SetMarkerColor(kRed);
    hist_bg->Draw("SAME HIST");

    script.cd();
    TH1* signal_wo_bg = bin_hist_subtract_bg(hist_signal, hist_bg, yields.GetBinArg(i));
    signal_wo_bg->SetTitle("Signal - BG");
    signal_wo_bg->SetMinimum(0);
    signal_wo_bg->Draw();
    signal_wo_bg->Fit("gaus", "", "", -cuts.theta_width/2, cuts.theta_width/2);
    // double sn = draw_bin_sn(yields.GetBinArg(i), signal_wo_bg, cuts);
    double sn = 1;
    sn_ratios.push_back(sn);
    // s13::ana::draw_marker_line_ysym(cuts.theta_width / 2,
    //                                 signal_wo_bg->GetMinimum(), signal_wo_bg->GetMaximum());
    signals_wo_bg.push_back(signal_wo_bg);
  }

  return bins_sn_vector_as_yield(sn_ratios);
}

Yield
scale_carbon_bg(DsType ds_type, Yield carbon_yield) {
  assert(ds_type != DsType::carbon && "Only he4 & he6 datasets allowed here!");
  s13::misc::MessageLogger logger("scale_carbon_bg()");
  double fact = ds_type == DsType::he6 ?
    s13::gk_carbon_bg_norm_factor : s13::gk_carbon_to_he4_bg_norm_factor;
  logger.info("Will use carbon BG scaling factor for %s dataset",
              ds_type == DsType::he6 ? "He6" : "He4");
  logger.info("Scaling factor value: %.2f.", fact);
  if (cli_opts.pol_direction() != PolDir::both) {
    logger.info("Carbon data scaling factor is reduced by 2 due to "
                "up/down polarization direction selected.");
    fact /= 2;
  }
  auto scaled_yields = carbon_yield.multiply(fact);
  return scaled_yields;
}

Yield
draw_bin_signal_bg_levels_carbon(Script& script,
                                 DsType ds_type,
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
    TH1* hist_bg = scale_carbon_bg(ds_type, bgs.at(i)).AsHist("BG");
    // TH1* hist_bg = bgs.at(i).AsHist("BG");

    script.cd();
    hist_signal->Draw("HIST");
    hist_bg->SetLineColor(kRed);
    hist_bg->Draw("SAME HIST");

    script.cd();
    TH1* signal_wo_bg = static_cast<TH1*>(hist_signal->Clone());
    signal_wo_bg->SetName(TString::Format("%s_nobg", signal_wo_bg->GetName()));
    // only remove C BG if configured in cuts && at BW angles
    // (due to low stats at FW)
    // if (cuts.subtract_he6_carbon_bg && yields.GetBinArg(i) > 60.3) {
    //   signal_wo_bg->Add(hist_bg, -1);
    // }
    if (cuts.subtract_he6_carbon_bg) {
      signal_wo_bg->Add(hist_bg, -1);
    }
    signal_wo_bg->SetTitle("Signal - BG");
    signal_wo_bg->SetLineColor(kGreen);
    // signal_wo_bg->Fit("gaus", "", "", -cuts.theta_width/2, cuts.theta_width/2);
    // signal_wo_bg->SetMinimum(-6);
    signal_wo_bg->Draw("HIST");
    signals_wo_bg.push_back(signal_wo_bg);
    double sn = draw_bin_sn(yields.GetBinArg(i), signal_wo_bg, cuts);
    // double sn = 1;
    sn_ratios.push_back(sn);
  }

  return bins_sn_vector_as_yield(sn_ratios);
}

void
draw_total_kin_corr_out_phi_bg(std::vector<Yield> yields_signal,
                        std::vector<Yield> yields_bg,
                        const Cuts& cuts) {
  auto bg_norm_fact_interp =
    s13::ana::make_bg_norm_factor_interp2(DsType::he4,
                                         cuts.dsdc_selector);
  double lab_angle = s13::gk_cs_range_start;
  double step = (s13::gk_cs_range_end - s13::gk_cs_range_start) / s13::gk_cs_bin_number;
  for (size_t i = 0; i < yields_signal.size(); i++) {
    double scaling_f =
      bg_norm_fact_interp.eval(lab_angle,
                               cuts.GetPhiCutWidth(lab_angle)) - 1;
    yields_bg[i] = yields_bg[i].multiply(scaling_f);
    lab_angle += step;
  }

  using SY = Yield;
  SY total_yield_signal = std::accumulate(yields_signal.begin(), yields_signal.end(),
                                          yields_signal[0].Empty());
  SY total_yield_bg = std::accumulate(yields_bg.begin(), yields_bg.end(),
                                      yields_bg[0].Empty());
  SY clean_yield(total_yield_signal - total_yield_bg);

  TH1* hist_tot =
    total_yield_signal.AsHist(-4,4, "TOT",
                              "#theta_{{}^{6}He} - #theta_{{}^{6}He-theor}");
  hist_tot->SetMaximum(hist_tot->GetMaximum() * 1.2);
  TH1* hist_bg =
    total_yield_bg.AsHist(-4,4, "BG",
                          "#theta_{{}^{6}He} - #theta_{{}^{6}He-theor}");
  TH1* hist_nobg =
    clean_yield.AsHist(-4,4, "TOT - BG",
                          "#theta_{{}^{6}He} - #theta_{{}^{6}He-theor}");

  s13::ana::draw_thstack({hist_tot, hist_bg, hist_nobg},
                         "p-4He kin. corr.", "NOSTACK HIST");
  gPad->BuildLegend(0.6, 0.65, 0.85, 0.85);
  s13::ana::draw_text(cuts.GetTotalCutName().Data(), 0.1, 0.8, 0.6, 0.9);
  gPad->Update();
}

void
draw_total_kin_corr_out_phi_bg(Script& script,
                        Cuts cuts) {
  script.NewPage(2,1);

  cuts.apply_theta_cut = false;
  cuts.espri_selector = Espri::left;
  auto snmon_yields = load_snmon_data(cuts, DsType::he4);
  script.cd();
  draw_total_kin_corr_out_phi_bg(snmon_yields->signal_yields(),
                          snmon_yields->bgs_yields(),
                          cuts);

  cuts.espri_selector = Espri::right;
  snmon_yields = load_snmon_data(cuts, DsType::he4);
  script.cd();
  draw_total_kin_corr_out_phi_bg(snmon_yields->signal_yields(),
                          snmon_yields->bgs_yields(),
                          cuts);
}

void
draw_total_kin_corr_he6(Script& script,
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
    total_yield_carbon = scale_carbon_bg(DsType::he6, total_yield_carbon);
  }

  TH1* hist_he =
    total_yield_he.AsHist(-4,4, "Total", "#theta_{{}^{6}He} - #theta_{{}^{6}He-theor}");
  hist_he->SetMaximum(hist_he->GetMaximum() * 1.2);
  TH1* hist_carbon =
    total_yield_carbon.AsHist(-4,4, "Carbon", "#theta_{exp} - #theta_{theor}");
  TH1* hist_he_wo_carbon = static_cast<TH1*>(hist_he->Clone());
  // TODO: (low-p) should only be subtracted if configured in cuts
  hist_he_wo_carbon->Add(hist_carbon, -1);
  hist_he_wo_carbon->SetTitle("Elastic");

  script.NewPage(1,1).cd();
  s13::ana::draw_thstack({hist_he_wo_carbon, hist_he, hist_carbon},
                         "He & Carbon yields", "NOSTACK HIST");
  gPad->BuildLegend(0.6, 0.65, 0.85, 0.85);
  s13::ana::draw_text(cut_conditions_str, 0.1, 0.8, 0.6, 0.9);
  gPad->Update();
}

Yield
extract_out_phi_bg(DsType ds_type, const Cuts& cuts, Yield total_yields, Yield bg_yields) {
  auto bg_in_yields = s13::ana::convert_to_in_phi_bg(ds_type, cuts, bg_yields);
  return total_yields - bg_in_yields;
}

double norm_dist_p_norm_fact(double dist_sigma, double cut_width) {
  double hw = cut_width / 2;
  double p_in = TMath::Erf(hw * (TMath::Sqrt(1/2.) / dist_sigma));
  return 1 / p_in;
}

double
phi_cut_width_norm_fact(double phi_cut_width, double sigma) {
  return norm_dist_p_norm_fact(sigma, phi_cut_width);
}

double
vertz_cut_width_norm_fact(double vertz_cut_width) {
  return norm_dist_p_norm_fact(s13::gk_vertz_distribution_sigma, vertz_cut_width);
}

double
e_cut_width_norm_fact(double e_cut_width_stdev) {
  return norm_dist_p_norm_fact(1.0, e_cut_width_stdev * 2);
}

double
theta_cut_width_norm_fact(double theta_cut_width, double sigma) {
  return norm_dist_p_norm_fact(sigma, theta_cut_width);
}

/*
  Converts elastic scattering signal_yield into abosolute units:
  counts -> mbarns.
*/
Yield
normalize_ela_yields(DsType ds_type, PolDir pol_dir, Yield elastic_yields, const Cuts& cuts) {
  s13::misc::MessageLogger logger("normalize_ela_yields()");
  // KLUDGE: Enum type carbon is for target, but beam actually is
  // he6. Just fix this issue by remapping in this function carbon -> He6
  DsType beam_ds_type = ds_type;
  if (ds_type == DsType::carbon) {
    beam_ds_type = DsType::he6;
  }
  double num_incident = s13::ana::get_num_beam(beam_ds_type, pol_dir,
                                               cuts.vertexXY_radius);
  logger.info("Total number of incident particles: %.9e", num_incident);

  s13::ana::RdcSolidAngle2
    solid_angle(s13::ana::CoordinateFrameType::lab,
                cuts.espri_selector,
                cuts.espri_y_width * 2,
                s13::gk_cs_bin_number, s13::gk_cs_range_start, s13::gk_cs_range_end);

  double target_density = 0;
  if (ds_type == DsType::he6 || ds_type == DsType::he4) {
    target_density = s13::gk_pol_tgt_areal_number_density;
  } else if (ds_type == DsType::carbon) {
    target_density = s13::gk_carbon_tgt_areal_number_density;
  } else {
    throw std::invalid_argument("Invalid dataset type!");
  }
  logger.info("Target areal density: %.9e", target_density);

  auto yields = s13::ana::
    cs_from_yields(elastic_yields, solid_angle,
                   num_incident, target_density);

  // P-norm for different cuts
  if (cuts.apply_phi_cut) {
    const double sigma = cuts.dsdc_selector == Dsdcs::s1dc ?
      s13::gk_phi_distribution_sigma_s1dc :
      s13::gk_phi_distribution_sigma_fdc0;
    const double norm_fact =
      phi_cut_width_norm_fact(cuts.phi_width * 2 * sigma, sigma);
    logger.info("dPhi cut P-norm factor: %.4f", norm_fact);
    yields = yields.multiply(norm_fact);
  }
  if (cuts.apply_espri_vertz_cut) {
    yields = yields.multiply(vertz_cut_width_norm_fact(cuts.espri_vertz_width));
  }
  if (cuts.apply_espri_e_cut) {
    logger.info("E cut P-norm factor: %.4f",
                e_cut_width_norm_fact(cuts.espri_e_width_stdev));
    yields = yields.multiply(e_cut_width_norm_fact(cuts.espri_e_width_stdev));
  }
  if (cuts.apply_theta_cut) {
    double sigma = cuts.dsdc_selector == Dsdcs::s1dc ?
      s13::gk_theta_distribution_sigma_s1dc :
      s13::gk_theta_distribution_sigma_fdc0;
    yields = yields.multiply(theta_cut_width_norm_fact(cuts.theta_width, sigma));
  }

  // normalization factor for solid angle
  // ratio of bin solid angle NaIs/RDC
  if (cuts.apply_espri_e_cut || cuts.apply_espri_min_e_de_cut) {
    if (cuts.espri_selector == Espri::left) {
      yields = yields.multiply(1/s13::gk_rdc2nai_esl_solid_angle_factor);
    } else if (cuts.espri_selector == Espri::right) {
      yields = yields.multiply(1/s13::gk_rdc2nai_esr_solid_angle_factor);
    } else {
      throw std::invalid_argument("ESPRI selector 'both' is invalid for "
                                  "yields normalization!");
    }
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
    Double_t p_theta_lab = cs.GetBinArg(i) * s13::gk_d2r;
    Double_t edeg_z = compute_edeg_thickness(p_theta_lab);
    Double_t att_factor = 1 + (0.896 - 1) / 20 * edeg_z;
    logger.debug("For proton theta lab: %.2f [deg.] "
                 "edeg_z: %.2f [mm] att. factor: %.2f",
                 p_theta_lab/s13::gk_d2r, edeg_z, att_factor);
    corr_facts.SetBinValue(cs.GetBinArg(i), 1 / att_factor);
    corr_facts.SetBinError(cs.GetBinArg(i), 0);
  }
  logger.info("Computed correction factors for E-deg reaction loss:");
  std::cout << corr_facts;
  return cs * corr_facts;
}

Yield
correct_det_efficiencies(Yield cs,
                         Espri espri_selector,
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
  lab_cs = correct_det_efficiencies(lab_cs, cuts.espri_selector, cuts);
  lab_cs = correct_edeg_losses(lab_cs);
  auto cm_cs = s13::ana::lab_cs_to_cm(lab_cs, ds_type);
  std::cout << "CM CS: " << cm_cs << std::endl;
  return cm_cs;
}

TGraph* draw_sn_ratios_graph(Script& script, Yield sn_ratios,
                             DsType ds_type) {
  sn_ratios = s13::ana::lab_yield_to_cm(sn_ratios, ds_type);
  auto title = s13::ana::tstrfmt("Bin data: (S/N)");
  auto sn_graph = s13::ana::
    BuildTGraphErrors(sn_ratios, title, "CM angle [deg]", "(S/N)");
  script.NewPage(1,1).cd();
  sn_graph->Draw("ACP");
  // sn_graph->SetMinimum(0);

  title = s13::ana::tstrfmt("Bin data: 1 / (S/N)");
  auto inv_sn_graph = s13::ana::
    BuildTGraphErrors(sn_ratios.invert(), title,
                      "CM angle [deg]", "1 / (S/N)");
  script.NewPage(1,1).cd();
  inv_sn_graph->Draw("ACP");
  inv_sn_graph->SetMinimum(-0.0);
  // inv_sn_graph->SetMaximum(0.5);

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

TMultiGraph* make_sn_mgraph(TString title, std::vector<TGraph*> graphs) {
  TMultiGraph* mg = new TMultiGraph();
  mg->SetTitle(TString::Format("%s;%s;%s", title.Data(),
                               "CM angle [deg.]", "1/(S/N)"));
  for (auto* el : graphs) {
    mg->Add(el, "PC");
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

void write_cs_to_file(TString yields_fname, Yield cs) {
  std::ofstream ofs(yields_fname);
  cs.Serialize(ofs);
}

Yield
compute_sn_out_phi_bg(const Cuts& cuts, DsType ds_type) {
  Script sn_script(s13::ana::tstrfmt("%s_%s_sn",
                                     cli_opts.dataset_type_str().data(),
                                     cuts.GetTotalCutName().Data()),
                   "cs_data/pdf");
  auto snmon_yields = load_snmon_data(cuts, ds_type);
  auto sn_ratios =
    draw_bin_signal_bg_levels(sn_script,
                              ds_type,
                              cuts,
                              snmon_yields->signal_yields(),
                              snmon_yields->bgs_yields(),
                              snmon_yields->total_yields());
  draw_sn_ratios_graph(sn_script, sn_ratios, ds_type);
  return sn_ratios;
}

Yield
compute_sn_carbon_bg(const Cuts& cuts, DsType ds_type) {
  Script sn_script(s13::ana::tstrfmt("%s_%s_sn",
                                     cli_opts.dataset_type_str().data(),
                                     cuts.GetTotalCutName().Data()),
                   "cs_data/pdf");
  auto snmon_yields = load_snmon_data(cuts, ds_type);
  using SnMonData = std::result_of<decltype(&load_snmon_data)(Cuts, DsType)>::type;
  SnMonData carbon_snmon_yields;
  if (ds_type == DsType::he6) {
    carbon_snmon_yields = load_snmon_data(cuts, DsType::carbon);
  } else if (ds_type == DsType::he4) {
    auto bg_cuts = cuts;
    bg_cuts.apply_hodf_he4_cut = true;
    carbon_snmon_yields = load_snmon_data(bg_cuts, DsType::carbon);
  }

  auto sn_ratios =
    draw_bin_signal_bg_levels_carbon(sn_script,
                                     ds_type,
                                     cuts,
                                     snmon_yields->signal_yields(),
                                     carbon_snmon_yields->signal_yields(),
                                     snmon_yields->total_yields());
  draw_sn_ratios_graph(sn_script, sn_ratios, ds_type);
  return sn_ratios;
}

using LrSnsRangePair = std::pair< std::vector< Yield >,
                                  std::vector< Yield > >;

template <typename SnProc>
LrSnsRangePair
draw_sn_range(Script& script, DsType ds_type,
              std::vector<Cuts>& cuts_range,
              SnProc sn_proc) {
  LrSnsRangePair res({}, {});
  for (auto& c : cuts_range) {
    auto loc_c = c;
    loc_c.apply_theta_cut = false;
    loc_c.espri_selector = Espri::left;
    auto esl_sns = sn_proc(loc_c, ds_type);
    std::get<0>(res).push_back(esl_sns);
    loc_c.espri_selector = Espri::right;
    auto esr_sns = sn_proc(loc_c, ds_type);
    std::get<1>(res).push_back(esr_sns);
  }
  return res;
}

LrSnsRangePair
draw_sn_range_out_phi_bg(Script& script, DsType ds_type,
                         std::vector<Cuts>& cuts_range) {
  return draw_sn_range(script, ds_type, cuts_range, compute_sn_out_phi_bg);
}

LrSnsRangePair
draw_sn_range_carbon_bg(Script& script, DsType ds_type,
                        std::vector<Cuts>& cuts_range) {
  return draw_sn_range(script, ds_type, cuts_range, compute_sn_carbon_bg);
}

Yield compute_he4_cs_out_phi_bg(const Cuts& cuts) {
  assert(cuts.espri_selector != Espri::both);
  auto he4 = load_yields_data(cuts, DsType::he4);
  auto signal = he4->total_yields();
  auto bg = he4->bg_yields();
  auto elastic = extract_out_phi_bg(DsType::he4, cuts,
                                    signal, bg);
  auto he_cm_cs = yields_to_cm_cs(elastic,
                                  cuts, DsType::he4);
  return he_cm_cs;
}

Yield compute_he4_cs_carbon_bg(const Cuts& cuts) {
  assert(cuts.espri_selector != Espri::both);
  auto he4 = load_yields_data(cuts, DsType::he4);
  auto signal = he4->total_yields();
  auto bg = he4->bg_yields();

  auto bg_cuts = cuts;
  bg_cuts.apply_hodf_he4_cut = true;
  auto carbon = load_yields_data(bg_cuts, DsType::carbon);
  auto carbon_bg = carbon->total_yields();
  auto carbon_bg_norm = scale_carbon_bg(DsType::he4, carbon_bg);
  auto signal_wo_bg = signal - carbon_bg_norm;
  auto he_cm_cs = yields_to_cm_cs(signal_wo_bg,
                                  cuts, DsType::he4);
  return he_cm_cs;
}

Yield compute_he6_cs(const Cuts& cuts) {
  assert(cuts.espri_selector != Espri::both);
  auto he6 = load_yields_data(cuts, DsType::he6);
  auto signal = he6->total_yields();
  auto bg = he6->bg_yields();
  auto carbon = load_yields_data(cuts, DsType::carbon);
  auto carbon_bg = carbon->total_yields();
  auto carbon_bg_norm = scale_carbon_bg(DsType::he6, carbon_bg);
  auto signal_wo_bg = signal;
  if (cuts.subtract_he6_carbon_bg) {
    // for (int i = 0; i < carbon_bg_norm.GetBinsNumber(); i++) {
    //   double arg = carbon_bg_norm.GetBinArg(i);
    //   // do not subtract C_bg at FW lab angles due to low stats.
    //   if (arg < 60.3) {
    //     carbon_bg_norm.SetBinValue(arg, 0);
    //     carbon_bg_norm.SetBinError(arg, 0);
    //   }
    // }
    signal_wo_bg -= carbon_bg_norm;
  }

  auto he_cm_cs = yields_to_cm_cs(signal_wo_bg,
                                  cuts, DsType::he6);
  return he_cm_cs;
}

Yield compute_carbon_cs(const Cuts& cuts) {
  assert(cuts.espri_selector != Espri::both);
  auto carbon = load_yields_data(cuts, DsType::carbon);
  auto carbon_bg = carbon->total_yields();
  auto carbon_bg_norm = scale_carbon_bg(DsType::he6, carbon_bg);
  auto carbon_cm_cs = yields_to_cm_cs(carbon_bg_norm,
                                      cuts, DsType::carbon);
  return carbon_cm_cs;
}

template <typename CSProc>
std::vector<std::vector<Yield> >
compute_cses_range(DsType ds_type, std::vector<Cuts>& cuts, CSProc cs_proc) {
  std::vector<std::vector<Yield> > cses_range;
  for (auto& c : cuts) {
    c.espri_selector = Espri::left;
    Yield left_cs = cs_proc(c);
    auto tit = s13::ana::tstrfmt("cs_data/%s_%s_%s_cs.csv",
                                 s13::ana::ds_type_to_str(ds_type).data(),
                                 cli_opts.pol_direction_str().data(),
                                 c.GetTotalCutName().Data());
    write_cs_to_file(tit, left_cs);

    c.espri_selector = Espri::right;
    Yield right_cs = cs_proc(c);
    tit = s13::ana::tstrfmt("cs_data/%s_%s_%s_cs.csv",
                            s13::ana::ds_type_to_str(ds_type).data(),
                            cli_opts.pol_direction_str().data(),
                            c.GetTotalCutName().Data());
    write_cs_to_file(tit, right_cs);

    c.espri_selector = Espri::both;
    // Yield total_cs = (left_cs + right_cs).multiply(1/2.);
    Yield total_cs =
      s13::ana::weighted_average(left_cs, right_cs);
    tit = s13::ana::tstrfmt("cs_data/%s_%s_%s_cs.csv",
                            s13::ana::ds_type_to_str(ds_type).data(),
                            cli_opts.pol_direction_str().data(),
                            c.GetTotalCutName().Data());
    write_cs_to_file(tit, total_cs);

    cses_range.push_back({left_cs, right_cs, total_cs});
  }
  return cses_range;
}

std::vector<std::vector<Yield> >
compute_he4_cses_range_out_phi_bg(std::vector<Cuts>& cuts) {
  return compute_cses_range(DsType::he4, cuts, compute_he4_cs_out_phi_bg);
}

std::vector<std::vector<Yield> >
compute_he4_cses_range_carbon_bg(std::vector<Cuts>& cuts) {
  return compute_cses_range(DsType::he4, cuts, compute_he4_cs_carbon_bg);
}

std::vector<std::vector<Yield> >
compute_he6_cses_range(std::vector<Cuts>& cuts) {
  return compute_cses_range(DsType::he6, cuts, compute_he6_cs);
}

std::vector<std::vector<Yield> >
compute_carbon_cses_range(std::vector<Cuts>& cuts) {
  return compute_cses_range(DsType::carbon, cuts, compute_carbon_cs);
}

using DrawAxesP = s13::ana::NamedType<bool, struct draw_axes_>;

void draw_cses_vector(const std::vector<Yield>& cses,
                      const std::vector<TString>& titles,
                      DrawAxesP draw_axes = DrawAxesP(true),
                      TString draw_opts = "P") {
  // s13::misc::MessageLogger logger("draw_cses_vector()");
  int cr = 1;
  for (size_t i = 0; i < cses.size(); i++) {
    auto cs_graph = s13::ana::
      BuildTGraphErrors(cses[i], titles[i], "CM angle [deg]", "DCS [mb/sr]");
    assert(cs_graph != nullptr && "Failed to construct TGraph");
    cs_graph->SetMarkerColor(cr++);
    if (i == 0 && draw_axes.get() == true) {
      cs_graph->Draw("A " + draw_opts);
    } else {
      cs_graph->Draw("SAME " + draw_opts);
    }
  }
}

std::vector<TGraph*>
make_sn_graphs(std::vector<Yield> sns, std::vector<TString> titles) {
  assert(sns.size() == titles.size());
  int cr = 1;
  std::vector<TGraph*> graphs;
  for (size_t i = 0; i < sns.size(); i++) {
    auto sn_graph = s13::ana::
      BuildTGraphErrors(sns[i], titles[i], "Lab angle [deg]", "S/N [mb/sr]");
    sn_graph->SetMarkerColor(cr);
    sn_graph->SetLineColor(cr++);
    sn_graph->SetMinimum(0.0);
    graphs.push_back(sn_graph);
  }
  return graphs;
}

void draw_sns_range(Script& script,
                    LrSnsRangePair& lr_sns_range,
                    const std::vector<Cuts>& cuts,
                    DsType ds_type) {
  auto left_sns_vector = std::get<0>(lr_sns_range);
  auto right_sns_vector = std::get<1>(lr_sns_range);
  assert(left_sns_vector.size() == right_sns_vector.size());
  assert(left_sns_vector.size() == cuts.size());

  std::vector<Yield> left_sns, right_sns;
  std::vector<TString> left_titles, right_titles;
  auto transform_sn = [ds_type](Yield& y) -> Yield {
    auto y_inv = y.invert();
    return s13::ana::lab_yield_to_cm(y_inv, ds_type);
  };
  for (size_t i = 0; i < left_sns_vector.size(); i++) {
    auto cur_cuts = cuts[i];

    left_sns.push_back(transform_sn(left_sns_vector[i]));
    cur_cuts.espri_selector = Espri::left;
    left_titles.push_back(cur_cuts.GetTotalCutName());

    right_sns.push_back(transform_sn(right_sns_vector[i]));
    cur_cuts.espri_selector = Espri::right;
    right_titles.push_back(cur_cuts.GetTotalCutName());
  }

  script.NewPage(2,1);

  script.cd();
  auto mg_left =
    make_sn_mgraph("ESL 1/(S/N)",
                   make_sn_graphs(left_sns, left_titles));
  mg_left->Draw("AP");
  mg_left->SetMinimum(0.0);
  mg_left->SetMaximum(0.7);
  gPad->BuildLegend(0.50, 0.75, 0.95, 0.9, "");

  script.cd();
  auto mg_right =
    make_sn_mgraph("ESR 1/(S/N)",
                   make_sn_graphs(right_sns, right_titles));
  mg_right->Draw("AP");
  mg_right->SetMinimum(0.0);
  mg_right->SetMaximum(0.7);
  gPad->BuildLegend(0.50, 0.75, 0.95, 0.9, "");
}

void draw_cses_range(Script& script,
                     std::vector<std::vector<Yield> > cses_range,
                     const std::vector<Cuts>& cuts,
                     TMultiGraph* base_mg) {
  s13::misc::MessageLogger logger("draw_cses_range()");
  assert(cuts.size() == cses_range.size());
  std::vector<Yield> left_cses, right_cses, total_cses;
  std::vector<TString> left_titles, right_titles, total_titles;

  logger.debug("separating LRT CSes...");
  for (size_t i = 0; i < cses_range.size(); i++) {
    auto& lrt_cses = cses_range[i];
    auto cur_cuts = cuts[i];
    assert(lrt_cses.size() == 3);

    left_cses.push_back(lrt_cses[0]);
    cur_cuts.espri_selector = Espri::left;
    left_titles.push_back(cur_cuts.GetTotalCutName());

    right_cses.push_back(lrt_cses[1]);
    cur_cuts.espri_selector = Espri::right;
    right_titles.push_back(cur_cuts.GetTotalCutName());

    total_cses.push_back(lrt_cses[2]);
    cur_cuts.espri_selector = Espri::both;
    total_titles.push_back(cur_cuts.GetTotalCutName());
  }

  logger.debug("drawing LR CSes...");
  script.NewPage(2,1);

  script.cd();
  gPad->SetLogy();
  base_mg->Draw("AL");
  draw_cses_vector(left_cses, left_titles, DrawAxesP(false));
  gPad->BuildLegend(0.50, 0.75, 0.95, 1.0, "");

  script.cd();
  gPad->SetLogy();
  base_mg->Draw("AL");
  draw_cses_vector(right_cses, right_titles, DrawAxesP(false));
  gPad->BuildLegend(0.50, 0.75, 0.95, 1.0, "");
}

void draw_cs_range() {
  bool use_he4_carbon_bg = false;
  Script script("draw_cs_range", "cs_data/pdf");
  std::vector<Cuts> cuts =
    s13::io::parse_cuts_range_config(std::fstream(cli_opts.cuts_config()));
  // KLUDGE: won't use theta cut
  // cuts.PrepareVarThetaCut(cli_opts.dataset_type());

  if (cli_opts.dataset_type() == DsType::he4) {
    std::vector<std::vector<Yield> > cses_range;
    LrSnsRangePair lr_sns_range;
    if (use_he4_carbon_bg) {
      lr_sns_range = draw_sn_range_carbon_bg(script, DsType::he4, cuts);
      cses_range = compute_he4_cses_range_carbon_bg(cuts);
      // draw_total_kin_corr_carbon_bg(script, cuts[0]);
    } else {
      lr_sns_range = draw_sn_range_out_phi_bg(script, DsType::he4, cuts);
      cses_range = compute_he4_cses_range_out_phi_bg(cuts);
      draw_total_kin_corr_out_phi_bg(script, cuts[0]);
    }
    auto* moss_mg = make_cs_mgraph("", {s13::ana::
          BuildMossCmsCsGraph("data/he4_200mev_moss_cs.csv")});
    draw_sns_range(script, lr_sns_range, cuts, DsType::he4);
    draw_cses_range(script, cses_range, cuts, moss_mg);
  } else if (cli_opts.dataset_type() == DsType::he6) {
    auto lr_sns_range = draw_sn_range_carbon_bg(script, DsType::he6, cuts);
    auto cses_range = compute_he6_cses_range(cuts);
    auto* mg = make_cs_mgraph("He6 CS", {});
    s13::ana::AddHe6TheorCses(mg);
    draw_sns_range(script, lr_sns_range, cuts, DsType::he6);
    draw_cses_range(script, cses_range, cuts, mg);
  } else if (cli_opts.dataset_type() == DsType::carbon) {
    auto cses_range = compute_carbon_cses_range(cuts);
    auto* mg = make_cs_mgraph("He6-C(p) CS", {});
    s13::ana::AddHe6TheorCses(mg);
    draw_cses_range(script, cses_range, cuts, mg);
  } else {
    throw std::invalid_argument("Invalid argument for CS drawing");
  }
}


int main(int argc, char** argv) {
  s13::misc::MessageLogger logger("main()");
  // s13::ana::SetRootStyle(s13::ana::RootStyle::publication);
  gStyle->SetLegendTextSize(0.023);
  gStyle->SetOptFit();

  try {
    cli_opts.parse_options(argc, argv);
  } catch (std::exception& ex) {
    logger.error("Invalid argument provided: %s", ex.what());
    return -1;
  }

  // script.SetName(cli_opts.output_file_basename("cs_"));
  draw_cs_range();
}
