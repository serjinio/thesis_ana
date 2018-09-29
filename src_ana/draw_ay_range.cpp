/**
 *   \file draw_ay_range.cpp
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
get_ay_out_fname(DsType ds_type,
                 const Cuts& cuts) {
  auto ds_type_str = ds_type == DsType::carbon ? "carbon" :
    cli_opts.dataset_type_str();
  TString fname = s13::ana::tstrfmt("cs_data/%s_%s_ay.csv",
                                    ds_type_str.data(),
                                    cuts.GetTotalCutName().Data());
  return fname.Data();
}

std::shared_ptr<s13::ana::ElasticScatteringYieldsAlg>
load_yields_data(PolDir pol_dir, const Cuts& cuts, DsType ds_type) {
  s13::misc::MessageLogger logger("load_yields_data()");
  auto ds_type_str = ds_type == DsType::carbon ? "carbon" :
    cli_opts.dataset_type_str();
  auto pol_dir_str = ds_type != DsType::carbon ?
    s13::ana::pol_dir_to_str(pol_dir) : "pol_both";
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

Yield
scale_carbon_bg(DsType ds_type, Yield carbon_yield) {
  assert(ds_type != DsType::carbon && "Only he4 & he4 datasets allowed here!");
  s13::misc::MessageLogger logger("scale_carbon_bg()");
  double fact = ds_type == DsType::he6 ?
    s13::gk_carbon_bg_norm_factor : s13::gk_carbon_to_he4_bg_norm_factor;
  logger.info("Will use carbon BG scaling factor for %s dataset",
              ds_type == DsType::he6 ? "He6" : "He4");
  logger.info("Carbon data scaling factor is reduced by 2 as we compute "
              "asymmetry (for up/down runs).");
  fact /= 2;
  logger.info("Scaling factor value: %.2f.", fact);
  auto scaled_yields = carbon_yield.multiply(fact);
  return scaled_yields;
}

Yield
extract_out_phi_bg(DsType ds_type, const Cuts& cuts, Yield total_yields, Yield bg_yields) {
  auto bg_in_yields = s13::ana::convert_to_in_phi_bg(ds_type, cuts, bg_yields);
  return total_yields - bg_in_yields;
}

TMultiGraph* make_ay_mgraph(TString title, std::vector<TGraph*> graphs) {
  TMultiGraph* mg = new TMultiGraph();
  mg->SetTitle(TString::Format("%s;%s;%s", title.Data(),
                               "CM angle [deg.]", "Ay x Pt"));
  for (auto* el : graphs) {
    mg->Add(el);
  }
  return mg;
}

void write_yields_to_file(TString yields_fname, Yield cs) {
  std::ofstream ofs(yields_fname);
  cs.Serialize(ofs);
}

Yield normalize_ela_yields(DsType ds_type, PolDir pol_dir,
                           Yield elastic_yields, const Cuts& cuts) {
  s13::misc::MessageLogger logger("normalize_ela_yields()");
  logger.debug("Will normalize yields...");
  double num_incident = s13::ana::get_num_beam(ds_type, pol_dir, cuts.vertexXY_radius);
  std::string pol_dir_str = s13::ana::pol_dir_to_str(pol_dir);
  logger.debug("Number of beam particles for %s: %.6f.",
               pol_dir_str.data(), num_incident);
  auto norm_yield = elastic_yields.multiply(1/num_incident);
  // logger.debug("Yield before normalization:");
  // std::cout << elastic_yields;
  // logger.debug("Yield after normalization:");
  // std::cout << norm_yield;
  return norm_yield;
}

Yield compute_he4_ela_yields_out_phi_bg(PolDir pol_dir, const Cuts& cuts) {
  assert(cuts.espri_selector != Espri::both);
  auto he4 = load_yields_data(pol_dir, cuts, DsType::he4);
  auto signal = he4->total_yields();
  auto bg = he4->bg_yields();
  auto elastic = extract_out_phi_bg(DsType::he4, cuts,
                                    signal, bg);
  auto elastic_cm = s13::ana::lab_yield_to_cm(elastic, DsType::he4);
  return elastic_cm;
}

Yield compute_he4_ela_yields_carbon_bg(PolDir pol_dir, const Cuts& cuts) {
  assert(cuts.espri_selector != Espri::both);
  auto he4 = load_yields_data(pol_dir, cuts, DsType::he4);
  auto signal = he4->total_yields();
  auto bg = he4->bg_yields();

  auto bg_cuts = cuts;
  bg_cuts.apply_hodf_he4_cut = true;
  auto carbon = load_yields_data(pol_dir, bg_cuts, DsType::carbon);
  auto carbon_bg = carbon->total_yields();
  auto carbon_bg_norm = scale_carbon_bg(DsType::he4, carbon_bg);
  auto elastic = signal - carbon_bg_norm;
  auto elastic_cm = s13::ana::lab_yield_to_cm(elastic, DsType::he4);
  return elastic_cm;
}

Yield compute_he6_ela_yields(PolDir pol_dir, const Cuts& cuts) {
  assert(cuts.espri_selector != Espri::both);
  auto he6 = load_yields_data(pol_dir, cuts, DsType::he6);
  auto signal = he6->total_yields();
  auto bg = he6->bg_yields();
  auto carbon = load_yields_data(pol_dir, cuts, DsType::carbon);
  auto carbon_bg = carbon->total_yields();
  auto carbon_bg_norm = scale_carbon_bg(DsType::he6, carbon_bg);
  auto elastic = signal;
  if (cuts.subtract_he6_carbon_bg) {
    for (int i = 0; i < carbon_bg_norm.GetBinsNumber(); i++) {
      double arg = carbon_bg_norm.GetBinArg(i);
      // do not subtract C_bg at FW lab angles due to low stats.
      if (arg < 60.3) {
        carbon_bg_norm.SetBinValue(arg, 0);
        carbon_bg_norm.SetBinError(arg, 0);
      }
    }
    elastic -= carbon_bg_norm;
  }

  auto elastic_cm = s13::ana::lab_yield_to_cm(elastic, DsType::he4);
  return elastic_cm;
}

template <typename ComputeProc>
std::vector<std::vector<Yield> >
compute_yields_range(PolDir pol_dir, std::vector<Cuts>& cuts, ComputeProc compute_yields) {
  std::vector<std::vector<Yield> > yields_range;
  std::string pol_dir_str = s13::ana::pol_dir_to_str(pol_dir);
  for (auto& c : cuts) {
    c.espri_selector = Espri::left;

    Yield left_yield = compute_yields(pol_dir, c);
    auto tit = s13::ana::tstrfmt("cs_data/%s_%s_ela_yield.csv",
                                 pol_dir_str.data(),
                                 c.GetTotalCutName().Data());
    write_yields_to_file(tit, left_yield);

    c.espri_selector = Espri::right;
    Yield right_yield = compute_yields(pol_dir, c);
    tit = s13::ana::tstrfmt("cs_data/%s_%s_ela_yield.csv",
                            pol_dir_str.data(),
                            c.GetTotalCutName().Data());
    write_yields_to_file(tit, right_yield);

    yields_range.push_back({left_yield, right_yield});
  }
  return yields_range;
}

std::vector<std::vector<Yield> >
compute_he4_yields_range_out_phi_bg(PolDir pol_dir, std::vector<Cuts>& cuts) {
  return compute_yields_range(pol_dir, cuts, compute_he4_ela_yields_out_phi_bg);
}

std::vector<std::vector<Yield> >
compute_he4_yields_range_carbon_bg(PolDir pol_dir, std::vector<Cuts>& cuts) {
  return compute_yields_range(pol_dir, cuts, compute_he4_ela_yields_carbon_bg);
}

std::vector<std::vector<Yield> >
compute_he6_yields_range(PolDir pol_dir, std::vector<Cuts>& cuts) {
  return compute_yields_range(pol_dir, cuts, compute_he6_ela_yields);
}

Yield
compute_ay(Yield left_pol_up, Yield right_pol_up,
           Yield left_pol_down, Yield right_pol_down) {
  auto left = (left_pol_up * right_pol_down).sqrt();
  auto right = (right_pol_up * left_pol_down).sqrt();
  auto ay = (left - right) / (left + right);

  // std::cout << "lu: " << left_pol_up;
  // std::cout << "ru: " << right_pol_up;
  // std::cout << "ld: " << left_pol_down;
  // std::cout << "rd: " << right_pol_down;

  // orig err
  auto num = (left_pol_up * right_pol_down * right_pol_up * left_pol_down).sqrt();
  auto den = ((left_pol_up * right_pol_down).sqrt() +
              (right_pol_up * left_pol_down).sqrt()).pow(2);
  auto fact = (left_pol_up.invert() + right_pol_up.invert() +
               left_pol_down.invert() + right_pol_down.invert()).sqrt();
  auto err = (num / den) * fact;

  // my err
  // auto num = (( left_pol_up * right_pol_down ) + (right_pol_up * left_pol_down)).sqrt();
  // auto den = ((left_pol_up * right_pol_down).sqrt() +
  //             (right_pol_up * left_pol_down).sqrt()).pow(2);
  // auto fact = (right_pol_down + left_pol_up + left_pol_down + right_pol_up).sqrt();
  // fact = fact.multiply(1/std::sqrt(2));
  // auto err = (num / den) * fact;

  for (int i = 0; i < ay.GetBinsNumber(); i++) {
    double proton_cm_theta = ay.GetBinArg(i);
    ay.SetBinError(proton_cm_theta, err.GetBinValue(i));
    // if (proton_cm_theta > 45 && proton_cm_theta < 50) {
    //   ay.SetBinValue(proton_cm_theta, -0.25);
    //   ay.SetBinError(proton_cm_theta, 0.2);
    // }
  }

  std::cout << "Ay: " << ay;
  return ay;
}

std::vector<Yield>
compute_ays_range(DsType ds_type,
                  std::vector<std::vector<Yield> >& yields_range,
                  std::vector<std::vector<Yield> >& pol_down_yields_range,
                  std::vector<Cuts> cuts_range) {
  assert(yields_range.size() == pol_down_yields_range.size());
  std::vector<Yield> ays_range;
  for (size_t i = 0; i < yields_range.size(); i++) {
    auto& yield_left_pol_up = yields_range[i][0];
    auto& yield_right_pol_up = yields_range[i][1];
    auto& yield_left_pol_down = pol_down_yields_range[i][0];
    auto& yield_right_pol_down = pol_down_yields_range[i][1];
    auto ay = compute_ay(yield_left_pol_up, yield_right_pol_up,
                         yield_left_pol_down, yield_right_pol_down);
    std::string out_fname = get_ay_out_fname(ds_type, cuts_range[i]);
    auto ofs = std::ofstream(out_fname);
    ay.Serialize(ofs);
    ays_range.push_back(ay);
  }
  return ays_range;
}

double
fit_he4_ay_with_moss(Yield& ay, TGraph* ay_graph) {
  auto moss_ay_data = s13::io::load_csv("data/he4_200mev_moss_ay.csv", 2, 2);
  auto moss_ay_interp =
    new s13::ana::TInterp(moss_ay_data[0].size(), ROOT::Math::Interpolation::kLINEAR);
  moss_ay_interp->SetData(moss_ay_data[0], moss_ay_data[1]);
  s13::ana::TInterpPtr moss_ay_interp_ptr{moss_ay_interp};
  auto fit_fn = [moss_ay_interp_ptr](Double_t* x, Double_t* params) {
    Double_t fit_val = moss_ay_interp_ptr->Eval(x[0]) * params[0];
    return fit_val;
  };

  TF1* fit_tfn = new TF1("ay_fit", fit_fn, 33, 62, 1);
  fit_tfn->SetParameters(0.12, 0.0);
  fit_tfn->SetParNames("Pt");

  ay_graph->Fit("ay_fit", "r");
  double* fit_fn_params = fit_tfn->GetParameters();
  return fit_fn_params[0];
}

std::vector<double>
fit_he4_ays_range_with_moss(std::vector<Yield> ays_range,
                            std::vector<TGraph*> ay_graphs) {
  std::vector<double> tgt_pol_consts;
  for (size_t i = 0; i < ays_range.size(); i++) {
    double fit_const = fit_he4_ay_with_moss(ays_range[i], ay_graphs[i]);
    tgt_pol_consts.push_back(fit_const);
  }
  return tgt_pol_consts;
}

using DrawAxesP = s13::ana::NamedType<bool, struct draw_axes_>;

std::vector<TGraph*>
draw_ays_vector(const std::vector<Yield>& ays,
                const std::vector<TString>& titles,
                DrawAxesP draw_axes = DrawAxesP(true),
                double ymin = -1, double ymax = 1,
                TString draw_opts = "P") {
  // s13::misc::MessageLogger logger("draw_ays_vector()");
  std::vector<TGraph*> graphs;
  int cr = 1;
  for (size_t i = 0; i < ays.size(); i++) {
    auto ay_graph = s13::ana::
      BuildTGraphErrors(ays[i], titles[i], "CM angle [deg]", "Ay");
    assert(ay_graph != nullptr && "Failed to construct TGraph");
    ay_graph->SetMarkerColor(cr++);
    if (i == 0 && draw_axes.get() == true) {
      ay_graph->Draw("A " + draw_opts);
    } else {
      ay_graph->Draw("SAME " + draw_opts);
    }
    ay_graph->SetMaximum(ymax);
    ay_graph->SetMinimum(ymin);
    graphs.push_back(ay_graph);
  }
  return graphs;
}

std::vector<TGraph*>
draw_ays_range(Script& script,
               std::vector<Yield> ays_range,
               const std::vector<Cuts>& cuts,
               TMultiGraph* base_mg,
               double ymax = 1) {
  s13::misc::MessageLogger logger("draw_ays_range()");
  assert(cuts.size() == ays_range.size());
  std::vector<TString> ay_titles;
  for (auto& c : cuts) {
    TString tit = s13::ana::tstrfmt("%s_ay", c.GetTotalCutName().Data());
    ay_titles.push_back(tit);
  }

  logger.debug("drawing Ays...");
  script.NewPage(1,1);

  script.cd();
  // base_mg->Draw("AP");
  auto graphs = draw_ays_vector(ays_range, ay_titles, DrawAxesP(true),
                                -ymax, ymax);
  gPad->BuildLegend(0.1, 0.80, 0.6, 0.90, "");
  return graphs;
}

void draw_ay_range() {
  bool use_he4_carbon_bg = false;
  Script script("draw_ay_range", "cs_data/pdf");
  std::vector<Cuts> cuts =
    s13::io::parse_cuts_range_config(std::fstream(cli_opts.cuts_config()));
  // KLUDGE: won't use vartheta cut
  // cuts.PrepareVarThetaCut(cli_opts.dataset_type());

  if (cli_opts.dataset_type() == DsType::he4) {
    std::vector<Yield> ays_range;
    if (use_he4_carbon_bg) {
      // TODO: Maybe implement Ay for he4 - C_bg
      throw std::invalid_argument("Not supported yet!");

      // lr_sns_range = draw_sn_range_carbon_bg(script, DsType::he4, cuts);
      // pol_up_yields_range = compute_he4_yields_range_carbon_bg(cuts);
      // draw_total_kin_corr_carbon_bg(script, cuts[0]);
    } else {
      auto pol_up_yields_range = compute_he4_yields_range_out_phi_bg(PolDir::up, cuts);
      auto pol_down_yields_range = compute_he4_yields_range_out_phi_bg(PolDir::down, cuts);
      ays_range = compute_ays_range(cli_opts.dataset_type(),
                                    pol_up_yields_range, pol_down_yields_range, cuts);
    }

    auto* theor_mg = make_ay_mgraph("", {s13::ana::
          BuildMossAyGraph("data/he4_200mev_moss_ay.csv")});
    auto ay_graphs = draw_ays_range(script, ays_range, cuts, theor_mg, 0.25);
    auto* ays_mg = make_ay_mgraph("Asymmetries & Moss fit (red)", ay_graphs);
    auto target_pol_consts = fit_he4_ays_range_with_moss(ays_range, ay_graphs);
    std::vector<Yield> ays_fitted_range;
    for (size_t i = 0; i < ays_range.size(); i++) {
      auto title = s13::ana::tstrfmt("%s (P_tgt: %.3f)",
                                     ay_graphs[i]->GetTitle(),
                                     target_pol_consts[i]);
      ay_graphs[i]->SetTitle(title);
    }
    ays_mg->Draw("AP");
    ays_mg->SetMinimum(-0.5);
    ays_mg->SetMaximum(0.5);
    gPad->BuildLegend(0.13, 0.13, 0.60, 0.20, "");
    // draw_ays_range(script, ays_fitted_range, cuts, theor_mg);
  } else if (cli_opts.dataset_type() == DsType::he6) {
    auto pol_up_yields_range = compute_he6_yields_range(PolDir::up, cuts);
    auto pol_down_yields_range = compute_he6_yields_range(PolDir::down, cuts);
    auto ays_range = compute_ays_range(cli_opts.dataset_type(),
                                       pol_up_yields_range, pol_down_yields_range, cuts);
    auto* mg = make_ay_mgraph("He6 CS", {});
    // s13::ana::AddHe6TheorAys(mg);
    draw_ays_range(script, ays_range, cuts, mg, 0.2);
  } else {
    throw std::invalid_argument("Invalid argument for CS drawing");
  }
}


int main(int argc, char** argv) {
  s13::misc::MessageLogger logger("main()");
  // s13::ana::SetRootStyle(s13::ana::RootStyle::publication);
  gStyle->SetLegendTextSize(0.023);
  gStyle->SetOptFit(1111);

  try {
    cli_opts.parse_options(argc, argv);
  } catch (std::exception& ex) {
    logger.error("Invalid argument provided: %s", ex.what());
    return -1;
  }

  // script.SetName(cli_opts.output_file_basename("cs_"));
  draw_ay_range();
}
