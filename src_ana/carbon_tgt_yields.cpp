
#include <functional>

#include "TROOT.h"

#include "common.hpp"
#include "cli.hpp"
#include "init_ds.hpp"
#include "consts.hpp"
#include "drawing.hpp"
#include "scattyield.hpp"
#include "bgsubtract.hpp"
#include "treewalk.hpp"
#include "cuts_conf.hpp"
#include "rootscript.hpp"
#include "csutil.hpp"


static s13::misc::CliOptions cli_opts;

using Script = s13::ana::RootScript;
using Yield = s13::ana::ScatteringYield;
using DsType = s13::ana::DatasetType;
using PolDir = s13::ana::PolarizationDirection;
using Dsdcs = s13::ana::Dsdcs;
using ElaCuts = s13::ana::ElasticScatteringCuts;
using Evt = s13::ana::ScatteringEvent;
using FitParams = s13::ana::FitParams;
using Espri = s13::ana::Espri;
using LinearInterp = s13::ana::LinearInterp;


class TgtYieldsPolAngle : public s13::ana::SimpleTTreeAlgorithmBase {

public:

  static constexpr double k_min_theta = s13::gk_cs_range_start;
  static constexpr double k_max_theta = s13::gk_cs_range_end;

  TH1* tgt_yields_hist;
  TH2* tgt_yields_2d_hist;

private:

  s13::misc::MessageLogger logger_;

  ElaCuts cuts_;
  DsType ds_type_;

  size_t bin_num_;

  s13::ana::TInterpPtr he4_kin_interp_;
  s13::ana::TInterpPtr lab2cm_angle_interp_;

public:

  TgtYieldsPolAngle(const ElaCuts& cuts, const DsType ds_type, size_t bin_num = 15) :
    logger_{"TgtYieldsPolAngle"}, cuts_{cuts}, ds_type_{ds_type}, bin_num_{bin_num} {
    init_hists();
    init_kin_interp();
  }

  TgtYieldsPolAngle(const TgtYieldsPolAngle&) = delete;
  TgtYieldsPolAngle(TgtYieldsPolAngle&&) = delete;
  TgtYieldsPolAngle operator=(const TgtYieldsPolAngle&) = delete;
  TgtYieldsPolAngle operator=(TgtYieldsPolAngle&&) = delete;

  virtual void setup() {
  }

  virtual void process_event(Evt& evt) {
    fill_hists(evt);
  }

  virtual void finalize() {
  }

  const ElaCuts& cuts() {
    return cuts_;
  }

  double bin_width() {
    return (k_max_theta - k_min_theta) / bin_num_;
  }

  double bin_min_theta(int idx) {
    double minth = k_min_theta + bin_width() * idx;
    return minth;
  }

  double bin_max_theta(int idx) {
    double maxth = k_min_theta + bin_width() * (idx + 1);
    return maxth;
  }

  double bin_center_theta(int idx) {
    return (bin_min_theta(idx) + bin_max_theta(idx)) / 2;
  }

private:

  TH1* init_kin_corr_hist(TString title) {
    auto th1 = &s13::ana::make_th1;
    auto hist = th1(80, -12, 12, title, "#theta_{exp} - #theta_{theor}", "Counts");
    return hist;
  }

  TH1* init_pol_angle_corr_hist(TString title) {
    auto th1 = &s13::ana::make_th1;
    auto hist = th1(16, 56, 71, title, "#theta_{p} (deg.)", "Counts");
    return hist;
  }

  TH2* init_2d_kin_corr_hist(TString title) {
    auto th2 = &s13::ana::make_th2;
    auto hist = th2(200, 0, 12, 200, 55, 70, title, "#theta_{He}", "#theta_{p}");
    return hist;
  }

  int bin_idx(Evt& evt) {
    if (evt.p_theta_eff - k_min_theta < 0) {
      return -1;
    }

    size_t idx = (evt.p_theta_eff - k_min_theta) / bin_width();
    if (idx < bin_num_) {
      return idx;
    } else {
      return -1;
    }
  }

  void init_hists() {
    std::string title = "Pol. target (C_{10}H_{8})";
    if (ds_type_ == DsType::carbon) {
      title = "C target";
    }
    tgt_yields_hist = init_pol_angle_corr_hist(title);
    tgt_yields_2d_hist = init_2d_kin_corr_hist(title);
  }

  void init_kin_interp() {
    he4_kin_interp_ = s13::ana::LoadKinematicsData("data/he4_kin.csv");
    lab2cm_angle_interp_ = make_lab2cm_angle_interp(ds_type_);
  }

  void fill_hists(Evt& evt) {
    Double_t he_theta = cuts_.dsdc_selector == Dsdcs::s1dc
      ? evt.s1dc_theta : evt.fdc0_theta;
    double p_theta = evt.p_theta_eff;

    if (!is_well_defined(p_theta)
        || !is_well_defined(he_theta) ||
        !evt.IsEspriEvent()) {
      return;
    }

    double he_theta_diff = he_theta - he4_kin_interp_->Eval(p_theta);

    if (cuts_.IsSignalEvent(evt)) {
      double p_theta_cm = lab2cm_angle_interp_->Eval(evt.p_theta);
      // logger_.info("Got event with lab angle %.2f; CM angle: %.2f",
      // p_theta, p_theta_cm);
      tgt_yields_hist->Fill(p_theta);
    }
  }

};


void
draw_hists(s13::ana::RootScript& script,
           std::shared_ptr<TgtYieldsPolAngle> alg_signal,
           std::shared_ptr<TgtYieldsPolAngle> alg_bg) {
  s13::misc::MessageLogger logger("draw_hists()");
  logger.info("beginning to draw hists...");
  script.NewPage(1).cd();
  s13::ana::draw_thstack({
      alg_signal->tgt_yields_hist,
      alg_bg->tgt_yields_hist}, "", "NOSTACK");
  // s13::ana::draw_text(alg_signal->cuts().GetTotalCutTitle());
  gPad->BuildLegend(0.3, 0.75, 0.65, 0.86);
  // gPad->SetLogy();
  logger.info("function finished");
}

std::vector<TChain*>
init_carbon_dataset(int num_of_runs = -1) {
  if (num_of_runs == -1) {
    return init_carbon_segmented_dataset(1, s13::gk_carbon_num_of_runs);
  }

  return init_carbon_segmented_dataset(1, num_of_runs);
}

void carbon_tgt_yields() {
  auto cuts = s13::io::parse_cuts_config(std::fstream(cli_opts.cuts_config()));
  auto alg_carbon =
    s13::ana::walk_alg2<TgtYieldsPolAngle>(init_carbon_dataset(cli_opts.dataset_size()),
                                           cli_opts.use_mt(),
                                           cuts,
                                           DsType::carbon);
  auto alg =
    s13::ana::walk_alg2<TgtYieldsPolAngle>(init_dataset_from_cli(cli_opts),
                                           cli_opts.use_mt(),
                                           cuts,
                                           cli_opts.dataset_type());

  auto script_fname =
    s13::ana::tstrfmt("%s_%s_carbon_tgt_yields",
                      cli_opts.output_file_basename().data(),
                      s13::ana::espri_selector_to_str(cuts.espri_selector).data());
  s13::ana::RootScript script(script_fname);
  gStyle->SetOptStat(1111111);
  // gStyle->SetOptFit();
  draw_hists(script, alg, alg_carbon);
  script.Close();
}


int main(int argc, char** argv) {
  s13::misc::MessageLogger log("main()");
  s13::ana::SetRootStyle(s13::ana::RootStyle::publication);

  TThread::Initialize();
  cli_opts.require_cuts_conf_opt(true);
  try {
    cli_opts.parse_options(argc, argv);
  } catch (std::exception& ex) {
    log.error("Invalid argument provided: %s", ex.what());
    return -1;
  }

  carbon_tgt_yields();
}
