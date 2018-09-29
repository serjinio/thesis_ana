

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


// bins number to use in computtions
const int k_bin_num = 15;


class He4CarbonBgYields : public s13::ana::SimpleTTreeAlgorithmBase {

public:

  static constexpr double k_min_theta = s13::gk_cs_range_start;
  static constexpr double k_max_theta = s13::gk_cs_range_end;

  std::vector<TH1*> sn_in_hists;
  std::vector<TH1*> sn_out_hists;
  std::vector<TH1*> sn_tot_hists;
  TH2* in_kin_corr_hist;
  TH2* out_kin_corr_hist;
  TH1* sn_in_hist;
  TH1* sn_out_hist;

private:

  s13::misc::MessageLogger logger_;

  ElaCuts cuts_;

  size_t bin_num_;

  s13::ana::TInterpPtr he4_kin_interp_;

public:

  He4CarbonBgYields(const ElaCuts& cuts, size_t bin_num = 15) :
    logger_{"He4CarbonBgYields"}, cuts_{cuts}, bin_num_{bin_num} {
    init_hists();
    init_kin_interp();
  }

  He4CarbonBgYields(const He4CarbonBgYields&) = delete;
  He4CarbonBgYields(He4CarbonBgYields&&) = delete;
  He4CarbonBgYields operator=(const He4CarbonBgYields&) = delete;
  He4CarbonBgYields operator=(He4CarbonBgYields&&) = delete;

  virtual void setup() {
  }

  virtual void process_event(Evt& evt) {
    fill_hists(evt);
  }

  virtual void finalize() {
    for (size_t i = 0; i < sn_in_hists.size(); i++) {
      // TODO: compute more precisely for p-4He vs. C(p)-(4He + 2n)
      sn_in_hists[i]->Scale(s13::gk_carbon_to_he4_bg_norm_factor);
    }
    sn_in_hist->Scale(s13::gk_carbon_to_he4_bg_norm_factor);
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
    for (size_t i = 0; i < bin_num_; i++) {
      auto ang_range = s13::ana::tstrfmt("(%.2f<p_{#theta}<%.2f)",
                                         bin_min_theta(i), bin_max_theta(i));
      sn_in_hists.push_back(init_kin_corr_hist("IN C+4He " + ang_range));
      sn_out_hists.push_back(init_kin_corr_hist("OUT C+4He " + ang_range));
      sn_tot_hists.push_back(init_kin_corr_hist("TOT C+4He " + ang_range));
    }

    in_kin_corr_hist = init_2d_kin_corr_hist("IN C+4He");
    out_kin_corr_hist = init_2d_kin_corr_hist("OUT C+4He");
    sn_in_hist = init_kin_corr_hist("IN C+4He");
    sn_out_hist = init_kin_corr_hist("OUT C+4He");
  }

  void init_kin_interp() {
    he4_kin_interp_ = s13::ana::LoadKinematicsData("data/he4_kin.csv");
  }

  void fill_hists(Evt& evt) {
    if (!cuts_.IsHodfHe4Cut(evt)) {
      return;
    }
    int evt_bin_idx =  bin_idx(evt);
    if (evt_bin_idx == -1) {
      return;
    }

    Double_t he_theta = cuts_.dsdc_selector == Dsdcs::s1dc
      ? evt.s1dc_theta : evt.fdc0_theta;
    double p_theta = evt.p_theta_eff;

    if (!is_well_defined(p_theta)
        || !is_well_defined(he_theta) ||
        !evt.IsEspriEvent()) {
      return;
    }

    double he_theta_diff = he_theta - he4_kin_interp_->Eval(p_theta);

    sn_tot_hists[evt_bin_idx]->Fill(he_theta_diff);
    if (cuts_.IsSignalEvent(evt)) {
      in_kin_corr_hist->Fill(he_theta, p_theta);
      sn_in_hist->Fill(he_theta_diff);
      sn_in_hists[evt_bin_idx]->Fill(he_theta_diff);
    } else if (cuts_.IsBackgroundEvent(evt)){
      out_kin_corr_hist->Fill(he_theta, p_theta);
      sn_out_hist->Fill(he_theta_diff);
      sn_out_hists[evt_bin_idx]->Fill(he_theta_diff);
    }
  }

};


class He4Yields : public s13::ana::SimpleTTreeAlgorithmBase {

public:

  static constexpr double k_min_theta = s13::gk_cs_range_start;
  static constexpr double k_max_theta = s13::gk_cs_range_end;

  std::vector<TH1*> sn_in_hists;
  std::vector<TH1*> sn_out_hists;
  std::vector<TH1*> sn_tot_hists;
  TH2* in_kin_corr_hist;
  TH2* out_kin_corr_hist;
  TH1* sn_in_hist;
  TH1* sn_out_hist;

private:

  s13::misc::MessageLogger logger_;

  ElaCuts cuts_;

  size_t bin_num_;

  s13::ana::TInterpPtr he4_kin_interp_;

public:

  He4Yields(const ElaCuts& cuts, size_t bin_num = 15) :
    logger_{"He4Yields"}, cuts_{cuts}, bin_num_{bin_num} {
    init_hists();
    init_kin_interp();
  }

  He4Yields(const He4Yields&) = delete;
  He4Yields(He4Yields&&) = delete;
  He4Yields operator=(const He4Yields&) = delete;
  He4Yields operator=(He4Yields&&) = delete;

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
    for (size_t i = 0; i < bin_num_; i++) {
      auto ang_range = s13::ana::tstrfmt("(%.2f<p_{#theta}<%.2f)",
                                         bin_min_theta(i), bin_max_theta(i));
      sn_in_hists.push_back(init_kin_corr_hist("IN p+4He " + ang_range));
      sn_out_hists.push_back(init_kin_corr_hist("OUT p+4He " + ang_range));
      sn_tot_hists.push_back(init_kin_corr_hist("TOT p+4He " + ang_range));
    }

    in_kin_corr_hist = init_2d_kin_corr_hist("IN p+4he");
    out_kin_corr_hist = init_2d_kin_corr_hist("OUT p+4he.");
    sn_in_hist = init_kin_corr_hist("IN p+4He");
    sn_out_hist = init_kin_corr_hist("OUT p+4He");
  }

  void init_kin_interp() {
    he4_kin_interp_ = s13::ana::LoadKinematicsData("data/he4_kin.csv");
  }

  void fill_hists(Evt& evt) {
    int evt_bin_idx =  bin_idx(evt);
    if (evt_bin_idx == -1) {
      return;
    }

    Double_t he_theta = cuts_.dsdc_selector == Dsdcs::s1dc
      ? evt.s1dc_theta : evt.fdc0_theta;
    double p_theta = evt.p_theta_eff;

    if (!is_well_defined(p_theta)
        || !is_well_defined(he_theta) ||
        !evt.IsEspriEvent()) {
      return;
    }

    double he_theta_diff = he_theta - evt.he_theta_theor;

    sn_tot_hists[evt_bin_idx]->Fill(he_theta_diff);
    if (cuts_.IsSignalEvent(evt)) {
      in_kin_corr_hist->Fill(he_theta, p_theta);
      sn_in_hist->Fill(he_theta_diff);
      sn_in_hists[evt_bin_idx]->Fill(he_theta_diff);
    } else if (cuts_.IsBackgroundEvent(evt)){
      out_kin_corr_hist->Fill(he_theta, p_theta);
      sn_out_hist->Fill(he_theta_diff);
      sn_out_hists[evt_bin_idx]->Fill(he_theta_diff);
    }
  }

};


void draw_cut_marker_lines(double cut_width, TH1* hist) {
  s13::ana::draw_marker_line_ysym(cut_width / 2,
                                  hist->GetMinimum(),
                                  hist->GetMaximum());
}

using SnsVector = std::vector< std::pair<double, double> >;

void write_sns(std::string fname, const SnsVector& sns) {
  std::ofstream ofs(fname);
  ofs << "#p_theta,sn" << std::endl;
  for (size_t i = 0; i < sns.size(); ++i) {
    double theta = std::get<0>(sns[i]);
    double sn = std::get<1>(sns[i]);
    ofs << theta << "," << sn << std::endl;
  }
}

TH1* scale_out_phi_events(TH1* bg_hist, double proton_angle, const ElaCuts& cuts) {
  auto bg_norm_fact_interp =
    s13::ana::make_bg_norm_factor_interp2(DsType::he4,
                                          cuts.dsdc_selector);
  auto h = static_cast<TH1*>(bg_hist->Clone());
  auto tit = s13::ana::tstrfmt("%s (scaled)", bg_hist->GetTitle());
  h->SetTitle(tit);
  double r = bg_norm_fact_interp.eval(proton_angle,
                                      cuts.GetPhiCutWidth(proton_angle));
  h->Scale(r - 1);

  return h;
}

void
draw_hists(s13::ana::RootScript& script,
           std::shared_ptr<He4Yields> alg,
           std::shared_ptr<He4CarbonBgYields> alg_carbon) {
  auto kin_graph = s13::ana::BuildHe4ProtonKinematicsGraph();
  script.NewPage(1,1).cd();
  alg->in_kin_corr_hist->Draw("colz");
  kin_graph->Draw("SAME");
  script.NewPage(1,1).cd();
  alg_carbon->in_kin_corr_hist->Draw("colz");
  kin_graph->Draw("SAME");

  script.NewPage(1,1).cd();
  auto sn_wo_bg = static_cast<TH1*>(alg->sn_in_hist->Clone());
  sn_wo_bg->SetTitle(s13::ana::tstrfmt("%s w/o C BG", sn_wo_bg->GetTitle()));
  sn_wo_bg->Add(alg_carbon->sn_in_hist, -1);
  auto s1 = s13::ana::draw_thstack({
      alg->sn_in_hist, alg_carbon->sn_in_hist, sn_wo_bg},
    "S: p+(4He), N_{scaled}: C(p)+(4He+2n)", "NOSTACK");
  s1->SetMaximum(alg->sn_in_hist->GetMaximum() * 1.5);
  gPad->BuildLegend(0.15, 0.75, 0.60, 0.90);


  SnsVector sns;
  for (size_t i = 0; i < alg->sn_tot_hists.size(); ++i) {
    if (i % 3 == 0) {
      script.NewPage(3, 3).cd(s13::ana::PadSeq::column);
    }

    auto out_phi_bg_scaled =
      scale_out_phi_events(alg->sn_out_hists[i],
                           alg->bin_center_theta(i), alg->cuts());
    out_phi_bg_scaled->SetTitle("out #phi BG");
    auto he4_wo_out_phi_bg = static_cast<TH1*>(alg->sn_in_hists[i]->Clone());
    he4_wo_out_phi_bg->SetTitle(s13::ana::tstrfmt("%s w/o out #phi BG",
                                         he4_wo_out_phi_bg->GetTitle()));
    he4_wo_out_phi_bg->Add(out_phi_bg_scaled, -1);
    auto he4_wo_carbon_bg = static_cast<TH1*>(alg->sn_in_hists[i]->Clone());
    he4_wo_carbon_bg->SetTitle(s13::ana::tstrfmt("%s w/o C BG",
                                          he4_wo_carbon_bg->GetTitle()));
    he4_wo_carbon_bg->Add(alg_carbon->sn_in_hists[i], -1);

    script.cd();
    auto s1 = s13::ana::draw_thstack({
        alg->sn_in_hists[i], alg_carbon->sn_in_hists[i], out_phi_bg_scaled},
      "S: p+(4He), N_{scaled}: C(p)+(4He+2n)", "NOSTACK");
    s1->SetMaximum(alg->sn_in_hists[i]->GetMaximum() * 1.5);
    gPad->BuildLegend(0.15, 0.65, 0.70, 0.85);

    script.cd();
    he4_wo_carbon_bg->Draw();
    draw_cut_marker_lines(alg->cuts().theta_width, alg->sn_in_hists[i]);
    double sn =
      s13::ana::draw_kin_corr_hist_sn(alg->bin_center_theta(i),
                                      he4_wo_carbon_bg,
                                      alg->cuts().dsdc_selector,
                                      alg->cuts().theta_width,
                                      s13::ana::UseAbsVal(false),
                                      s13::ana::Degrees(-11));
    sns.push_back(std::pair<double, double>(alg->bin_center_theta(i), sn));

    script.cd();
    he4_wo_out_phi_bg->Draw();
    draw_cut_marker_lines(alg->cuts().theta_width, alg->sn_in_hists[i]);
    s13::ana::draw_kin_corr_hist_sn(alg->bin_center_theta(i),
                                    he4_wo_out_phi_bg,
                                    alg->cuts().dsdc_selector,
                                    alg->cuts().theta_width,
                                    s13::ana::UseAbsVal(false),
                                    s13::ana::Degrees(-11));
  }

  auto espri_str =
    s13::ana::espri_selector_to_str(alg->cuts().espri_selector).data();
  auto fname =
    s13::ana::tstrfmt("out/%s_%s_he4_carbon_bg_subtraction_sns.csv",
                      cli_opts.output_file_basename().data(),
                      espri_str);
  write_sns(fname.Data(), sns);
}

std::shared_ptr<He4CarbonBgYields>
compute_carbon_yields(ElaCuts& cuts) {
  auto num_runs = cli_opts.dataset_size();
  int carbon_num_runs = s13::gk_carbon_num_of_runs;
  if (num_runs != -1) {
    carbon_num_runs = 0.1971 * num_runs * 3;
    if (carbon_num_runs == 0) carbon_num_runs = 1;
  }

  auto carbon_ds = init_carbon_segmented_dataset(1, carbon_num_runs);
  auto carbon_yields_alg =
    s13::ana::walk_alg2<He4CarbonBgYields>(carbon_ds,
                                           cli_opts.use_mt(),
                                           cuts, k_bin_num);
  return carbon_yields_alg;
}

void he4_c_bg_subtraction() {
  auto cuts = s13::io::parse_cuts_config(std::fstream(cli_opts.cuts_config()));
  auto carbon_yields_alg = compute_carbon_yields(cuts);
  auto alg =
    s13::ana::walk_alg2<He4Yields>(init_dataset_from_cli(cli_opts),
                                   cli_opts.use_mt(),
                                   cuts, k_bin_num);

  auto script_fname =
    s13::ana::tstrfmt("%s_%s_he4_carbon_bg_subtraction",
                      cli_opts.output_file_basename().data(),
                      s13::ana::espri_selector_to_str(cuts.espri_selector).data());
  s13::ana::RootScript script(script_fname);
  gStyle->SetOptStat(1111111);
  gStyle->SetOptFit();
  draw_hists(script, alg, carbon_yields_alg);
}


int main(int argc, char** argv) {
  s13::misc::MessageLogger log("main()");
  // s13::ana::SetRootStyle(s13::ana::RootStyle::publication);

  TThread::Initialize();
  cli_opts.require_cuts_conf_opt(true);
  try {
    cli_opts.parse_options(argc, argv);
  } catch (std::exception& ex) {
    log.error("Invalid argument provided: %s", ex.what());
    return -1;
  }

  he4_c_bg_subtraction();
}
