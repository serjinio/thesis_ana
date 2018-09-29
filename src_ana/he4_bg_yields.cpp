

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


class He4BgYields : public s13::ana::SimpleTTreeAlgorithmBase {

public:

  static constexpr double k_min_theta = s13::gk_cs_range_start;
  static constexpr double k_max_theta = s13::gk_cs_range_end;

  std::vector<TH1*> sn_in_hists;
  std::vector<TH1*> sn_out_hists;
  std::vector<TH1*> sn_tot_hists;
  TH2* in_kin_corr_hist;
  TH2* out_kin_corr_hist;

private:

  s13::misc::MessageLogger logger_;

  ElaCuts cuts_;

  size_t bin_num_;

  s13::ana::TInterpPtr he4_kin_interp_;

public:

  He4BgYields(const ElaCuts& cuts, size_t bin_num = 15) :
    logger_{"He4BgYields"}, cuts_{cuts}, bin_num_{bin_num} {
    init_hists();
    init_kin_interp();
  }

  He4BgYields(const He4BgYields&) = delete;
  He4BgYields(He4BgYields&&) = delete;
  He4BgYields operator=(const He4BgYields&) = delete;
  He4BgYields operator=(He4BgYields&&) = delete;

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
      sn_in_hists.push_back(init_kin_corr_hist("IN " + ang_range));
      sn_out_hists.push_back(init_kin_corr_hist("OUT " + ang_range));
      sn_tot_hists.push_back(init_kin_corr_hist("TOT " + ang_range));
    }

    in_kin_corr_hist = init_2d_kin_corr_hist("IN kin. corr.");
    out_kin_corr_hist = init_2d_kin_corr_hist("OUT kin. corr.");
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
      sn_in_hists[evt_bin_idx]->Fill(he_theta_diff);
    } else if (cuts_.IsBackgroundEvent(evt)){
      out_kin_corr_hist->Fill(he_theta, p_theta);
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

std::pair<TH1*, TH1*>
cleanup_he4_bg_events(TH1* bg_hist, TH1* signal_hist,
                      double theta_width) {
  TH1* clean_bg = static_cast<TH1*>(bg_hist->Clone());
  auto tit = s13::ana::tstrfmt("%s (w/o IN fraction)", bg_hist->GetTitle());
  auto name = s13::ana::tstrfmt("%s_s", bg_hist->GetName());
  clean_bg->SetTitle(tit);
  clean_bg->SetName(name);

  TH1* scaled_signal = static_cast<TH1*>(signal_hist->Clone());
  tit = s13::ana::tstrfmt("%s (scaled)", signal_hist->GetTitle());
  name = s13::ana::tstrfmt("%s_s", signal_hist->GetName());
  scaled_signal->SetTitle(tit);
  scaled_signal->SetName(name);

  double theta_min = -theta_width/2;
  double theta_max = theta_width/2;

  double bg_int =
    bg_hist->Integral(bg_hist->GetXaxis()->FindBin(theta_min),
                      bg_hist->GetXaxis()->FindBin(theta_max));
  double signal_int =
    signal_hist->Integral(signal_hist->GetXaxis()->FindBin(theta_min),
                          signal_hist->GetXaxis()->FindBin(theta_max));
  double scaling_f = bg_int / signal_int;
  scaled_signal->Scale(scaling_f);
  clean_bg->Add(scaled_signal, -1);

  return std::pair<TH1*, TH1*>(clean_bg, scaled_signal);
}

void
draw_bg_hists(s13::ana::RootScript& script,
              std::shared_ptr<He4BgYields> alg) {
  auto kin_graph = s13::ana::BuildHe4ProtonKinematicsGraph();
  script.NewPage(1,1).cd();
  alg->in_kin_corr_hist->Draw("colz");
  kin_graph->Draw("SAME");
  script.NewPage(1,1).cd();
  alg->out_kin_corr_hist->Draw("colz");
  kin_graph->Draw("SAME");

  SnsVector sns;
  for (size_t i = 0; i < alg->sn_tot_hists.size(); ++i) {
    if (i % 3 == 0) {
      script.NewPage(3, 2).cd(s13::ana::PadSeq::column);
    }

    // auto scaled_out_cut_hist =
    //   scale_out_phi_events(alg->sn_out_hists[i],
    //                        alg->bin_center_theta(i),
    //                        alg->cuts());
    // auto bg_cleanup_hists =
    //   cleanup_he4_bg_events(alg->sn_in_hists[i],
    //                         alg->sn_out_hists[i],
    //                         alg->cuts().theta_width/3);
    // auto clean_bg_hist = std::get<0>(bg_cleanup_hists);
    // auto scaled_signal_hist = std::get<1>(bg_cleanup_hists);

    script.cd();
    auto s1 = s13::ana::draw_thstack({
        alg->sn_in_hists[i], alg->sn_out_hists[i]},
      "IN #phi, OUT #phi", "NOSTACK");
    s1->SetMaximum(alg->sn_out_hists[i]->GetMaximum() * 1.5);
    gPad->BuildLegend(0.15, 0.65, 0.70, 0.85);

    // script.cd();
    // auto s2 = s13::ana::draw_thstack({
    //     alg->sn_in_hists[i], scaled_out_cut_hist},
    //   "IN #phi, OUT #phi (scaled)", "NOSTACK");
    // s2->SetMaximum(alg->sn_in_hists[i]->GetMaximum() * 1.5);
    // gPad->BuildLegend(0.15, 0.65, 0.70, 0.85);

    script.cd();
    alg->sn_in_hists[i]->Draw();
    draw_cut_marker_lines(alg->cuts().theta_width, alg->sn_in_hists[i]);
    double sn =
      s13::ana::draw_kin_corr_hist_sn(alg->bin_center_theta(i),
                                      alg->sn_in_hists[i],
                                      alg->cuts().dsdc_selector,
                                      alg->cuts().theta_width,
                                      s13::ana::UseAbsVal(false),
                                      s13::ana::Degrees(-11));
    sns.push_back(std::pair<double, double>(alg->bin_center_theta(i), sn));
  }

  auto espri_str =
    s13::ana::espri_selector_to_str(alg->cuts().espri_selector).data();
  auto fname =
    s13::ana::tstrfmt("out/%s_%s_he4_bg_yields_sns.csv",
                      cli_opts.output_file_basename().data(),
                      espri_str);
  write_sns(fname.Data(), sns);
}

void he4_bg_yields() {
  auto cuts = s13::io::parse_cuts_config(std::fstream(cli_opts.cuts_config()));
  auto alg =
    s13::ana::walk_alg2<He4BgYields>(init_dataset_from_cli(cli_opts),
                                     cli_opts.use_mt(),
                                     cuts);

  auto script_fname =
    s13::ana::tstrfmt("%s_%s_he4_bg_yields",
                      cli_opts.output_file_basename().data(),
                      s13::ana::espri_selector_to_str(cuts.espri_selector).data());
  s13::ana::RootScript script(script_fname);
  gStyle->SetOptStat(1111111);
  gStyle->SetOptFit();
  draw_bg_hists(script, alg);
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

  he4_bg_yields();
}
