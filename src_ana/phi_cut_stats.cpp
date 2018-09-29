
#include "TROOT.h"

#include "cli.hpp"
#include "init_ds.hpp"
#include "consts.hpp"
#include "csalg.hpp"
#include "drawing.hpp"
#include "scattyield.hpp"
#include "treewalk.hpp"
#include "cuts_conf.hpp"
#include "rootscript.hpp"
#include "msgoutput.hpp"
#include "common.hpp"


static s13::misc::CliOptions cli_opts;


using ElaCuts = s13::ana::ElasticScatteringCuts;
using Evt = s13::ana::ScatteringEvent;


class PhiCutStatsAlg : public s13::ana::SimpleTTreeAlgorithmBase {

public:

  static constexpr double k_min_theta = 55;
  static constexpr double k_max_theta = 70;

  PhiCutStatsAlg(const ElaCuts& cuts) :
    logger_{"PhiCutStatsAlg"}, cuts_{cuts} {
      init_hists();
    }

  virtual void process_event(Evt& evt) {
    if (cuts_.IsSignalEvent(evt)) {
      show_phi_dist(evt);
    }

    if (cuts_.IsBackgroundEvent(evt)) {
      show_phi_bg_dist(evt);
    }
  }

  virtual void finalize() {
  }

  ElaCuts cuts() {
    return cuts_;
  }

  TH1* hist_phi_cut_total_range;
  TH1* hist_phi_cut_zero_region;
  TH1* hist_phi_cut_esl_region;
  TH1* hist_phi_cut_esr_region;
  TH1* hist_phi_cut_abs_range;
  TH1* hist_phi_cut_bg_total_range;

  TH2* hist_kin_corr_phi_zero_region;
  TH2* hist_kin_corr_phi_espri_region;

private:

  void init_hists() {
    auto th1 = &s13::ana::make_th1;

    hist_phi_cut_total_range =
      th1(720, -360, 360, "Phi dist. total", "#phi_p - #phi_{He}", "Counts");
    hist_phi_cut_bg_total_range =
      th1(720, -360, 360, "Phi dist. BG total", "#phi_p - #phi_{He}", "Counts");
    hist_phi_cut_zero_region =
      th1(200, -100, 100, "Phi dist. around zero region", "#phi_p - #phi_{He}", "Counts");
    hist_phi_cut_esl_region =
      th1(160, -260, -100, "Phi dist. ESL region", "#phi_p - #phi_{He}", "Counts");
    hist_phi_cut_esr_region =
      th1(160, 100, 260, "Phi dist. ESR region", "#phi_p - #phi_{He}", "Counts");
    hist_phi_cut_abs_range =
      th1(200, 0, 360, "Phi dist. total (absolute diff.)", "|#phi_p - #phi_{He}|", "Counts");

    hist_kin_corr_phi_zero_region =
      s13::ana::make_th2(100, 0, 10, 100, 50, 75,
                         "Kin. corr. (|#phi_{p} - #phi_{He}| < 100)",
                         "#theta_{{}^{6}He} [deg.]", "#theta_{p} [deg.]");
    hist_kin_corr_phi_espri_region =
      s13::ana::make_th2(100, 0, 10, 100, 50, 75,
                         "Kin. corr. (|#phi_{p} - #phi_{He}| > 100)",
                         "#theta_{{}^{6}He} [deg.]", "#theta_{p} [deg.]");
  }

  bool check_common_cuts(Evt& evt) {
    if (!is_well_defined(evt.p_phi) || !is_well_defined(evt.he_phi)) {
      return false;
    }
    if (std::abs(evt.p_e - evt.p_e_sim) > 0.1 * evt.p_e_sim) {
      return false;
    }

    // // for eff. evaluation cut on s1dc region
    if (!(std::abs(evt.s1dc_xpos) > 32 * 0.66
          && std::abs(evt.s1dc_xpos) < 160 * 0.66
          && std::abs(evt.s1dc_ypos) < 60)) {
      return false;
    }

    // // the same cut but using angular data of S1DC
    // if (evt.s1dc_theta < 4) {
    //   return false;
    // }
    // if (std::abs(evt.s1dc_phi) < 70 && std::abs(evt.s1dc_phi) > 110) {
    //   return false;
    // }

    return true;
  }

  void show_phi_dist(Evt& evt) {
    if (!check_common_cuts(evt)) {
      return;
    }
    total_events_ += 1;

    // logger_.debug("p_phi: %.1f; he_phi: %.1f", evt.p_phi, evt.he_phi);
    double copl_angle = evt.p_phi - evt.he_phi;
    hist_phi_cut_total_range->Fill(copl_angle);
    hist_phi_cut_abs_range->Fill(std::abs(copl_angle));
    hist_phi_cut_zero_region->Fill(copl_angle);
    if (copl_angle < 0) {
      hist_phi_cut_esl_region->Fill(copl_angle);
    } else {
      hist_phi_cut_esr_region->Fill(copl_angle);
    }

    if (std::abs(copl_angle) < 100) {
      hist_kin_corr_phi_zero_region->Fill(evt.he_theta, evt.p_theta_eff);
    }
    if (std::abs(copl_angle) > 100) {
      hist_kin_corr_phi_espri_region->Fill(evt.he_theta, evt.p_theta_eff);
    }
  }

  void show_phi_bg_dist(Evt& evt) {
    if (!check_common_cuts(evt)) {
      return;
    }

    double copl_angle = evt.p_phi - evt.he_phi;
    hist_phi_cut_bg_total_range->Fill(copl_angle);
  }

  s13::misc::MessageLogger logger_;

  ElaCuts cuts_;

  int total_events_ = 0;
};


void draw_phi_cut_stats(s13::ana::RootScript& script, std::shared_ptr<PhiCutStatsAlg> alg) {
  script.NewPage(1).cd();
  alg->hist_phi_cut_total_range->Draw();
  s13::ana::draw_text(alg->cuts().GetTotalCutTitle() + " react_beam@S1DC");
  script.NewPage(1).cd();
  alg->hist_phi_cut_bg_total_range->Draw();
  s13::ana::draw_text(alg->cuts().GetTotalCutTitle() + " react_beam@S1DC");

  script.NewPage(1).cd();
  alg->hist_phi_cut_abs_range->Draw();
  script.NewPage(1).cd();
  alg->hist_phi_cut_zero_region->Draw();
  script.NewPage(2, 1).cd();
  alg->hist_phi_cut_esl_region->Draw();
  script.cd();
  alg->hist_phi_cut_esr_region->Draw();

  script.NewPage().cd();
  alg->hist_kin_corr_phi_zero_region->Draw("col");
  s13::ana::draw_text(alg->cuts().GetTotalCutTitle() + " |#phi_{p} - #phi_{He}| < 100");
  script.NewPage().cd();
  alg->hist_kin_corr_phi_espri_region->Draw("col");
  s13::ana::draw_text(alg->cuts().GetTotalCutTitle() + " |#phi_{p} - #phi_{He}| > 100");
}

void phi_cut_stats() {
  auto cuts = s13::io::parse_cuts_config(std::fstream(cli_opts.cuts_config()));
  auto alg =
    s13::ana::walk_alg<PhiCutStatsAlg>(cuts,
                                  init_dataset_from_cli(cli_opts),
                                  cli_opts.use_mt());

  s13::ana::RootScript script("phi_cut_stats");
  gStyle->SetOptStat(1111111);
  // gStyle->SetOptFit();
  draw_phi_cut_stats(script, alg);
}


int main(int argc, char** argv) {
  s13::misc::MessageLogger log("main()");
  gROOT->SetStyle("Modern");
  TThread::Initialize();
  cli_opts.require_cuts_conf_opt(true);
  try {
    cli_opts.parse_options(argc, argv);
  } catch (std::exception& ex) {
    log.error("Invalid argument provided: %s", ex.what());
    return -1;
  }

  phi_cut_stats();
  std::cout << "exiting main\n" << std::flush;
}
