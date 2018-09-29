
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


static s13::misc::CliOptions cli_opts;


using ElaCuts = s13::ana::ElasticScatteringCuts;
using Evt = s13::ana::ScatteringEvent;


class PhiCutAlg : public s13::ana::SimpleTTreeAlgorithmBase {

public:

  static constexpr double k_min_theta = 55;
  static constexpr double k_max_theta = 70;

  PhiCutAlg(const ElaCuts& cuts) :
    logger_{"PhiCutAlg"}, cuts_{cuts} {
      init_hists();
    }

  virtual void process_event(Evt& evt) {
    show_phi_in_out_events(evt);
    if (cuts_.IsSignalEvent(evt)) {
      show_phi_dist(evt);
    }
  }

  virtual void finalize() {
    logger_.info("In range number of events: %d", total_events_);
    logger_.info("Out of range p_phi events: %d", out_of_range_p_phi);
    logger_.info("Out of range he_phi events: %d \n", out_of_range_he_phi);

    logger_.info("Input events (for eff. computation): %d", input_evts_);
    logger_.info("\t Untracked events by RDC: %d", untracked_rdc_);
    logger_.info("\t Untracked events by S1DC: %d", untracked_s1dc_);
    logger_.info("\t RDC eff.: %.1f%%", ((double)input_evts_ - untracked_rdc_) / input_evts_ * 100);
  }

  ElaCuts cuts() {
    return cuts_;
  }

  TH1* hist_phi_out_cut_total_range;
  TH1* hist_phi_in_cut_total_range;
  TH1* hist_phi_cut_local_range;
  TH1* hist_kin_corr_in_cut;
  TH1* hist_kin_corr_out_cut;
  TH1* hist_kin_corr_total;
  TH1* hist_rdcx;
  TH1* hist_rdcx_nai_cut;

private:

  void init_hists() {
    auto th1 = &s13::ana::make_th1;

    hist_phi_out_cut_total_range =
      th1(720, -360, 360, "Phi dist. total", "#phi_p - #phi_{He}", "Counts");
    hist_phi_in_cut_total_range =
      th1(720, -360, 360, "Phi dist. in cut", "#phi_p - #phi_{He}", "Counts");
    hist_phi_cut_local_range =
      th1(200, 180 - 32, 180 + 32, "Phi dist.", "#phi_p - #phi_{He}", "Counts");
    hist_kin_corr_total =
      th1(100, -8, 8, "Total", "#theta_{exp} - #theta_{theor}", "Counts");
    hist_kin_corr_in_cut =
      th1(100, -8, 8, "In cut", "#theta_{exp} - #theta_{theor}", "Counts");
    hist_kin_corr_out_cut =
      th1(100, -8, 8, "Out cut", "#theta_{exp} - #theta_{theor}", "Counts");

    TString title = s13::ana::tstrfmt("%s RDC X", cuts_.espri_selector ==
                                      s13::ana::Espri::left ?
                                      "ESPRI left" : "ESPRI right");
    hist_rdcx =
      th1(200, -250, 250, title, "X [mm]", "Counts");
    hist_rdcx_nai_cut =
      th1(200, -250, 250, title, "X [mm]", "Counts");
  }

  bool is_well_defined(double value) {
    if (value == -99999) {
      return false;
    } else {
      return true;
    }
  }

  void show_phi_dist(Evt& evt) {
    total_events_ += 1;
    if (std::abs(evt.p_phi) > 360) {
      out_of_range_p_phi += 1;
    }
    if (std::abs(evt.he_phi) > 360) {
      out_of_range_he_phi += 1;
    }

    if (!is_well_defined(evt.s1dc_phi)) {
      untracked_s1dc_ += 1;
      return;
    }

    if (evt.IsLeftEspriEvent() && evt.esl_e_raw < 100) {
      return;
    }
    if (evt.IsRightEspriEvent() && evt.esr_e_raw < 100) {
      return;
    }
    if (!(std::abs(evt.s1dc_xpos) > 32 * 0.66
          && std::abs(evt.s1dc_xpos) < 160 * 0.66
          && std::abs(evt.s1dc_ypos) < 60)) {
      return;
    }
    input_evts_ += 1;

    if (!is_well_defined(evt.p_phi)) {
      untracked_rdc_ += 1;
    }
  }

  void show_phi_in_out_events(Evt& evt) {
    ElaCuts phi_cuts = cuts_;
    phi_cuts.apply_phi_cut = false;
    Double_t he_theta = phi_cuts.dsdc_selector == s13::ana::Dsdcs::fdc0
      ? evt.fdc0_theta : evt.s1dc_theta;
    Double_t he_phi = phi_cuts.dsdc_selector == s13::ana::Dsdcs::fdc0
      ? evt.fdc0_phi : evt.s1dc_phi;
    if (!phi_cuts.IsSignalEvent(evt)) {
      return;
    }
    if (!is_well_defined(he_phi)) {
      return;
    }
    // if (evt.IsLeftEspriEvent() && evt.esl_e_raw < 100) {
    //   return;
    // }
    // if (evt.IsRightEspriEvent() && evt.esr_e_raw < 100) {
    //   return;
    // }

    // if (!(std::abs(evt.s1dc_xpos) > 32 * 0.66
    //       && std::abs(evt.s1dc_xpos) < 160 * 0.66
    //       && std::abs(evt.s1dc_ypos) < 60)) {
    //   return;
    // }

    if (std::abs(he_theta - evt.he_theta_theor) > 1.6) {
      return;
    }

    hist_kin_corr_total->Fill(he_theta - evt.he_theta_theor);
    hist_phi_out_cut_total_range->Fill(evt.p_phi - he_phi);
    hist_phi_cut_local_range->Fill(std::abs(evt.p_phi - he_phi));
    if (phi_cuts.IsPhiCut(evt)) {
      hist_phi_in_cut_total_range->Fill(evt.p_phi - he_phi);
      hist_kin_corr_in_cut->Fill(he_theta - evt.he_theta_theor);
    } else {
      hist_kin_corr_out_cut->Fill(he_theta - evt.he_theta_theor);
    }
  }

  s13::misc::MessageLogger logger_;

  ElaCuts cuts_;

  int total_events_ = 0;
  int out_of_range_p_phi = 0;
  int out_of_range_he_phi = 0;
  int untracked_rdc_ = 0;
  int untracked_s1dc_ = 0;
  int input_evts_ = 0;
};


void draw_phi_cut(s13::ana::RootScript& script, std::shared_ptr<PhiCutAlg> alg) {
  script.NewPage(1, 2).cd(1);
  alg->hist_phi_out_cut_total_range->Draw("");
  s13::ana::draw_text(alg->cuts().GetTotalCutTitle());
  gPad->SetLogy(true);
  script.cd(2);
  alg->hist_phi_in_cut_total_range->Draw("");
  s13::ana::draw_text(alg->cuts().GetTotalCutTitle());
  gPad->SetLogy(true);

  script.NewPage(1).cd();
  alg->hist_phi_cut_local_range->Draw("");

  script.NewPage(1).cd();
  s13::ana::draw_thstack({alg->hist_kin_corr_total, alg->hist_kin_corr_in_cut,
        alg->hist_kin_corr_out_cut}, "In & out phi cut events", "NOSTACK HIST");
  s13::ana::draw_text(alg->cuts().GetTotalCutTitle());
  gPad->BuildLegend();

  script.NewPage(1).cd();
  alg->hist_kin_corr_total->Draw();
  s13::ana::draw_text(alg->cuts().GetTotalCutTitle());

  script.NewPage(1).cd();
  alg->hist_kin_corr_in_cut->Draw();
  script.NewPage(1).cd();
  alg->hist_kin_corr_out_cut->Draw();

  script.NewPage(1, 2).cd(1);
  alg->hist_rdcx->Draw();
  s13::ana::draw_text(alg->cuts().GetTotalCutTitle());
  script.cd(2);
  alg->hist_rdcx_nai_cut->Draw();
  s13::ana::draw_text(alg->cuts().GetTotalCutTitle() + " NaI 'fired'");

}

void phi_cut() {
  auto cuts = s13::io::parse_cuts_config(std::fstream(cli_opts.cuts_config()));
  auto alg =
    s13::ana::walk_alg<PhiCutAlg>(cuts,
                                  init_dataset_from_cli(cli_opts),
                                  cli_opts.use_mt());

  s13::ana::RootScript script("phi_cut");
  gStyle->SetOptStat(1111111);
  // gStyle->SetOptFit();
  draw_phi_cut(script, alg);
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

  phi_cut();
  std::cout << "exiting main\n" << std::flush;
}
