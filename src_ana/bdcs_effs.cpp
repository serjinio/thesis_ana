
#include <fstream>

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


class BdcsEffsAlg : public s13::ana::SimpleTTreeAlgorithmBase {

public:

  BdcsEffsAlg(const ElaCuts& cuts) :
    logger_{"BdcsEffsAlg"}, cuts_{cuts} {
      init();
    }

  virtual void process_event(Evt& evt) {
    if (cuts_.IsSignalEvent(evt)) {
      compute_bdcs_effs(evt);
    }
  }

  virtual void finalize() {
    double bdc1_eff = bdcs_eff_hits_[0] / (double)bdcs_ref_hits_[1] * 100;
    double bdc2_eff = bdcs_eff_hits_[1] / (double)bdcs_ref_hits_[0] * 100;
    logger_.info("Total number of events: %d", total_events_);
    logger_.info("BDC1 in pos. cut hits: %d", bdcs_ref_hits_[0]);
    logger_.info("BDC2 in pos. cut hits: %d", bdcs_ref_hits_[1]);
    logger_.info("BDC2*BDC1 hits: %d", bdcs_eff_hits_[0]);
    logger_.info("BDC1*BDC2 hits: %d", bdcs_eff_hits_[1]);
    logger_.info("BDC1 efficiency: %.1f%%", bdc1_eff);
    logger_.info("BDC2 efficiency: %.1f%%", bdc2_eff);
    double bdc1_eff_to_sbt = sbt_bdc1_hits_ / (double)sbt_hits_ * 100;

    logger_.info("\nBDC1 efficiency relative to SBT:");
    logger_.info("\tSBT hits: %d", sbt_hits_);
    logger_.info("\tSBT*BDC1 hits: %d", sbt_bdc1_hits_);
    logger_.info("\tBDC1 efficiency: %.1f%%", bdc1_eff_to_sbt);
  }

  ElaCuts cuts() {
    return cuts_;
  }

private:

  void init() {
    // auto th1 = &s13::ana::make_th1;
    bdcs_ref_hits_ = { { 0, 0 } };
    bdcs_eff_hits_ = { { 0, 0 } };
  }

  bool is_bdc1_pos_cut(Evt& evt) {
    double evt_r = std::sqrt(std::pow(evt.bdc1_xpos, 2) + std::pow(evt.bdc1_ypos, 2));
    if (evt_r < 1) {
      return true;
    }
    return false;
  }

  bool is_bdc2_pos_cut(Evt& evt) {
    double evt_r = std::sqrt(std::pow(evt.bdc2_xpos, 2) + std::pow(evt.bdc2_ypos, 2));
    if (evt_r < 1) {
      return true;
    }
    return false;
  }

  bool is_bdc1_hit(Evt& evt) {
    if (is_well_defined(evt.bdc1_xpos) && is_well_defined(evt.bdc1_ypos)) {
      return true;
    }
    return false;
  }

  bool is_bdc2_hit(Evt& evt) {
    if (is_well_defined(evt.bdc2_xpos) && is_well_defined(evt.bdc2_ypos)) {
      return true;
    }
    return false;
  }

  bool is_sbt_hit(Evt& evt) {
    return true;
  }

  void compute_bdcs_effs(Evt& evt) {
    total_events_ += 1;

    if (is_bdc1_pos_cut(evt)) {
      ++bdcs_ref_hits_[0];
      if (is_bdc2_hit(evt)) {
        ++bdcs_eff_hits_[1];
      }
    }

    if (is_bdc2_pos_cut(evt)) {
      ++bdcs_ref_hits_[1];
      if (is_bdc1_hit(evt)) {
        ++bdcs_eff_hits_[0];
      }
    }

    if (is_sbt_hit(evt)) {
      ++sbt_hits_;
      if (is_bdc1_hit(evt)) {
        ++sbt_bdc1_hits_;
      }
    }
  }

  s13::misc::MessageLogger logger_;

  ElaCuts cuts_;

  int total_events_ = 0;
  int sbt_hits_ = 0;
  int sbt_bdc1_hits_ = 0;
  std::array<int, 2> bdcs_ref_hits_;
  std::array<int, 2> bdcs_eff_hits_;
};

void draw_bdcs_effs(s13::ana::RootScript& script, std::shared_ptr<BdcsEffsAlg> alg) {
  script.NewPage(1).cd();
}

void bdcs_effs() {
  auto cuts = s13::io::parse_cuts_config(std::fstream(cli_opts.cuts_config()));
  auto alg =
    s13::ana::walk_alg<BdcsEffsAlg>(cuts,
                                    init_dataset_from_cli(cli_opts),
                                    cli_opts.use_mt());

  s13::ana::RootScript script("bdcs_effs");
  gStyle->SetOptStat(1111111);
  // gStyle->SetOptFit();
  draw_bdcs_effs(script, alg);
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

  bdcs_effs();
}
