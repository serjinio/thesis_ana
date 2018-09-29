

#include "TROOT.h"
#include "TPDF.h"

#include "cli.hpp"
#include "init_ds.hpp"
#include "consts.hpp"
#include "csalg.hpp"
#include "drawing.hpp"
#include "scattyield.hpp"
#include "protonede.hpp"
#include "treewalk.hpp"
#include "cuts_conf.hpp"
#include "rootscript.hpp"


static s13::misc::CliOptions cli_opts;


using ElaCuts = s13::ana::ElasticScatteringCuts;
using Evt = s13::ana::ScatteringEvent;


class DsdcsEffsAlg : public s13::ana::SimpleTTreeAlgorithmBase {

public:

  static constexpr double k_min_theta = 55;
  static constexpr double k_max_theta = 70;

  DsdcsEffsAlg(const ElaCuts& cuts) :
    logger_{"DsdcsEffsAlg"}, cuts_{cuts} {
      init_hists();
    }

  DsdcsEffsAlg(const DsdcsEffsAlg&) = delete;
  DsdcsEffsAlg(DsdcsEffsAlg&&) = delete;
  DsdcsEffsAlg operator=(const DsdcsEffsAlg&) = delete;
  DsdcsEffsAlg operator=(DsdcsEffsAlg&&) = delete;

  virtual void setup() {
    hist_s1dc_yields->Reset();
    hist_fdc0_yields->Reset();
  }

  virtual void process_event(Evt& evt) {
    compute_kin_corr(evt);
    if (cuts_.IsSignalEvent(evt)) {
      compute_s1dc_eff(evt);
    }
  }

  virtual void finalize() {
    hist_fdc0_eff = static_cast<TH1*>(hist_fdc0_yields->Clone());
    hist_fdc0_eff->Divide(hist_s1dc_yields);
    hist_fdc0_eff->GetYaxis()->SetTitle("FDC0 eff. (FDC0/S1DC) [x100%]");
    hist_fdc0_eff->SetTitle("FDC0 efficiency");

    hist_s1dc_eff = static_cast<TH1*>(hist_s1dc_yields->Clone());
    hist_s1dc_eff->Divide(hist_fdc0_yields);
    hist_s1dc_eff->GetYaxis()->SetTitle("S1DC eff. (S1DC/FDC0) [x100%]");
    hist_s1dc_eff->SetTitle("S1DC efficiency");

    logger_.info("Stats for run: \n");
    logger_.info("Good FDC0 events: %i", good_fdc0_events_);
    logger_.info("Good S1DC events: %i", good_s1dc_events_);
  }

  virtual void merge(const TTreeAlgorithmBase& other) {
    // auto other_alg = static_cast<const DsdcsEffsAlg*>(&other);
    throw std::logic_error("not implemented");
  }

  virtual DsdcsEffsAlg* clone() const {
    DsdcsEffsAlg* alg = new DsdcsEffsAlg(cuts_);
    alg->merge(*this);
    return alg;
  }

  ElaCuts cuts() {
    return cuts_;
  }

  TH1* hist_fdc0_eff;
  TH1* hist_s1dc_eff;
  TH1* hist_s1dc_yields;
  TH1* hist_fdc0_yields;
  TH2* hist_kin_corr;

private:

  void init_hists() {
    auto th1 = &s13::ana::make_th1;
    // auto th2 = &s13::ana::make_th2;

    hist_s1dc_yields = th1(27, -13, 13,
                           "S1DC yields", "aang [deg.]",
                           "Counts");
    hist_fdc0_yields = th1(27, -13, 13,
                           "FDC0 yields", "aang [deg.]",
                           "Counts");
    hist_kin_corr =
      s13::ana::make_th2(80, 0, 13, 80, 50, 75,
                         "Kin. corr.", "#theta_{6He} [deg.]", "#theta_{p} [deg.]");
  }

  void compute_kin_corr(Evt& evt) {
    auto kin_corr_cuts = cuts_;
    kin_corr_cuts.apply_theta_cut = false;
    if (kin_corr_cuts.IsSignalEvent(evt)) {
      hist_kin_corr->Fill(evt.fdc0_theta, evt.p_theta_eff);
    }
  }

  bool hodf_z2(Evt& evt) {
    for (int i = 0; i < 23; i++) {
      if (evt.hodf_q[i] > 7) return true;
    }
    return false;
  }

  void compute_s1dc_eff(Evt& evt) {
    // HODF cut on Z=2 particles
    if (!hodf_z2(evt)) {
      return;
    }
    // cut for "reference" DSDC - should have good data
    if (!is_well_defined(evt.fdc0_phi)) {
      // logger_.debug("bad event returning!");
      return;
    }

    if (is_well_defined(evt.s1dc_phi)) {
      good_s1dc_events_ += 1;
      hist_s1dc_yields->Fill(evt.s1dc_aang);
    }

    if (is_well_defined(evt.fdc0_phi)) {
      good_fdc0_events_ += 1;
      hist_fdc0_yields->Fill(evt.fdc0_aang);
    }
  }


  s13::misc::MessageLogger logger_;

  ElaCuts cuts_;

  size_t good_fdc0_events_ = 0;
  size_t good_s1dc_events_ = 0;

};


void draw_dsdcs_effs(s13::ana::RootScript& script, std::shared_ptr<DsdcsEffsAlg> alg) {
  script.NewPage(1).cd();
  alg->hist_kin_corr->Draw("COL");
  s13::ana::draw_text(alg->cuts().GetTotalCutTitle());

  script.NewPage(2,1);
  script.cd(); alg->hist_s1dc_yields->Draw(); gPad->SetLogy();
  script.cd(); alg->hist_fdc0_yields->Draw(); gPad->SetLogy();

  script.NewPage(1).cd();
  alg->hist_fdc0_eff->SetMaximum(2.0);
  alg->hist_fdc0_eff->Draw();

  script.NewPage(1).cd();
  alg->hist_s1dc_eff->SetMaximum(2.0);
  alg->hist_s1dc_eff->Draw();

}

void dsdcs_effs() {
  auto cuts = s13::io::parse_cuts_config(std::fstream(cli_opts.cuts_config()));
  auto alg =
    s13::ana::walk_alg<DsdcsEffsAlg>(cuts,
                                     init_dataset_from_cli(cli_opts),
                                     cli_opts.use_mt());

  s13::ana::RootScript script("dsdcs_effs");
  gStyle->SetOptStat(1111111);
  // gStyle->SetOptFit();
  draw_dsdcs_effs(script, alg);
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

  dsdcs_effs();
}
