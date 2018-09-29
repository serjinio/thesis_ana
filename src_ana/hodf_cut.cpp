
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


class HodfCutAlg : public s13::ana::SimpleTTreeAlgorithmBase {

public:

  HodfCutAlg(const ElaCuts& cuts) :
    logger_{"HodfCutAlg"}, cuts_{cuts} {
      init_hists();
    }

  virtual void process_event(Evt& evt) {
    if (cuts_.IsSignalEvent(evt)) {
      compute_signal_bg(evt);
      show_hodf_de_spectra(evt);
    }
  }

  virtual void finalize() {
    hist_in_out_sum = static_cast<TH2*>(hist_in_cut->Clone());
    hist_in_out_sum->Add(hist_out_cut);
    hist_in_out_sum->SetTitle("In + out HODF cut");
  }

  ElaCuts cuts() {
    return cuts_;
  }

  TH1* hist_in_cut;
  TH1* hist_out_cut;
  TH1* hist_in_out_sum;
  TH1* hist_total;
  TH2* hist_hodf_de;

private:

  void init_hists() {
    hist_in_cut =
      s13::ana::make_th1(100, -8, 8,
                         "In HODF cut",
                         "#theta_{exp} - #theta_{theor}",
                         "Counts");
    hist_out_cut =
      s13::ana::make_th1(100, -8, 8,
                         "Out HODF cut",
                         "#theta_{exp} - #theta_{theor}",
                         "Counts");
    hist_total =
      s13::ana::make_th1(100, -8, 8,
                         "Total",
                         "#theta_{exp} - #theta_{theor}",
                         "Counts");
    hist_hodf_de =
      s13::ana::make_th2(24, 0, 24, 600, 0, 40,
                         "HODF dE spectra", "Plastic ID", "#DeltaE [MeV]");
  }

  bool check_common_cuts(Evt& evt) {
    if (!is_well_defined(evt.he_theta) || !is_well_defined(evt.p_phi)) {
      return false;
    }

    return true;
  }
  void compute_signal_bg(Evt& evt) {
    if (!check_common_cuts(evt)) {
      return;
    }

    double he_theta_diff = evt.he_theta - evt.he_theta_theor;
    hist_total->Fill(he_theta_diff);
    if (cuts_.IsHodfHe6Cut(evt)) {
      hist_in_cut->Fill(he_theta_diff);
    } else {
      hist_out_cut->Fill(he_theta_diff);
    }
  }

  void show_hodf_de_spectra(Evt& evt) {
    int hodf_id = cuts_.get_hodf_hit_id(evt);
    double hodf_q = cuts_.get_hodf_q(evt);

    if (std::abs(evt.fdc0_xpos) > 20 && std::abs(evt.fdc0_xpos) < 50) {
      hist_hodf_de->Fill(hodf_id, hodf_q);
    }
  }


  s13::misc::MessageLogger logger_;

  ElaCuts cuts_;
};


void draw_hodf_cut(s13::ana::RootScript& script, std::shared_ptr<HodfCutAlg> alg) {
  script.NewPage(1).cd();
  s13::ana::draw_thstack({alg->hist_in_cut,
        alg->hist_out_cut, alg->hist_in_out_sum, alg->hist_total}, "HODF cut", "NOSTACK");
  s13::ana::draw_text(alg->cuts().GetTotalCutTitle());
  gPad->BuildLegend(0.65, 0.75, 0.9, 0.9);

  script.NewPage(1).cd();
  s13::ana::draw_thstack({alg->hist_in_cut,
        alg->hist_out_cut, alg->hist_in_out_sum}, "HODF cut", "NOSTACK");
  s13::ana::draw_text(alg->cuts().GetTotalCutTitle());
  gPad->BuildLegend(0.65, 0.75, 0.9, 0.9);

  script.NewPage(1).cd();
  s13::ana::draw_thstack({alg->hist_in_cut}, "HODF cut", "NOSTACK");
  s13::ana::draw_text(alg->cuts().GetTotalCutTitle());
  gPad->BuildLegend(0.65, 0.75, 0.9, 0.9);

  script.NewPage(1, 2).cd(1);
  alg->hist_in_cut->Draw();
  s13::ana::draw_text(alg->cuts().GetTotalCutTitle());
  script.cd(2);
  alg->hist_out_cut->Draw();

  script.NewPage(1).cd();
  gPad->SetLogz();
  alg->hist_hodf_de->SetNdivisions(5, "xyz");
  alg->hist_hodf_de->Draw("col");

}

void hodf_cut() {
  auto cuts = s13::io::parse_cuts_config(std::fstream(cli_opts.cuts_config()));
  cuts.apply_reaction_trigger_cut = false;
  cuts.apply_beam_trigger_cut = false;
  auto alg =
    s13::ana::walk_alg<HodfCutAlg>(cuts,
                                   init_dataset_from_cli(cli_opts),
                                   cli_opts.use_mt());

  s13::ana::RootScript script("hodf_cut");
  gStyle->SetOptStat(1111111);
  // gStyle->SetOptFit();
  draw_hodf_cut(script, alg);
}


int main(int argc, char** argv) {
  s13::ana::SetRootStyle(s13::ana::RootStyle::publication);
  s13::misc::MessageLogger log("main()");

  TThread::Initialize();
  cli_opts.require_cuts_conf_opt(true);
  try {
    cli_opts.parse_options(argc, argv);
  } catch (std::exception& ex) {
    log.error("Invalid argument provided: %s", ex.what());
    return -1;
  }

  hodf_cut();
  std::cout << "exiting main\n" << std::flush;
}
