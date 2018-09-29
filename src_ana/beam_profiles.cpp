

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


class BeamProfileAlg : public s13::ana::SimpleTTreeAlgorithmBase {

public:

  BeamProfileAlg(const ElaCuts& cuts) :
    logger_{"BeamProfileAlg"}, cuts_{cuts}
  {
    init_hists();
  }

  virtual void process_event(Evt& evt) {
    if (cuts_.IsSignalEvent(evt)) {
      compute_beam_profiles(evt);
    }
  }

  virtual void finalize() {
  }

  ElaCuts cuts() {
    return cuts_;
  }

  // upstream info
  TH1* us_beam_x;
  TH1* us_beam_y;
  TH1* us_beam_aang;
  TH1* us_beam_bang;
  TH2* us_beam_x_vs_aang;
  TH2* us_beam_y_vs_bang;

  // downstream info
  TH1* ds_beam_x;
  TH1* ds_beam_y;
  TH1* ds_beam_aang;
  TH1* ds_beam_bang;
  TH2* ds_beam_x_vs_aang;
  TH2* ds_beam_y_vs_bang;

private:

  TH1* beam_pos_hist(TString hist_type) {
    TString title = s13::ana::tstrfmt("Beam %s position", hist_type.Data());
    return s13::ana::make_th1(400, -100, 100, title, hist_type);
  }

  TH1* beam_angle_hist(TString hist_type, Double_t anglim = 3) {
    TString title = s13::ana::tstrfmt("Beam %s angle", hist_type.Data());
    return s13::ana::make_th1(400, -anglim, anglim, title, hist_type);
  }

  TH2* beam_pos_angle_hist(TString hist_type, Double_t anglim = 3) {
    TString title = s13::ana::tstrfmt("Beam %s", hist_type.Data());
    return s13::ana::make_th2(400, -100, 100, 400, -anglim, anglim,
                              title, hist_type);
  }

  void init_hists() {
    us_beam_x = beam_pos_hist("upstream X");
    us_beam_y = beam_pos_hist("upstream Y");
    us_beam_aang = beam_angle_hist("upstream aang", 2);
    us_beam_bang = beam_angle_hist("upstream bang", 2);
    us_beam_x_vs_aang = beam_pos_angle_hist("upstream: X vs. aang", 2);
    us_beam_y_vs_bang = beam_pos_angle_hist("upstream: Y vs. bang", 2);

    ds_beam_x = beam_pos_hist("downstream X");
    ds_beam_y = beam_pos_hist("downstream Y");
    ds_beam_aang = beam_angle_hist("downstream aang", 12);
    ds_beam_bang = beam_angle_hist("downstream bang", 12);
    ds_beam_x_vs_aang = beam_pos_angle_hist("downstream: X vs. aang", 12);
    ds_beam_y_vs_bang = beam_pos_angle_hist("downstream: Y vs. bang", 12);
  }

  void compute_beam_profiles(Evt& evt) {
    // if (evt.tgt_up_aang > -9999) {
    //   logger_.debug("us aang: %.2f", evt.tgt_up_aang);
    //   logger_.debug("us aang [deg]: %.2f", evt.tgt_up_aang * 1/s13::gk_d2r);
    // }

    us_beam_x->Fill(evt.vert_xpos);
    us_beam_y->Fill(evt.vert_ypos);
    us_beam_aang->Fill(evt.tgt_up_aang * 1/s13::gk_d2r);
    us_beam_bang->Fill(evt.tgt_up_bang * 1/s13::gk_d2r);
    us_beam_x_vs_aang->Fill(evt.vert_xpos, evt.tgt_up_aang * 1/s13::gk_d2r);
    us_beam_y_vs_bang->Fill(evt.vert_ypos, evt.tgt_up_bang * 1/s13::gk_d2r);

    ds_beam_x->Fill(evt.tgt_dwn_vert_xpos);
    ds_beam_y->Fill(evt.tgt_dwn_vert_ypos);
    ds_beam_aang->Fill(evt.tgt_dwn_aang * 1/s13::gk_d2r);
    ds_beam_bang->Fill(evt.tgt_dwn_bang * 1/s13::gk_d2r);
    ds_beam_x_vs_aang->Fill(evt.tgt_dwn_vert_xpos, evt.tgt_dwn_aang * 1/s13::gk_d2r);
    ds_beam_y_vs_bang->Fill(evt.tgt_dwn_vert_ypos, evt.tgt_dwn_bang * 1/s13::gk_d2r);
  }

  s13::misc::MessageLogger logger_;

  ElaCuts cuts_;

};


void draw_beam_profiles(s13::ana::RootScript& script, std::shared_ptr<BeamProfileAlg> alg) {
  script.NewPage(2,2).cd(); alg->us_beam_x->Draw();
  script.cd(); alg->us_beam_y->Draw();
  script.cd(); alg->us_beam_aang->Draw();
  script.cd(); alg->us_beam_bang->Draw();
  script.NewPage(2,1).cd(); alg->us_beam_x_vs_aang->Draw("colz");
  gPad->SetLogz();
  script.cd(); alg->us_beam_y_vs_bang->Draw("colz");
  gPad->SetLogz();

  script.NewPage(2,2).cd(); alg->ds_beam_x->Draw();
  script.cd(); alg->ds_beam_y->Draw();
  script.cd(); alg->ds_beam_aang->Draw();
  script.cd(); alg->ds_beam_bang->Draw();
  script.NewPage(2,1).cd(); alg->ds_beam_x_vs_aang->Draw("colz");
  gPad->SetLogz();
  script.cd(); alg->ds_beam_y_vs_bang->Draw("colz");
  gPad->SetLogz();
}

void beam_profiles() {
  auto cuts = s13::io::parse_cuts_config(std::fstream(cli_opts.cuts_config()));
  auto alg =
    s13::ana::walk_alg<BeamProfileAlg>(cuts,
                                       init_dataset_from_cli(cli_opts),
                                       cli_opts.use_mt());

  s13::ana::RootScript script("beam_profiles");
  gStyle->SetOptStat(1111111);
  // gStyle->SetOptFit();
  draw_beam_profiles(script, alg);
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

  beam_profiles();
  std::cout << "exiting main\n" << std::flush;
}
