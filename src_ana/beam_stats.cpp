
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


class BeamStatsAlg : public s13::ana::SimpleTTreeAlgorithmBase {

public:

  BeamStatsAlg(const ElaCuts& cuts) :
    logger_{"BeamStatsAlg"}, cuts_{cuts}
  {
    init_hists();
  }

  virtual void process_event(Evt& evt) {
    if (cuts_.IsSignalEvent(evt)) {
      compute_beam_stats(evt);
    }
  }

  virtual void finalize() {
    for (size_t i = 0; i < target_radii.size(); i++) {
      num_beam_particles[i] *= s13::gk_beam_dsc_factor;
    }
    make_beam_stats_graph();
    write_stats_to_file("data/he6_beam_pol_up_vs_tgtR.csv");
  }

  ElaCuts cuts() {
    return cuts_;
  }

  std::vector<double> num_beam_particles;
  std::vector<double> target_radii;
  TGraph* beam_tgtr_dep;

private:

  void make_beam_stats_graph() {
    double tgt_radii[target_radii.size()];
    double beam_part[target_radii.size()];
    for (size_t i = 0; i < target_radii.size(); i++) {
      tgt_radii[i] = target_radii[i];
      beam_part[i] = num_beam_particles[i];
    }

    beam_tgtr_dep = new TGraph(target_radii.size(), tgt_radii, beam_part);
    beam_tgtr_dep->SetLineColor(kRed);
    beam_tgtr_dep->SetLineWidth(5);
    beam_tgtr_dep->SetMarkerColor(5);
    beam_tgtr_dep->SetMarkerStyle(22);
    beam_tgtr_dep->SetTitle("Number of beam particles vs. target radius");
    beam_tgtr_dep->GetXaxis()->SetTitle("Target radius [mm]");
    beam_tgtr_dep->GetYaxis()->SetTitle("Num. beam [#]");
  }

  void write_stats_to_file(std::string filename) {
    std::ofstream ofs(filename);
    ofs << "# tgtR,Nbeam" << std::endl;
    for (size_t i = 0; i < target_radii.size(); i++) {
      ofs << target_radii[i] << "," << num_beam_particles[i] << std::endl;
    }
  }

  void init_hists() {
    // auto th1 = &s13::ana::make_th1;

    for (int i = 1; i < 25; i += 1) {
      target_radii.push_back(i);
      num_beam_particles.push_back(0);
    }
  }

  bool check_common_cuts(Evt& evt) {
    if (!is_well_defined(evt.vert_xpos) || !is_well_defined(evt.vert_ypos)) {
      return false;
    }
    return true;
  }

  void compute_beam_stats(Evt& evt) {
    if (!check_common_cuts(evt)) {
      return;
    }
    total_events_ += 1;

    for (size_t idx = 0; idx < target_radii.size(); idx++) {
      cuts_.vertexXY_radius = target_radii[idx];
      if (cuts_.IsVertexXYCut(evt)) {
        num_beam_particles[idx] += 1;
      }
    }
  }

  s13::misc::MessageLogger logger_;

  ElaCuts cuts_;

  int total_events_ = 0;
};

void draw_beam_stats(s13::ana::RootScript& script, std::shared_ptr<BeamStatsAlg> alg) {
  script.NewPage(1).cd();
  alg->beam_tgtr_dep->Draw("ACP");
}

void beam_stats() {
  auto cuts = s13::io::parse_cuts_config(std::fstream(cli_opts.cuts_config()));
  auto alg =
    s13::ana::walk_alg<BeamStatsAlg>(cuts,
                                  init_dataset_from_cli(cli_opts),
                                  cli_opts.use_mt());

  s13::ana::RootScript script("beam_stats");
  gStyle->SetOptStat(1111111);
  // gStyle->SetOptFit();
  draw_beam_stats(script, alg);
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

  beam_stats();
  std::cout << "exiting main\n" << std::flush;
}
