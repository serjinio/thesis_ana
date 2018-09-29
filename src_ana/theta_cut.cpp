
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


class ThetaCutAlg : public s13::ana::SimpleTTreeAlgorithmBase {

public:

  static constexpr double k_min_theta = 55;
  static constexpr double k_max_theta = 70;

  ThetaCutAlg(const ElaCuts& cuts) :
    logger_{"ThetaCutAlg"}, cuts_{cuts}, bin_num_{4} {
      init_hists();
    }

  virtual void process_event(Evt& evt) {
    if (cuts_.IsSignalEvent(evt)) {
      compute_signal_bg(evt);
    }
  }

  void fit_kin_corr_hists(double bw_theta_width = 2.4) {
    constexpr double corr_width_factor = 0.07894736842105263;
    const double max_center_theta = bin_center_theta(hists_in_cut.size() - 1);
    int hist_idx = 0;
    for (auto el : hists_in_cut) {
      double hist_min_theta = -bw_theta_width/2
      + std::abs(max_center_theta - bin_center_theta(hist_idx)) * corr_width_factor;
      double hist_max_theta = std::abs(hist_min_theta);
      logger_.debug("Fitting kin. corr hist with central p_theta: %.1f in "
                    "He_theta region: +/-%.1f",
                    bin_center_theta(hist_idx), hist_max_theta);

      TF1 *fit_fn = new TF1("fit_fn","gaus(0)", hist_min_theta, hist_max_theta);
      el->Fit(fit_fn, "R");
      hist_idx++;
    }
  }

  ElaCuts cuts() {
    return cuts_;
  }

  std::vector<TH1*> hists_in_cut;
  std::vector<TH1*> hists_out_cut;
  std::vector<TH1*> hists_total;
  TH2* hist_kin_corr;

  void bin_num(int value) {
    bin_num_ = value;
  }

  size_t bin_num() {
    return bin_num_;
  }

  double bin_width() {
    return (k_max_theta - k_min_theta) / bin_num_;
  }

  double bin_center_theta(size_t bin_idx) {
    return k_min_theta + bin_width()/2 + bin_width()*bin_idx;
  }

  double bin_min_theta(size_t bin_idx) {
    return bin_center_theta(bin_idx) - bin_width()/2;
  }

  double bin_max_theta(size_t bin_idx) {
    return bin_center_theta(bin_idx) + bin_width()/2;
  }

private:

  std::vector<TH1*> init_1d_theta_hists(TString title) {
    std::vector<TH1*> res;
    for (size_t i = 0; i < bin_num_; ++i) {
      auto tit = s13::ana::tstrfmt("%s (#theta: %.1f)",
                                   title.Data(), bin_center_theta(i));
      auto hist = s13::ana::make_th1(100, -8, 8,
                                     tit,
                                     "#theta_{exp} - #theta_{theor}",
                                     "Counts");
      res.push_back(hist);
    }
    return res;
  }

  void init_hists() {
    hists_in_cut = init_1d_theta_hists("In cut");
    hists_out_cut = init_1d_theta_hists("Out cut");
    hists_total = init_1d_theta_hists("Total yield");
    hist_kin_corr =
      s13::ana::make_th2(100, 0, 10, 100, 55, 75,
                         "Kin. corr.", "#theta_{{}^{6}He} [deg.]", "#theta_{p} [deg.]");
  }


  int get_theta_bin_idx(double proton_theta) {
    double offset = proton_theta - k_min_theta;
    if (offset < 0) {
      return -1;
    }
    size_t bin_idx = offset / bin_width();
    if (bin_idx >= bin_num_) {
      return -1;
    }

    return bin_idx;
  }

  void compute_signal_bg(Evt& evt) {
    Double_t he_theta = cuts_.dsdc_selector == s13::ana::Dsdcs::s1dc
      ? evt.s1dc_theta : evt.fdc0_theta;

    hist_kin_corr->Fill(he_theta, evt.p_theta_eff);

    int theta_bin_idx = get_theta_bin_idx(evt.p_theta_eff);
    if (theta_bin_idx == -1) {
      return;
    }
    hists_in_cut.at(theta_bin_idx)->Fill(he_theta - evt.he_theta_theor);
  }


  s13::misc::MessageLogger logger_;

  ElaCuts cuts_;

  size_t bin_num_;
};


void draw_theta_cut(s13::ana::RootScript& script, std::shared_ptr<ThetaCutAlg> alg) {
  script.NewPage(2, 2).cd(s13::ana::PadSeq::row);
  for (size_t i = 0; i < alg->bin_num(); i++) {
    script.cd();
    alg->hists_in_cut.at(i)->Draw();
    if (i == 0) {
      s13::ana::draw_text(alg->cuts().GetTotalCutTitle(), 0.1, 0.8, 0.6, 0.9);
    }
  }

  script.NewPage(1).cd(1);
  alg->hist_kin_corr->SetNdivisions(6, "xyz");
  // gPad->SetLogz();
  // alg->hist_kin_corr->SetLineColor(kBlue);
  alg->hist_kin_corr->Draw("COL");
  auto kin_graph = s13::ana::BuildHe6ProtonKinematicsGraph();
  kin_graph->Draw("SAME");
}

void theta_cut() {
  auto cuts = s13::io::parse_cuts_config(std::fstream(cli_opts.cuts_config()));
  auto alg =
    s13::ana::walk_alg<ThetaCutAlg>(cuts,
                                   init_dataset_from_cli(cli_opts),
                                   cli_opts.use_mt());
  if (cuts.dsdc_selector == s13::ana::Dsdcs::s1dc) {
    alg->fit_kin_corr_hists(3.0);
  } else {
    alg->fit_kin_corr_hists(3.0);
  }

  s13::ana::RootScript script("theta_cut");
  gStyle->SetOptStat(1111111);
  gStyle->SetOptFit();
  draw_theta_cut(script, alg);
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

  theta_cut();
  std::cout << "exiting main\n" << std::flush;
}
