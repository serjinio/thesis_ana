

#include "TROOT.h"
#include "TPDF.h"
#include "TPaveStats.h"

#include "common.hpp"
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
#include "vertzalg.hpp"


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


class ESpectra : public s13::ana::SimpleTTreeAlgorithmBase {

public:

  static constexpr double k_min_theta = 55;
  static constexpr double k_max_theta = 70;

  TH1* sn_hist;
  TH2* signal_theta_vs_e_hist;
  TH2* bg_theta_vs_e_hist;
  TH1* signal_1d_e_hist;
  TH1* bg_1d_e_hist;

  std::vector<double> e_cut_widths = {0.5, 1, 3};
  std::vector<TH1*> e_cut_hists;
  std::vector<TH1*> e_cut_sn_in_hists;
  std::vector<TH1*> e_cut_sn_out_hists;
  std::vector<TH1*> e_cut_sn_tot_hists;

private:

  s13::misc::MessageLogger logger_;

  ElaCuts cuts_;

  size_t bin_num_;

public:

  ESpectra(const ElaCuts& cuts) :
    logger_{"ESpectra"}, cuts_{cuts}, bin_num_{8} {
      init_hists();
    }

  ESpectra(const ESpectra&) = delete;
  ESpectra(ESpectra&&) = delete;
  ESpectra operator=(const ESpectra&) = delete;
  ESpectra operator=(ESpectra&&) = delete;

  virtual void setup() {
  }

  virtual void process_event(Evt& evt) {
    if (cuts_.IsSignalEvent(evt)) {
      fill_e_spectra(evt);
      fill_e_cut_hists(evt);
    }
  }

  virtual void finalize() {
    // compute_bg_scaling_factors();
  }

  ElaCuts cuts() {
    return cuts_;
  }

  size_t bin_num() {
    return bin_num_;
  }

private:

  TH1* init_kin_corr_hist(TString title) {
    auto th1 = &s13::ana::make_th1;
    auto hist = th1(100, -8, 8, title, "#theta_{exp} - #theta_{theor}", "Counts");
    return hist;
  }

  TH2* init_2d_kin_corr_hist(TString title) {
    auto th2 = &s13::ana::make_th2;
    auto hist = th2(200, 0, 12, 200, 50, 75, title, "#theta_{He}", "#theta_{p}");
    return hist;
  }

  TH2* init_p_theta_vs_e_hist(TString title) {
    auto th2 = &s13::ana::make_th2;
    auto hist = th2(200, 52, 73, 200, 0, 100, title, "#theta_{p}", "E [MeV]");
    return hist;
  }

  TH1* init_1d_e_hist(TString title) {
    auto th1 = &s13::ana::make_th1;
    auto hist = th1(200, -1, 1, title, "#Delta E (E_{kin} - E_{exp}) / E_{kin}", "Counts");
    return hist;
  }

  void init_hists() {
    sn_hist = init_kin_corr_hist("1D kin. corr.");
    signal_theta_vs_e_hist = init_p_theta_vs_e_hist("Angle - energy corr. (IN CUT)");
    bg_theta_vs_e_hist = init_p_theta_vs_e_hist("Angle - energy corr. (OUT CUT)");
    signal_1d_e_hist = init_1d_e_hist("1-d E dist. (IN CUT)");
    bg_1d_e_hist = init_1d_e_hist("1-d E dist. (OUT CUT)");

    for (size_t i = 0; i < e_cut_widths.size(); i++) {
      auto title = s13::ana::tstrfmt("E cut: +/-%.2f", e_cut_widths[i]/2);
      e_cut_hists.push_back(init_1d_e_hist(title));
      e_cut_sn_in_hists.push_back(init_kin_corr_hist("IN"));
      e_cut_sn_out_hists.push_back(init_kin_corr_hist("OUT"));
      e_cut_sn_tot_hists.push_back(init_kin_corr_hist("TOT"));
    }
  }

  virtual void fill_e_spectra(Evt& evt) {
    if (!is_well_defined(evt.p_theta)
        || !is_well_defined(evt.he_theta)) {
      return;
    }
    double xpos = cuts_.IsLeftEspriEvent(evt) ? evt.esl_xpos : evt.esr_xpos;
    if (std::abs(xpos) > 9000) {
      return;
    }
    if (( cuts_.IsLeftEspriEvent(evt) && evt.esl_de_raw < 240 ) ||
        ( cuts_.IsRightEspriEvent(evt) && evt.esr_de_raw < 240 )) {
      return;
    }

    double he_theta_diff = evt.he_theta - evt.he_theta_theor;
    double p_theta = evt.p_theta_eff;
    double p_e = evt.p_e;
    double p_e_sim = evt.p_e_sim;
    double rel_delta_e = ((p_e - p_e_sim) / p_e_sim);
    // double rel_delta_e = ((p_e - p_e_sim) / p_e_sim);

    sn_hist->Fill(he_theta_diff);
    if (cuts_.IsThetaCut(evt)) {
      signal_theta_vs_e_hist->Fill(p_theta, p_e);
      signal_1d_e_hist->Fill(rel_delta_e);
    } else {
      bg_theta_vs_e_hist->Fill(p_theta, p_e);
      bg_1d_e_hist->Fill(rel_delta_e);
    }
  }

  virtual void fill_e_cut_hists(Evt& evt) {
    if (!is_well_defined(evt.p_theta)
        || !is_well_defined(evt.he_theta)) {
      return;
    }
    double xpos = cuts_.IsLeftEspriEvent(evt) ? evt.esl_xpos : evt.esr_xpos;
    if (std::abs(xpos) > 9000) {
      return;
    }
    if (( cuts_.IsLeftEspriEvent(evt) && evt.esl_de_raw < 240 ) ||
        ( cuts_.IsRightEspriEvent(evt) && evt.esr_de_raw < 240 )) {
      return;
    }

    double he_theta_diff = evt.he_theta - evt.he_theta_theor;
    // double p_theta = evt.p_theta_eff;
    double p_e = evt.p_e;
    double p_e_sim = evt.p_e_sim;
    double rel_delta_e = ((p_e - p_e_sim) / p_e_sim);
    ElaCuts loc_cuts = cuts_;

    for (size_t i = 0; i < e_cut_widths.size(); i++) {
      loc_cuts.espri_e_width_stdev = e_cut_widths[i];
      e_cut_hists[i]->Fill(rel_delta_e);
      e_cut_sn_tot_hists[i]->Fill(he_theta_diff);
      if (loc_cuts.IsEspriECut(evt)) {
        e_cut_sn_in_hists[i]->Fill(he_theta_diff);
      } else {
        e_cut_sn_out_hists[i]->Fill(he_theta_diff);
      }
    }
  }

};


void draw_cut_marker_lines(double cut_width, TH1* hist) {
  s13::ana::draw_marker_line_ysym(cut_width / 2,
                                  hist->GetMinimum(),
                                  hist->GetMaximum());
}

void
draw_e_spectra(s13::ana::RootScript& script, std::shared_ptr<ESpectra> alg) {
  script.NewPage(1,1).cd();
  alg->sn_hist->Draw();
  draw_cut_marker_lines(alg->cuts().theta_width, alg->sn_hist);

  auto pe_graph =
    s13::ana::make_proton_angle_energy_graph(cli_opts.dataset_type(),
                                             alg->cuts().espri_selector);
  script.NewPage(2,1).cd();
  alg->signal_theta_vs_e_hist->Draw("colz");
  gPad->SetLogz();
  pe_graph->Draw("SAME");
  script.cd();
  alg->bg_theta_vs_e_hist->Draw("colz");
  gPad->SetLogz();

  script.NewPage(2,1).cd();
  alg->signal_1d_e_hist->Fit("gaus", "", "");
  alg->signal_1d_e_hist->Draw();
  script.cd();
  alg->bg_1d_e_hist->Draw();
}

void
draw_e_cut_sn(s13::ana::RootScript& script, std::shared_ptr<ESpectra> alg) {
  for (size_t i = 0; i < alg->e_cut_widths.size(); i++) {
    if (i % 3 == 0) {
      script.NewPage(3,2);
    }

    script.cd();
    alg->e_cut_hists[i]->Draw();
    draw_cut_marker_lines(alg->e_cut_widths[i], alg->e_cut_hists[i]);

    script.cd();
    auto s2 = s13::ana::draw_thstack({alg->e_cut_sn_in_hists[i],
          alg->e_cut_sn_out_hists[i], alg->e_cut_sn_tot_hists[i]},
      "1-d kin. corrs.", "NOSTACK");
    s2->SetMaximum(s2->GetMaximum() * 0.8);
    gPad->BuildLegend(0.15, 0.65, 0.40, 0.85);
  }
}

void e_spectra() {
  auto cuts = s13::io::parse_cuts_config(std::fstream(cli_opts.cuts_config()));
  cuts.use_espri_aang_corr = true;
  auto alg =
    s13::ana::walk_alg<ESpectra>(cuts,
                                 init_dataset_from_cli(cli_opts),
                                 cli_opts.use_mt());

  s13::ana::RootScript script(cli_opts.output_file_basename());
  gStyle->SetOptStat(1111111);
  gStyle->SetOptFit();
  draw_e_spectra(script, alg);
  draw_e_cut_sn(script, alg);
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

  e_spectra();
}
