
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


class FwAngBg : public s13::ana::SimpleTTreeAlgorithmBase {

private:
  s13::misc::MessageLogger logger_;

  ElaCuts cuts_;

public:
  TH1* x_hist;
  TH1* y_hist_problem_rgn;
  TH1* y_hist_good;
  TH1* theta_hist;
  TH2* thetay_hist;
  TH2* xy_hist;

  TH2* total_thetas_hist;
  TH1* total_kin_corr_hist;
  TH2* bw_noise_thetas_hist;
  TH1* bw_noise_kin_corr_hist;

  TH2* thetae_hist;
  TH1* rel_e_hist;
  TH1* vertz_hist;

public:
  FwAngBg(const ElaCuts& cuts, int bin_num = 7) :
    logger_{"FwAngBg"}, cuts_{cuts} {
      assert(cuts_.espri_selector != Espri::both
             && "Can work only on L/R ESPRI separately!");
      init_hists();
    }

  FwAngBg(const FwAngBg& other) = delete;
  FwAngBg(FwAngBg&&) = delete;
  FwAngBg operator=(const FwAngBg&) = delete;
  FwAngBg operator=(FwAngBg&&) = delete;

  virtual void setup() {
  }

  virtual void process_event(Evt& evt) {
    if (cuts_.IsSignalEvent(evt)) {
      fill_hists(evt);
    }
  }

  virtual void finalize() {
  }

  ElaCuts cuts() const {
    return cuts_;
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

  void init_hists() {
    xy_hist = s13::ana::make_th2(200, -280, 280, 200, -280, 280, "RDC X/Y", "X [mm]", "Y [mm]");
    thetay_hist = s13::ana::make_th2(200, 50, 75, 200, -280, 280, "#theta_p vs. RDC Y",
                                     "#theta_p [lab. deg.]", "Y [mm]");
    x_hist = s13::ana::make_th1(65, -100, 100, "RDC X", "X [mm]", "Counts");
    y_hist_problem_rgn = s13::ana::make_th1(150, -220, 220,
                                            "RDC Y (62.5-63.5)", "Y [mm]", "Counts");
    y_hist_good = s13::ana::make_th1(150, -220, 220,
                                     "RDC Y (64.5-65.5)", "Y [mm]", "Counts");
    theta_hist = s13::ana::make_th1(60, 58, 68, "Proton theta", "#theta_p [lab. deg.]", "Counts");

    rel_e_hist = s13::ana::make_th1(200, -1, 1, "BW noise - proton E_rel", "Counts");
    vertz_hist = s13::ana::make_th1(200, -200, 200, "BW noise - vertZ", "Counts");
    thetae_hist = s13::ana::make_th2(200, 50, 75, 200, 0, 150,
                                     "BW noise - proton #theta vs. E",
                                     "#theta_{p}", "E [MeV]");
    total_thetas_hist = init_2d_kin_corr_hist("Total - kin corr");
    total_kin_corr_hist = init_kin_corr_hist("Total - kin corr");
    bw_noise_thetas_hist = init_2d_kin_corr_hist("BW noise - kin corr");
    bw_noise_kin_corr_hist = init_kin_corr_hist("BW noise - kin corr");
  }

  virtual void fill_hists(Evt& evt) {
    if (!is_well_defined(evt.p_theta)
        || !is_well_defined(evt.he_theta)) {
      return;
    }

    Espri espri = evt.IsLeftEspriEvent() ? Espri::left : Espri::right;
    double xpos = cuts_.IsLeftEspriEvent(evt) ? evt.esl_xpos : evt.esr_xpos;
    double ypos = cuts_.IsLeftEspriEvent(evt) ? evt.esl_ypos : evt.esr_ypos;
    double aang = cuts_.GetCorrectedEspriAang(evt, espri);
    // double aang = cuts_.IsLeftEspriEvent(evt) ? evt.esl_aang : evt.esr_aang;
    double e_rel = cuts_.ComputeCorrectedERel(evt);
    double vertz = s13::ana::compute_vertz(xpos, aang, espri);
    double p_theta = evt.p_theta_eff;
    double he_theta = cuts_.dsdc_selector == Dsdcs::s1dc ? evt.s1dc_theta : evt.fdc0_theta;
    double he_theta_diff = he_theta - evt.he_theta_theor;

    if (std::abs(xpos) > 9000) {
      return;
    }

    xy_hist->Fill(xpos, ypos);
    thetay_hist->Fill(p_theta, ypos);
    x_hist->Fill(xpos);
    if (p_theta >= 64.5 && p_theta <= 65.5) {
      y_hist_good->Fill(ypos);
    }
    if (p_theta >= 62.5 && p_theta <= 63.5) {
      y_hist_problem_rgn->Fill(ypos);
    }
    theta_hist->Fill(p_theta);
    thetae_hist->Fill(p_theta, evt.p_e);

    total_kin_corr_hist->Fill(he_theta_diff);
    total_thetas_hist->Fill(he_theta, p_theta);

    if ((evt.IsLeftEspriEvent() && xpos > 160) ||
        (evt.IsRightEspriEvent() && xpos < -160)) {
      vertz_hist->Fill(vertz);
      rel_e_hist->Fill(e_rel);
      bw_noise_kin_corr_hist->Fill(he_theta_diff);
      bw_noise_thetas_hist->Fill(he_theta, p_theta);
    }

  }

};

void
draw_hists(s13::ana::RootScript& script, std::shared_ptr<FwAngBg> alg) {
  // alg->theta_hist->Fit("gaus", "", "", 51.2, 55);

  script.NewPage(2,2).cd(s13::ana::PadSeq::column);
  // script.cd();
  // alg->xy_hist->Draw("colz");
  // script.cd();
  // alg->thetay_hist->Draw("colz");

  script.cd();
  alg->x_hist->Draw();
  gPad->SetLogy();
  script.cd();
  alg->y_hist_good->Draw();
  gPad->SetLogy();
  script.cd();
  alg->y_hist_problem_rgn->Draw();
  gPad->SetLogy();
  script.cd();
  alg->theta_hist->Draw();
  gPad->SetLogy();

  script.NewPage(2,2).cd(s13::ana::PadSeq::column);
  script.cd();
  alg->thetae_hist->Draw("colz");
  script.cd();
  alg->rel_e_hist->Draw();
  script.cd();
  alg->vertz_hist->Draw();
  gPad->SetLogy();

  script.NewPage(2,1);
  script.cd();
  // auto s1 = s13::ana::draw_thstack({alg->bw_noise_thetas_hist, alg->total_thetas_hist},
  //   "TOT & BW noise kin corr", "NOSTACK");
  // s1->SetMaximum(alg->total_thetas_hist->GetMaximum());
  alg->total_thetas_hist->Draw("");
  alg->bw_noise_thetas_hist->Draw("col SAME");
  // gPad->BuildLegend(0.65, 0.8, 0.95, 0.95);
  script.cd();
  auto s2 = s13::ana::draw_thstack({alg->bw_noise_kin_corr_hist, alg->total_kin_corr_hist},
                                   "TOT & BW noise kin corr", "NOSTACK");
  s2->SetMaximum(alg->total_kin_corr_hist->GetMaximum());
  gPad->BuildLegend(0.65, 0.8, 0.95, 0.95);

  std::cout << "Cut definition: " << alg->cuts().GetTotalCutTitle() << std::endl;
}

void fw_ang_bg() {
  auto cuts = s13::io::parse_cuts_config(std::fstream(cli_opts.cuts_config()));
  cuts.use_espri_aang_corr = true;
  auto alg =
    s13::ana::walk_alg2<FwAngBg>(init_dataset_from_cli(cli_opts),
                                 cli_opts.use_mt(),
                                 cuts, 7);
  auto script_fname =
    s13::ana::tstrfmt("%s_%s", cli_opts.output_file_basename().data(),
                      s13::ana::espri_selector_to_str(cuts.espri_selector).data());
  s13::ana::RootScript script(script_fname);
  gStyle->SetOptStat(1111111);
  // gStyle->SetOptFit();
  draw_hists(script, alg);
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

  fw_ang_bg();
}
