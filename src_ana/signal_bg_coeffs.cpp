
#include <iostream>
#include <fstream>

#include "commonalg.hpp"
#include "treewalker.hpp"
#include "cli.hpp"
#include "init_ds.hpp"
#include "consts.hpp"
#include "drawing.hpp"
#include "cuts_conf.hpp"
#include "rootscript.hpp"


using namespace s13::ana;

static s13::misc::CliOptions cli_opts;


class OutPhiSignalBgCoeffs : public s13::ana::SimpleTTreeAlgorithmBase {

private:

  using Yield = s13::ana::ScatteringYield;
  using Cuts = s13::ana::ElasticScatteringCuts;
  using Espri = s13::ana::Espri;
  using DatasetType = s13::ana::DatasetType;
  using PolDir = s13::ana::PolarizationDirection;
  using Dsdcs = s13::ana::Dsdcs;
  using Evt = s13::ana::ScatteringEvent;

  static const int k_yaxis_limit_he4 = 150;
  static const int k_yaxis_limit_he6 = 150;

public:

  OutPhiSignalBgCoeffs(DatasetType ds_type,
                       Cuts cuts,
                       double theta_range_start = 50,
                       double theta_range_end = 70,
                       size_t bin_num = 15) :
    ds_type_{ds_type},
    cuts_{cuts},
    theta_range_start_{theta_range_start},
    theta_range_end_{theta_range_end},
    bin_num_{bin_num},
    logger_{"OutPhiSignalBgCoeffs"}
  {
    init();
  }

  double bin_width() {
    return (theta_range_end_ - theta_range_start_) / bin_num_;
  }

  double bin_center(int bin_idx) {
    return theta_range_start_ + bin_width() / 2 + bin_width() * bin_idx;
  }

  double bin_min_theta(int bin_idx) {
    return bin_center(bin_idx) - bin_width() / 2;
  }

  double bin_max_theta(int bin_idx) {
    return bin_center(bin_idx) + bin_width() / 2;
  }

  virtual void process_event(Evt& evt) {
    compute_bg_coeffs(evt);
  }

  virtual void finalize() {
    logger_.info("Computed phi distributions histograms, fitting results...");
    for (size_t i = 0; i < hist_phi_corrs_total.size(); ++i) {
      FitParams fit_params_bg = fit_phi_corr(hist_phi_corrs_out_cut.at(i));
      FitParams fit_params_signal = fit_phi_corr(hist_phi_corrs_in_cut.at(i));

      Double_t fit_params_overall[6];
      std::copy(fit_params_bg.params.begin(), fit_params_bg.params.end(),
                fit_params_overall);
      std::copy(fit_params_signal.params.begin(), fit_params_signal.params.end(),
                fit_params_overall + 3);
      auto total_fit_fn = new TF1("mstotal", "gaus(0)+gaus(3)",
                                  ElasticScatteringCuts::k_phi_bg_min_angle,
                                  ElasticScatteringCuts::k_phi_bg_max_angle);
      total_fit_fn->SetParameters(fit_params_overall);
      FitParams fit_params_total = fit_phi_corr_total(hist_phi_corrs_total.at(i),
                                                      total_fit_fn);

      // write fit params of phi corr to output files
      signal_fit_params_file_ << bin_center(i) << "," << fit_params_signal;
      bg_fit_params_file_ << bin_center(i) << "," << fit_params_bg;
      total_fit_params_file_ << bin_center(i) << "," << fit_params_total;

      // compute and write normalization factor
      fit_params_bg = FitParams(fit_params_total, 0 , 3);
      fit_params_signal = FitParams(fit_params_total, 3, 6);
      Double_t r_norm = compute_bg_norm_factor(bin_center(i),
                                               fit_params_bg, fit_params_signal);
      r_norm_file_ << bin_center(i) << "," << r_norm << std::endl;
    }

    bg_fit_params_file_.close();
    signal_fit_params_file_.close();
    total_fit_params_file_.close();
    r_norm_file_.close();
  }

  std::vector<TH1*> hist_phi_corrs_total;
  std::vector<TH1*> hist_phi_corrs_in_cut;
  std::vector<TH1*> hist_phi_corrs_out_cut;

private:

  std::string out_fname_base() {
    std::string fname_prefix =
      ds_type_ == s13::ana::DatasetType::he6 ? "he6" : "he4";
    fname_prefix += cuts_.use_var_phi_cut ? "_varphi" : "_constphi";
    if (!cuts_.use_var_phi_cut) {
      TString phiWstr = s13::ana::tstrfmt("_phiW%.1f", cuts_.phi_width);
      fname_prefix += phiWstr;
    }
    fname_prefix += cuts_.dsdc_selector == Dsdcs::fdc0 ? "_fdc0" : "_s1dc";
    if (cuts_.espri_selector == Espri::left) {
      fname_prefix += "_esl";
    } else if (cuts_.espri_selector == Espri::right) {
      fname_prefix += "_esr";
    }

    return "out/" + fname_prefix;
  }

  std::vector<TH1*> init_phicorr_hists(TString title) {
    std::vector<TH1*> res;

    for (size_t i = 0; i < bin_num_; ++i) {
      TString tit =
        TString::Format("%.0f<p_{#theta}<%.0f (%s)",
                        bin_min_theta(i), bin_max_theta(i),
                        title.Data());
      auto hist =
        s13::ana::make_th1(300, 120, 240, tit,
                           "p_{#phi} - He_{#phi} [deg]", "Counts");
      res.push_back(hist);
    }

    return res;
  }

  void init() {
    bg_fit_params_file_.open(out_fname_base() + "_phi_corr_bg_fit_params.csv");
    signal_fit_params_file_.open(out_fname_base() + "_phi_corr_signal_fit_params.csv");
    total_fit_params_file_.open(out_fname_base() + "_phi_corr_total_fit_params.csv");
    r_norm_file_.open(out_fname_base() + "_total_to_bg_norm_factor.csv");

    r_norm_file_ << "# computed at phi cut width: "
                 << cuts_.phi_width << std::endl;
    signal_fit_params_file_ << "# computed at phi cut width: "
                            << cuts_.phi_width << std::endl;
    bg_fit_params_file_ << "# computed at phi cut width: "
                        << cuts_.phi_width << std::endl;
    total_fit_params_file_ << "# computed at phi cut width: "
                           << cuts_.phi_width << std::endl;

    hist_phi_corrs_total = init_phicorr_hists("TOT");
    hist_phi_corrs_in_cut = init_phicorr_hists("IN cut");
    hist_phi_corrs_out_cut = init_phicorr_hists("OUT cut");
  }

  int event_bin_idx(Evt& evt) {
    size_t idx = (evt.p_theta_eff - theta_range_start_) / bin_width();
    if (idx < bin_num_) {
      return idx;
    } else {
      return -1;
    }
  }

  void fill_in_cut_hists(Evt& evt) {
    int bin_idx = event_bin_idx(evt);
    if (bin_idx == -1) {
      return;
    }

    bool do_phi_cut = cuts_.apply_phi_cut;
    cuts_.apply_phi_cut = true;
    if (!cuts_.IsSignalEvent(evt)) {
      return;
    }
    cuts_.apply_phi_cut = do_phi_cut;

    hist_phi_corrs_in_cut.at(bin_idx)->Fill(std::abs(evt.p_phi - cuts_.GetHePhi(evt)));
  }

  void fill_out_cut_hists(Evt& evt) {
    int bin_idx = event_bin_idx(evt);
    if (bin_idx == -1) {
      return;
    }

    if (!cuts_.IsBackgroundEvent(evt)) {
      return;
    }

    hist_phi_corrs_out_cut.at(bin_idx)->Fill(std::abs(evt.p_phi - cuts_.GetHePhi(evt)));
  }

  void fill_total_hists(Evt& evt) {
    int bin_idx = event_bin_idx(evt);
    if (bin_idx == -1) {
      return;
    }

    bool do_phi_cut = cuts_.apply_phi_cut;
    cuts_.apply_phi_cut = false;
    if (!cuts_.IsSignalEvent(evt)) {
      return;
    }
    cuts_.apply_phi_cut = do_phi_cut;

    hist_phi_corrs_total.at(bin_idx)->Fill(std::abs(evt.p_phi - cuts_.GetHePhi(evt)));
  }

  void compute_bg_coeffs(Evt& evt) {
    fill_in_cut_hists(evt);
    fill_out_cut_hists(evt);
    fill_total_hists(evt);
  }

  FitParams fit_phi_corr(TH1* hist) {
    FitParams fit_params{3};
    const Double_t *fit_params_raw;
    const Double_t *fit_errors_raw;

    if (hist->GetEntries() == 0) {
      std::cout << "WARN: Empty histogram passed for fitting!" << std::endl;
      const Double_t zero_params[] = {0,0,0};
      fit_params.SetParams(zero_params);
      fit_params.SetErrors(zero_params);
      return fit_params;
    }

    hist->Fit("gaus", "W");
    fit_params_raw = hist->GetFunction("gaus")->GetParameters();
    fit_errors_raw = hist->GetFunction("gaus")->GetParErrors();

    fit_params.SetParams(fit_params_raw);
    fit_params.SetErrors(fit_errors_raw);
    return fit_params;
  }

  FitParams fit_phi_corr_total(TH1* hist, TF1* fit_fn) {
    FitParams fit_params{6};
    const Double_t *fit_params_raw;
    const Double_t *fit_errors_raw;

    if (hist->GetEntries() == 0) {
      std::cout << "WARN: Empty histogram passed for fitting!" << std::endl;
      const Double_t zero_params[] = {0,0,0,0,0,0};
      fit_params.SetParams(zero_params);
      fit_params.SetErrors(zero_params);
      return fit_params;
    }

    hist->Fit(fit_fn, "R");
    fit_params_raw = hist->GetFunction(fit_fn->GetName())->GetParameters();
    fit_errors_raw = hist->GetFunction(fit_fn->GetName())->GetParErrors();

    fit_params.SetParams(fit_params_raw);
    fit_params.SetErrors(fit_errors_raw);
    return fit_params;
  }

  /*
    Function computes normalization factor to estimate contribution of background
    to the elastic scattering peak. R is defined as follows:

    R = Y_bg_tot / Y_bg_out

    During asymmetry/CS computation to get yield of background events when
    R_tot_in is known, just:

    Y_bg_in = Y_bg_out * (R - 1)

    Then elastic yield will be given by:

    Y_ela_in = Y_tot_in - Y_bg_in
  */
  Double_t compute_bg_norm_factor(Double_t p_theta,
                                  FitParams& fit_params_bg,
                                  FitParams& fit_params_signal) {
    Double_t phi_max = Cuts::k_phi_bg_max_angle;
    Double_t phi_min = Cuts::k_phi_bg_min_angle;
    Double_t phi_corr_center = Cuts::k_phi_dist_center_deg;

    TF1 bg_fit_fn("bg_fit_fn", "gaus(0)", phi_min, phi_max);
    bg_fit_fn.SetParameters(&fit_params_bg.params[0]);
    TF1 signal_fit_fn("signal_fit_fn", "gaus(0)", phi_min, phi_max);
    signal_fit_fn.SetParameters(&fit_params_signal.params[0]);

    Double_t a_lim = phi_corr_center - cuts_.phi_width / 2.;
    Double_t b_lim = phi_corr_center + cuts_.phi_width / 2.;
    if (cuts_.use_var_phi_cut) {
      Double_t phi_cut_sigma = cuts_.GetVarPhiCutSigma(p_theta);
      a_lim = phi_corr_center - 3 * phi_cut_sigma;
      b_lim = phi_corr_center + 3 * phi_cut_sigma;
    }
    Double_t bg_total_int = bg_fit_fn.Integral(phi_min, phi_max);
    Double_t bg_in_int = bg_fit_fn.Integral(a_lim, b_lim);
    Double_t bg_out_int = bg_total_int - bg_in_int;

    return bg_total_int / bg_out_int;
  }

  DatasetType ds_type_;
  Cuts cuts_;
  double theta_range_start_;
  double theta_range_end_;
  size_t bin_num_;

  std::ofstream bg_fit_params_file_;
  std::ofstream signal_fit_params_file_;
  std::ofstream total_fit_params_file_;
  std::ofstream r_norm_file_;

  s13::misc::MessageLogger logger_;
};


using Script = s13::ana::RootScript;
using Cuts = s13::ana::ElasticScatteringCuts;
using CoeffsAlgPtr = std::shared_ptr<OutPhiSignalBgCoeffs>;


void draw_signal_bg_coeffs(Script& script,
                           CoeffsAlgPtr alg) {
  for (size_t i = 0; i < alg->hist_phi_corrs_total.size(); ++i) {
    if (i % 4 == 0) {
      script.NewPage(4, 3).cd(s13::ana::PadSeq::column);
    }

    script.cd();
    alg->hist_phi_corrs_in_cut.at(i)->Draw();
    script.cd();
    alg->hist_phi_corrs_out_cut.at(i)->Draw();
    script.cd();
    alg->hist_phi_corrs_total.at(i)->Draw();
  }
}

CoeffsAlgPtr walk_tree(const Cuts& cuts,
                       DatasetType ds_type,
                       const std::vector<TChain*>& dataset,
                       bool use_mt = false) {
  using PtrTreeWalker = std::shared_ptr<s13::ana::ScatteringTreeWalker>;
  s13::misc::MessageLogger logger("s13::ana::walk_alg()");
  CoeffsAlgPtr alg(new OutPhiSignalBgCoeffs(ds_type, cuts,
                                            55, 71, 17));
  auto tree_walker = PtrTreeWalker(new ScatteringTreeWalker());
  tree_walker->add(alg);

  if (!use_mt) {
    logger.debug("Starting algorithm walk (single thread)...");
    tree_walker->walk(*dataset[0]);
  } else {
    logger.debug("Starting algorithm walk (multiple threads)...");
    s13::ana::ParallelTreeWalker ptree_walker(tree_walker, dataset);
    ptree_walker.walk();
  }

  return alg;
}
void signal_bg_coeffs() {
  auto cuts = s13::io::parse_cuts_config(std::fstream(cli_opts.cuts_config()));
  auto alg = walk_tree(cuts,
                       cli_opts.dataset_type(),
                       init_dataset_from_cli(cli_opts),
                       cli_opts.use_mt());

  s13::ana::RootScript script("signal_bg_coeffs");
  gStyle->SetOptStat(1111111);
  // gStyle->SetOptFit();
  draw_signal_bg_coeffs(script, alg);
}

int main(int argc, char** argv) {
  s13::misc::MessageLogger logger("main()");
  cli_opts.require_cuts_conf_opt(true);
  try {
    cli_opts.parse_options(argc, argv);
  } catch (std::exception& ex) {
    logger.error("Invalid argument provided: %s", ex.what());
    return -1;
  }

  if (cli_opts.dataset_type() == DatasetType::carbon) {
    logger.error("Carbon cannot be specified for this executable only he4|he6");
    exit(1);
  }

  signal_bg_coeffs();
}
