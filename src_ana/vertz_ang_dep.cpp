

#include <functional>

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


class CarbonYields : public s13::ana::SimpleTTreeAlgorithmBase {

public:

  static constexpr double k_min_theta = s13::gk_cs_range_start;
  static constexpr double k_max_theta = s13::gk_cs_range_end;

  // cut proc alias
  using CutProc = bool (ElaCuts::*)(const Evt& evt) const;
  // in & out components analysis on ang dependence
  using YieldData = std::vector< std::vector<TH1*> >;

  YieldData yields_in;
  YieldData yields_out;
  YieldData yields_tot;

  bool scale_carbon_bg = true;

private:

  s13::misc::MessageLogger logger_;

  std::vector<ElaCuts> cuts_range_;
  CutProc cut_proc_;

  size_t bin_num_;

public:

  CarbonYields(const std::vector<ElaCuts>& cuts_range,
               CutProc cut_proc,
               size_t bin_num = s13::gk_cs_bin_number) :
    logger_{"CarbonYields"}, cuts_range_{cuts_range}, cut_proc_{cut_proc},
    bin_num_{bin_num} {
      init_hists();
    }

  CarbonYields(const CarbonYields&) = delete;
  CarbonYields(CarbonYields&&) = delete;
  CarbonYields operator=(const CarbonYields&) = delete;
  CarbonYields operator=(CarbonYields&&) = delete;

  virtual void setup() {
  }

  virtual void process_event(Evt& evt) {
    assert(cuts_range_.size() > 0);
    if (cuts_range_[0].IsSignalEvent(evt)) {
      fill_carbon_yield_hists(evt);
    }
  }

  virtual void finalize() {
    if (scale_carbon_bg) {
      for (size_t i = 0; i < yields_in.size(); i++) {
        auto& v_in = yields_in[i];
        auto& v_out = yields_out[i];
        auto& v_tot = yields_tot[i];
        for (size_t c = 0; c < v_in.size(); c++) {
          v_in[c]->Scale(s13::gk_carbon_bg_norm_factor);
          v_out[c]->Scale(s13::gk_carbon_bg_norm_factor);
          v_tot[c]->Scale(s13::gk_carbon_bg_norm_factor);
        }
      }
    }
  }

  const std::vector<ElaCuts>& cuts_range() {
    return cuts_range_;
  }

  size_t bin_num() {
    return bin_num_;
  }

  double bin_width() {
    return (k_max_theta - k_min_theta) / bin_num_;
  }

  double bin_min_theta(int idx) {
    double minth = k_min_theta + bin_width() * idx;
    return minth;
  }

  double bin_max_theta(int idx) {
    double maxth = k_min_theta + bin_width() * (idx + 1);
    return maxth;
  }

  double bin_center_theta(int idx) {
    return (bin_min_theta(idx) + bin_max_theta(idx)) / 2;
  }

  int bin_idx(Evt& evt) {
    if (evt.p_theta_eff - k_min_theta < 0) {
      return -1;
    }

    size_t idx = (evt.p_theta_eff - k_min_theta) / bin_width();
    if (idx < bin_num_) {
      return idx;
    } else {
      return -1;
    }
  }

private:

  TH1* init_kin_corr_hist(TString title) {
    auto th1 = &s13::ana::make_th1;
    auto hist = th1(80, -8, 8, title, "#theta_{exp} - #theta_{theor}", "Counts");
    return hist;
  }

  void init_hists() {
    for (size_t i = 0; i < bin_num_; i++) {
      yields_in.push_back(std::vector<TH1*>());
      yields_out.push_back(std::vector<TH1*>());
      yields_tot.push_back(std::vector<TH1*>());
      for (size_t c = 0; c < cuts_range_.size(); c++) {
        yields_in[i].push_back(init_kin_corr_hist("IN"));
        yields_out[i].push_back(init_kin_corr_hist("OUT"));
        yields_tot[i].push_back(init_kin_corr_hist("TOT"));
      }
    }
  }

  void fill_carbon_yield_hists(Evt& evt) {
    Double_t he_theta = cuts_range()[0].dsdc_selector == Dsdcs::s1dc
      ? evt.s1dc_theta : evt.fdc0_theta;
    double p_theta = evt.p_theta_eff;
    double he_theta_diff = he_theta - evt.he_theta_theor;

    if (!is_well_defined(p_theta)
        || !is_well_defined(he_theta) ||
        !evt.IsEspriEvent()) {
      return;
    }

    int evt_bin_idx =  bin_idx(evt);
    if (evt_bin_idx == -1) {
      return;
    }

    for (size_t c = 0; c < cuts_range_.size(); c++) {
      auto lcuts = cuts_range_[c];
      yields_tot[evt_bin_idx][c]->Fill(he_theta_diff);
      if ((lcuts.*cut_proc_)(evt)) {
        yields_in[evt_bin_idx][c]->Fill(he_theta_diff);
      } else {
        yields_out[evt_bin_idx][c]->Fill(he_theta_diff);
      }
    }
  }

};


class VertzAngDep : public s13::ana::SimpleTTreeAlgorithmBase {

public:

  static constexpr double k_min_theta = s13::gk_cs_range_start;
  static constexpr double k_max_theta = s13::gk_cs_range_end;

  std::vector<TH1*> aang_hists;
  std::vector<TH1*> vertz_hists;
  std::vector<TH1*> vertz_sn_in_hists;
  std::vector<TH1*> vertz_sn_out_hists;
  std::vector<TH1*> vertz_sn_tot_hists;

  // in & out components analysis on ang dependence
  using AngularSnData = std::vector< std::vector<TH1*> >;

  std::vector<TH1*> ang_dep_vertz_dists;
  AngularSnData ang_dep_sn_in;
  AngularSnData ang_dep_sn_out;
  AngularSnData ang_dep_sn_tot;
  AngularSnData ang_dep_outphi;

private:

  s13::misc::MessageLogger logger_;

  std::vector<ElaCuts> cuts_range_;

  size_t bin_num_;

public:

  VertzAngDep(const std::vector<ElaCuts>& cuts,
              size_t bin_num = s13::gk_cs_bin_number) :
  logger_{"VertzAngDep"}, cuts_range_{cuts}, bin_num_{bin_num} {
    init_hists();
  }

  VertzAngDep(const VertzAngDep&) = delete;
  VertzAngDep(VertzAngDep&&) = delete;
  VertzAngDep operator=(const VertzAngDep&) = delete;
  VertzAngDep operator=(VertzAngDep&&) = delete;

  virtual void setup() {
  }

  virtual void process_event(Evt& evt) {
    assert(!cuts_range_.empty());
    if (cuts_range_[0].IsSignalEvent(evt)) {
      fill_aang_and_vertz_bin_hists(evt);
    }
    fill_vertz_ang_dep_hists(evt);
  }

  virtual void finalize() {
    auto bg_norm_fact_interp =
      s13::ana::make_bg_norm_factor_interp2(DsType::he4,
                                            cuts_range_[0].dsdc_selector);
    for (size_t i = 0; i < ang_dep_outphi.size(); i++) {
      auto& v = ang_dep_outphi.at(i);
      for (auto& h : v) {
        double theta = bin_center_theta(i);
        double phiW = cuts_range_[0].GetPhiCutWidth(theta);
        h->Scale(bg_norm_fact_interp.eval(theta, phiW) - 1);
      }
    }
  }

  const std::vector<ElaCuts>& cuts_range() {
    return cuts_range_;
  }

  size_t bin_num() {
    return bin_num_;
  }

  double bin_width() {
    return (k_max_theta - k_min_theta) / bin_num_;
  }

  double bin_min_theta(int idx) {
    double minth = k_min_theta + bin_width() * idx;
    return minth;
  }

  double bin_max_theta(int idx) {
    double maxth = k_min_theta + bin_width() * (idx + 1);
    return maxth;
  }

  double bin_center_theta(int idx) {
    return (bin_min_theta(idx) + bin_max_theta(idx)) / 2;
  }

  void write_aang_fit_params(std::string fname) {
    std::ofstream ofs(fname);
    ofs << "# p_theta,gaus_const,gaus_mean,gaus_sigma" << std::endl;
    for (size_t i = 0; i < bin_num_; ++i) {
      double theta = bin_center_theta(i);
      auto fit_fn = aang_hists[i]->GetFunction("gaus");
      assert(fit_fn != nullptr);
      double* params = aang_hists[i]->GetFunction("gaus")->GetParameters();
      ofs << theta << "," << params[0] << "," << params[1] << ","
          << params[2] << std::endl;
    }
  }

  void write_vertz_fit_params(std::string fname) {
    std::ofstream ofs(fname);
    ofs << "# p_theta,gaus_const,gaus_mean,gaus_sigma" << std::endl;
    for (size_t i = 0; i < bin_num_; ++i) {
      double theta = bin_center_theta(i);
      auto fit_fn = vertz_hists[i]->GetFunction("gaus");
      assert(fit_fn != nullptr);
      double* params = vertz_hists[i]->GetFunction("gaus")->GetParameters();
      ofs << theta << "," << params[0] << "," << params[1] << ","
          << params[2] << std::endl;
    }
  }

private:

  TH1* init_vertz_bin_hist(double bin_min_theta, double bin_max_theta,
                           double xlim = 350) {
    auto th1 = &s13::ana::make_th1;
    auto tit = s13::ana::tstrfmt("vertex Z pos. (%.2f<p_{#theta}<%.2f)",
                                 bin_min_theta, bin_max_theta);
    auto hist = th1(120, -xlim, xlim, tit, "Z [mm]", "Counts");
    return hist;
  }

  TH1* init_aang_bin_hist(double bin_min_theta, double bin_max_theta,
                          double xlim = 350) {
    auto th1 = &s13::ana::make_th1;
    auto tit = s13::ana::tstrfmt("#Delta aang (%.2f<p_{#theta}<%.2f)",
                                 bin_min_theta, bin_max_theta);
    auto hist = th1(200, -xlim, xlim, tit, "RDC aang [mrad]", "Counts");
    return hist;
  }

  TH1* init_aang_hist(TString title = "", double xlim = 1000) {
    auto th1 = &s13::ana::make_th1;
    auto tit = s13::ana::tstrfmt("#Delta aang. %s", title.Data());
    auto hist = th1(100, -xlim, xlim, tit, "RDC aang [mrad]", "Counts");
    return hist;
  }

  TH1* init_kin_corr_hist(TString title) {
    auto th1 = &s13::ana::make_th1;
    auto hist = th1(80, -8, 8, title, "#theta_{exp} - #theta_{theor}", "Counts");
    return hist;
  }

  TH2* init_2d_kin_corr_hist(TString title) {
    auto th2 = &s13::ana::make_th2;
    auto hist = th2(200, 0, 12, 200, 50, 75, title, "#theta_{He}", "#theta_{p}");
    return hist;
  }

  int bin_idx(Evt& evt) {
    if (evt.p_theta_eff - k_min_theta < 0) {
      return -1;
    }

    size_t idx = (evt.p_theta_eff - k_min_theta) / bin_width();
    if (idx < bin_num_) {
      return idx;
    } else {
      return -1;
    }
  }

  double correct_aang(double aang, double p_theta, Espri espri) {
    LinearInterp espri_computed_aang_mean{s13::gk_espri_computed_aang_a,
        s13::gk_espri_computed_aang_b};
    if (espri == Espri::left) {
      LinearInterp esl_aang_mean{s13::gk_espri_left_aang_a, s13::gk_espri_left_aang_b};
      auto corr = esl_aang_mean.eval(p_theta) - espri_computed_aang_mean.eval(p_theta);
      return aang - corr;
    } else if (espri == Espri::right) {
      LinearInterp esr_aang_mean{s13::gk_espri_right_aang_a, s13::gk_espri_right_aang_b};
      auto corr = esr_aang_mean.eval(p_theta) - espri_computed_aang_mean.eval(p_theta);
      return aang - corr;
    } else {
      throw std::logic_error("'BOTH' is an invalid ESPRI selector for this function!");
    }
  }

  void init_hists() {
    for (size_t i = 0; i < bin_num_; i++) {
      auto vz_hist = init_vertz_bin_hist(bin_min_theta(i),
                                         bin_max_theta(i), 300);
      auto aang_hist = init_aang_bin_hist(bin_min_theta(i),
                                          bin_max_theta(i), 300);
      vertz_sn_in_hists.push_back(init_kin_corr_hist("IN"));
      vertz_sn_out_hists.push_back(init_kin_corr_hist("OUT"));
      vertz_sn_tot_hists.push_back(init_kin_corr_hist("TOT"));
      vertz_hists.push_back(vz_hist);
      aang_hists.push_back(aang_hist);

      auto ad_vz_hist = init_vertz_bin_hist(bin_min_theta(i),
                                            bin_max_theta(i), 200);
      ang_dep_vertz_dists.push_back(ad_vz_hist);

      ang_dep_sn_in.push_back(std::vector<TH1*>());
      ang_dep_sn_out.push_back(std::vector<TH1*>());
      ang_dep_sn_tot.push_back(std::vector<TH1*>());
      ang_dep_outphi.push_back(std::vector<TH1*>());
      for (size_t c = 0; c < cuts_range_.size(); c++) {
        ang_dep_sn_in[i].push_back(init_kin_corr_hist("IN"));
        ang_dep_sn_out[i].push_back(init_kin_corr_hist("OUT"));
        ang_dep_sn_tot[i].push_back(init_kin_corr_hist("TOT"));
        ang_dep_outphi[i].push_back(init_kin_corr_hist("OUT_phi"));
      }
    }
  }

  void fill_aang_and_vertz_bin_hists(Evt& evt) {
    Double_t he_theta = cuts_range()[0].dsdc_selector == Dsdcs::s1dc
      ? evt.s1dc_theta : evt.fdc0_theta;
    double p_theta = evt.p_theta_eff;
    double he_theta_diff = he_theta - evt.he_theta_theor;

    if (!is_well_defined(p_theta)
        || !is_well_defined(he_theta) ||
        !evt.IsEspriEvent()) {
      return;
    }

    int evt_bin_idx =  bin_idx(evt);
    Espri espri = cuts_range_[0].IsLeftEspriEvent(evt) ? Espri::left : Espri::right;
    double xpos = cuts_range_[0].IsLeftEspriEvent(evt) ? evt.esl_xpos : evt.esr_xpos;

    if (evt_bin_idx == -1) {
      return;
    }

    double rdc_aang = espri == Espri::left ? evt.esl_aang : evt.esr_aang;
    double corr_aang = correct_aang(rdc_aang, evt.p_theta_eff, espri);
    double corr_vertz = s13::ana::compute_vertz(xpos, corr_aang, espri);
    double rdc_vertz = s13::ana::compute_vertz(xpos, rdc_aang, espri);
    double aang = cuts_range_[0].use_espri_aang_corr ? corr_aang : rdc_aang;
    double vertz = cuts_range_[0].use_espri_aang_corr ? corr_vertz : rdc_vertz;
    double tgt_aang = std::atan(xpos/s13::gk_dist_rdc_target) * 1e3;

    // we will not consider vertex Z positions far away from target chamber
    if (std::abs(vertz) > 300) {
      return;
    }

    vertz_hists[evt_bin_idx]->Fill(vertz);
    aang_hists[evt_bin_idx]->Fill(aang - tgt_aang);
    vertz_sn_tot_hists[evt_bin_idx]->Fill(he_theta_diff);
    if (cuts_range_[0].IsEspriVertzCut(evt)) {
      vertz_sn_in_hists[evt_bin_idx]->Fill(he_theta_diff);
    } else {
      vertz_sn_out_hists[evt_bin_idx]->Fill(he_theta_diff);
    }
  }

  void fill_vertz_ang_dep_hists(Evt& evt) {
    Double_t he_theta = cuts_range()[0].dsdc_selector == Dsdcs::s1dc
      ? evt.s1dc_theta : evt.fdc0_theta;
    double p_theta = evt.p_theta_eff;
    double he_theta_diff = he_theta - evt.he_theta_theor;

    if (!is_well_defined(p_theta)
        || !is_well_defined(he_theta) ||
        !evt.IsEspriEvent()) {
      return;
    }

    int evt_bin_idx =  bin_idx(evt);
    Espri espri = cuts_range_[0].IsLeftEspriEvent(evt) ? Espri::left : Espri::right;
    double xpos = cuts_range_[0].IsLeftEspriEvent(evt) ? evt.esl_xpos : evt.esr_xpos;

    if (std::abs(xpos) > 9000) {
      return;
    }
    if (( espri == Espri::left && evt.esl_de_raw < 240 ) ||
        ( espri == Espri::right && evt.esr_de_raw < 240 )) {
      return;
    }
    if (evt_bin_idx == -1) {
      return;
    }

    double rdc_aang = espri == Espri::left ? evt.esl_aang : evt.esr_aang;
    double corr_aang = correct_aang(rdc_aang, evt.p_theta_eff, espri);
    double corr_vertz = s13::ana::compute_vertz(xpos, corr_aang, espri);
    double rdc_vertz = s13::ana::compute_vertz(xpos, rdc_aang, espri);
    // double aang = cuts_range_.use_espri_aang_corr ? corr_aang : rdc_aang;
    double vertz = cuts_range_[0].use_espri_aang_corr ? corr_vertz : rdc_vertz;
    // double tgt_aang = std::atan(xpos/s13::gk_dist_rdc_target) * 1e3;

    for (size_t c = 0; c < cuts_range_.size(); c++) {
      auto loc_cut = cuts_range_[c];
      loc_cut.apply_espri_vertz_cut = false;
      if (loc_cut.IsSignalEvent(evt)) {
        ang_dep_vertz_dists[evt_bin_idx]->Fill(vertz);
        ang_dep_sn_tot[evt_bin_idx][c]->Fill(he_theta_diff);
        if (loc_cut.IsEspriVertzCut(evt)) {
          ang_dep_sn_in[evt_bin_idx][c]->Fill(he_theta_diff);
        } else {
          ang_dep_sn_out[evt_bin_idx][c]->Fill(he_theta_diff);
        }
      }

      if (loc_cut.IsBackgroundEvent(evt) && loc_cut.IsEspriVertzCut(evt)) {
        ang_dep_outphi[evt_bin_idx][c]->Fill(he_theta_diff);
      }
    }
  }

};

void draw_cut_marker_lines(double cut_width, TH1* hist) {
  s13::ana::draw_marker_line_ysym(cut_width / 2,
                                  hist->GetMinimum(),
                                  hist->GetMaximum());
}

double compute_vertz_bin_sn(double angle, TH1* hist, const ElaCuts& cuts) {
  s13::misc::MessageLogger logger("compute_vertz_bin_sn()");
  logger.debug("Computing vertZ bin S/N for bin at %.1f degrees...", angle);
  const double vertz_hw = cuts.espri_vertz_width / 2;
  const double bg_min1 = hist->GetXaxis()->GetXmin();
  const double bg_max1 = -vertz_hw;
  const double bg_min2 = vertz_hw;
  const double bg_max2 = hist->GetXaxis()->GetXmax();

  auto bin = [hist](double arg) -> int {
    return hist->GetXaxis()->FindBin(arg);
  };
  logger.debug("\tBG region (X arg): [%.1f, %.1f] - [%.1f, %.1f]",
               bg_min1, bg_max1, bg_min2, bg_max2);
  logger.debug("\tBG region (bin#): [%d, %d] - [%d, %d]",
               bin(bg_min1), bin(bg_max1), bin(bg_min2), bin(bg_max2) - 1);

  double signal = hist->Integral(bin(-vertz_hw), bin(vertz_hw));
  double bg = hist->Integral(bin(bg_min1), bin(bg_max1) );
  bg += hist->Integral(bin(bg_min2), bin(bg_max2) - 1);
  double sn = signal / bg;
  logger.debug("\tSignal/BG integrals: %.2f/%.2f; S/N: %.2f", signal, bg, sn);
  return sn;
}

double draw_bin_sn(double angle, TH1* bin_hist, const ElaCuts& cuts) {
  s13::misc::MessageLogger logger("draw_bin_sn()");
  double sigma = cuts.dsdc_selector == Dsdcs::s1dc ?
    s13::gk_theta_distribution_sigma_s1dc :
    s13::gk_theta_distribution_sigma_fdc0;
  const double bg_min_angle = -4;
  const double bg_max_angle = -3*sigma;

  double bg_int =
    bin_hist->Integral(bin_hist->GetXaxis()->FindBin(bg_min_angle),
                       bin_hist->GetXaxis()->FindBin(bg_max_angle));
  double signal_int =
    bin_hist->Integral(bin_hist->GetXaxis()->FindBin(-cuts.theta_width/2),
                       bin_hist->GetXaxis()->FindBin(cuts.theta_width/2));
  double sn = signal_int / bg_int;
  s13::ana::draw_text(s13::ana::tstrfmt("S/N: %.2f", sn));
  logger.debug("S/N for bin at %.2f deg.: %.2f (S: %.2f; N: %.2f)",
               angle, sn, signal_int, bg_int);

  double hist_ymin = bin_hist->GetMinimum();
  double hist_ymax = bin_hist->GetMaximum() * 0.7;
  s13::ana::draw_marker_line_ysym(cuts.theta_width / 2, hist_ymin, hist_ymax);
  s13::ana::draw_marker_line(bg_min_angle, hist_ymin, bg_min_angle, hist_ymax);
  s13::ana::draw_marker_line(bg_max_angle, hist_ymin, bg_max_angle, hist_ymax);
  return sn;
}

void draw_angdep_fits(s13::ana::RootScript& script, std::shared_ptr<VertzAngDep> alg) {
  for (size_t i = 0; i < alg->vertz_hists.size(); ++i) {
    if (i % 4 == 0) {
      script.NewPage(4, 2).cd(s13::ana::PadSeq::column);
    }

    script.cd();
    alg->aang_hists[i]->Draw();
    alg->aang_hists[i]->Fit("gaus", "", "", -60, 60);

    script.cd();
    alg->vertz_hists[i]->Draw();
    alg->vertz_hists[i]->Fit("gaus", "", "", -60, 60);
  }

  alg->write_aang_fit_params("out/vz_ang_dep_aang_fit_params.csv");
  alg->write_vertz_fit_params("out/vz_ang_dep_vertz_fit_params.csv");
}

void draw_angdep_sn2(s13::ana::RootScript& script,
                         std::shared_ptr<VertzAngDep> alg_vz,
                         std::shared_ptr<CarbonYields> alg_carbon) {
  for (size_t i = 0; i < alg_vz->ang_dep_vertz_dists.size(); i++) {
    script.NewPage(alg_vz->cuts_range().size(), 4);
    for (size_t c = 0; c < alg_vz->cuts_range().size(); c++) {
      script.cd();
      alg_vz->ang_dep_vertz_dists[i]->Draw();
      draw_cut_marker_lines(alg_vz->cuts_range()[c].espri_vertz_width,
                            alg_vz->ang_dep_vertz_dists[i]);
      auto cut_width_text =
        s13::ana::tstrfmt("vZ: +/-%.2f",
                          alg_vz->cuts_range()[c].espri_vertz_width/2);
      s13::ana::draw_text(cut_width_text);

      script.cd();
      auto s1 = s13::ana::draw_thstack({alg_vz->ang_dep_sn_tot[i][c],
            alg_vz->ang_dep_sn_out[i][c]},
        "TOT & OUT vertZ cut", "NOSTACK");
      s1->SetMaximum(alg_vz->ang_dep_sn_tot[i][c]->GetMaximum());
      gPad->BuildLegend(0.15, 0.65, 0.40, 0.85);

      if (cli_opts.dataset_type() == DsType::he6) {
        script.cd();
        alg_carbon->yields_in[i][c]->SetTitle("C IN");
        auto s2 = s13::ana::draw_thstack({alg_vz->ang_dep_sn_in[i][c],
              alg_carbon->yields_in[i][c]},
          "Aft. cuts vs. carbon BG", "NOSTACK");
        s2->SetMaximum(alg_vz->ang_dep_sn_in[i][c]->GetMaximum());
        gPad->BuildLegend(0.15, 0.65, 0.40, 0.85);

        script.cd();
        auto ela_yields = static_cast<TH1*>(alg_vz->ang_dep_sn_in[i][c]->Clone());
        ela_yields->SetTitle("Yields w/o C BG");
        ela_yields->Add(alg_carbon->yields_in[i][c], -1);
        ela_yields->Draw();
        draw_bin_sn(alg_vz->bin_center_theta(i),
                    ela_yields,
                    alg_vz->cuts_range()[c]);
      } else { // he4
        script.cd();
        auto s2 = s13::ana::draw_thstack({alg_vz->ang_dep_sn_in[i][c],
              alg_vz->ang_dep_outphi[i][c]},
          "Aft. cuts vs. OUT_phi BG", "NOSTACK");
        s2->SetMaximum(alg_vz->ang_dep_sn_in[i][c]->GetMaximum());
        gPad->BuildLegend(0.15, 0.65, 0.40, 0.85);

        script.cd();
        auto ela_yields = static_cast<TH1*>(alg_vz->ang_dep_sn_in[i][c]->Clone());
        ela_yields->SetTitle("Yields w/o OUT_phi BG");
        ela_yields->Add(alg_vz->ang_dep_outphi[i][c], -1);
        ela_yields->Draw();
        draw_bin_sn(alg_vz->bin_center_theta(i),
                    ela_yields,
                    alg_vz->cuts_range()[c]);
      }
    }
  }
}

std::vector<double>
draw_angdep_sn(s13::ana::RootScript& script, std::shared_ptr<VertzAngDep> alg) {
  std::vector<double> sn_ratios;
  for (size_t i = 0; i < alg->vertz_hists.size(); ++i) {
    if (i % 4 == 0) {
      script.NewPage(4, 2).cd(s13::ana::PadSeq::column);
    }

    script.cd();
    alg->vertz_hists[i]->Draw();
    gPad->SetLogy();
    draw_cut_marker_lines(alg->cuts_range()[0].espri_vertz_width, alg->aang_hists[i]);
    double vz_sn = compute_vertz_bin_sn(alg->bin_center_theta(i),
                                        alg->vertz_hists[i], alg->cuts_range()[0]);
    sn_ratios.push_back(vz_sn);
    s13::ana::draw_text(s13::ana::tstrfmt("IN/OUT: %.2f", vz_sn));

    script.cd();
    auto s1 = s13::ana::draw_thstack({alg->vertz_sn_in_hists[i],
          alg->vertz_sn_out_hists[i], alg->vertz_sn_tot_hists[i]},
      "1-d kin. corrs.", "NOSTACK");
    s1->SetMaximum(s1->GetMaximum() * 0.8);
    gPad->BuildLegend(0.15, 0.65, 0.40, 0.85);
    double sn = alg->vertz_sn_in_hists[i]->GetEntries()
      / alg->vertz_sn_out_hists[i]->GetEntries();
    double sel_ratio = alg->vertz_sn_in_hists[i]->GetEntries()
      / alg->vertz_sn_tot_hists[i]->GetEntries() * 100;
    s13::ana::draw_text(s13::ana::tstrfmt("IN/OUT: %.2f", sn),
                        0.6, 0.8, 0.9, 0.9);
    s13::ana::draw_text(s13::ana::tstrfmt("IN/TOT: %.2f%%", sel_ratio),
                        0.55, 0.7, 0.9, 0.8);

    // script.cd();
    // alg->vertz_sn_in_hists[i]->Draw();
    // script.cd();
    // alg->vertz_sn_out_hists[i]->Draw();
  }

  return sn_ratios;
}

void
draw_sn_graph(s13::ana::RootScript& script,
              std::shared_ptr<VertzAngDep> alg,
              std::vector<double> sn_ratios_data,
              DsType ds_type) {
  Yield sn_ratios(alg->bin_num(), s13::gk_cs_range_start, s13::gk_cs_range_end);
  Yield inv_sn_ratios(alg->bin_num(), s13::gk_cs_range_start, s13::gk_cs_range_end);
  for (size_t i = 0; i < alg->bin_num(); i++) {
    auto arg = inv_sn_ratios.GetBinArg(i);
    sn_ratios.SetBinValue(arg, sn_ratios_data[i]);
    sn_ratios.SetBinError(arg, 0);
    inv_sn_ratios.SetBinValue(arg, 1 / sn_ratios_data[i]);
    inv_sn_ratios.SetBinError(arg, 0);
  }
  sn_ratios = s13::ana::lab_yield_to_cm(sn_ratios, ds_type);
  inv_sn_ratios = s13::ana::lab_yield_to_cm(inv_sn_ratios, ds_type);

  auto title = s13::ana::tstrfmt("Bin data: (S/N)");
  auto sn_graph = s13::ana::
    BuildTGraphErrors(sn_ratios, title, "CM angle [deg]", "(S/N)");
  script.NewPage(1,1).cd();
  sn_graph->Draw("ACP");
  sn_graph->SetMinimum(0);

  title = s13::ana::tstrfmt("Bin data: 1 / (S/N)");
  auto inv_sn_graph = s13::ana::
    BuildTGraphErrors(inv_sn_ratios, title, "CM angle [deg]", "1 / (S/N)");
  script.NewPage(1,1).cd();
  inv_sn_graph->Draw("ACP");
  inv_sn_graph->SetMinimum(0);
}

std::shared_ptr<CarbonYields>
compute_carbon_yields(const std::vector<ElaCuts>& cuts_range) {
  auto num_runs = cli_opts.dataset_size();
  int carbon_num_runs = s13::gk_carbon_num_of_runs;
  if (num_runs != -1) {
    carbon_num_runs = 0.1971 * num_runs;
    if (carbon_num_runs == 0) carbon_num_runs = 1;
  }

  auto carbon_ds = init_carbon_segmented_dataset(1, carbon_num_runs);
  auto carbon_yields_alg =
    s13::ana::walk_alg2<CarbonYields>(carbon_ds,
                                      cli_opts.use_mt(),
                                      cuts_range,
                                      &ElaCuts::IsEspriVertzCut);
  return carbon_yields_alg;
}

void vertz_ang_dep() {
  auto cuts_range = s13::io::parse_cuts_range_config(std::fstream(cli_opts.cuts_config()));
  auto vz_alg =
    s13::ana::walk_alg2<VertzAngDep>(init_dataset_from_cli(cli_opts),
                                     cli_opts.use_mt(),
                                     cuts_range);

  auto script_fname =
    s13::ana::tstrfmt("%s_%s_vertz_ang_dep", cli_opts.output_file_basename().data(),
                      s13::ana::espri_selector_to_str(cuts_range[0].espri_selector).data());
  s13::ana::RootScript script(script_fname);
  gStyle->SetOptStat(1111111);
  gStyle->SetOptFit();
  // draw_angdep_fits(script, vz_alg);
  // auto sn_ratios = draw_angdep_sn(script, vz_alg);
  // draw_sn_graph(script, vz_alg, sn_ratios, cli_opts.dataset_type());
  if (cli_opts.dataset_type() == DsType::he6) {
    auto carbon_yields_alg = compute_carbon_yields(cuts_range);
    draw_angdep_sn2(script, vz_alg, carbon_yields_alg);
  } else {
    draw_angdep_sn2(script, vz_alg, nullptr);
  }
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

  vertz_ang_dep();
  std::cout << "exiting main\n" << std::flush;
}
