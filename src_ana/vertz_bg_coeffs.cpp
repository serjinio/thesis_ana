

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


using ElaCuts = s13::ana::ElasticScatteringCuts;
using Evt = s13::ana::ScatteringEvent;
using FitParams = s13::ana::FitParams;
using Espri = s13::ana::Espri;
using LinearInterp = s13::ana::LinearInterp;


class VertzBgCoeffs : public s13::ana::SimpleTTreeAlgorithmBase {

public:

  static constexpr double k_min_theta = 53;
  static constexpr double k_max_theta = 72;

  VertzBgCoeffs(const ElaCuts& cuts) :
    logger_{"VertzBgCoeffs"}, cuts_{cuts}, bin_num_{8} {
      init_hists();
    }

  VertzBgCoeffs(const VertzBgCoeffs&) = delete;
  VertzBgCoeffs(VertzBgCoeffs&&) = delete;
  VertzBgCoeffs operator=(const VertzBgCoeffs&) = delete;
  VertzBgCoeffs operator=(VertzBgCoeffs&&) = delete;

  virtual void setup() {
  }

  virtual void process_event(Evt& evt) {
    if (cuts_.IsSignalEvent(evt)) {
      compute_bg_scaling_coeff(evt);
      compute_varcut_sn(evt);
    }
  }

  virtual void finalize() {
    // compute_bg_scaling_factors();
  }

  ElaCuts cuts() {
    return cuts_;
  }

  ElaCuts cuts_no_aang() {
    auto tmp = cuts_;
    tmp.apply_espri_aang_cut = false;
    return tmp;
  }

  size_t bin_num() {
    return bin_num_;
  }


  std::vector<double> aang_cut_widths = {96 * 1/3., 96 * 2/3., 96, 96 * 2};
  std::vector<double> bg_signal_factors;
  std::vector<TH1*> vertz_cut_dists;
  std::vector<TH1*> vertz_aang_dists;
  std::vector<TH1*> vertz_in_cut;
  std::vector<TH1*> vertz_in_cut_aang_corr;
  std::vector<TH1*> vertz_out_cut;
  std::vector<TH1*> vertz_total;
  std::vector<TH1*> zpos_dists;
  std::vector<TH1*> zpos_in_cut_dists;
  std::vector<TH2*> kin_corrs_total_2d;
  std::vector<TH2*> kin_corrs_in_cut_2d;

  std::vector<TH1*> aang_dists_varcut;
  std::vector<TH1*> aang_dists_varcut_sn_in;
  std::vector<TH1*> aang_dists_varcut_sn_out;
  std::vector<TH1*> aang_dists_varcut_sn_tot;

  TH1* aang_dist;
  TH1* corr_aang_dist;
  TH1* rdc_zpos_dist;
  TH1* corr_zpos_dist;

private:

  TH1* init_vertz_cut_dist_hist(double min_theta, double max_theta,
                                double xlim = 350) {
    auto th1 = &s13::ana::make_th1;
    auto tit = s13::ana::tstrfmt("RDC aang difference (%.2f<#theta<%.2f)",
                                 min_theta, max_theta);
    auto hist = th1(200, -xlim, xlim, tit, "RDC (aang - aang_{tgt}) [mrad]", "Counts");
    return hist;
  }

  TH1* init_aang_dist_hist(TString title = "", double xlim = 1000) {
    auto th1 = &s13::ana::make_th1;
    auto tit = s13::ana::tstrfmt("#Delta aang. %s", title.Data());
    auto hist = th1(100, -xlim, xlim, tit, "RDC aang [mrad]", "Counts");
    return hist;
  }

  TH1* init_zpos_dist_hist(TString title = "Vert Z", double xlim = 400) {
    auto tit = s13::ana::tstrfmt("RDC vertZ distribution");
    auto th1 = &s13::ana::make_th1;
    auto hist = th1(200, -xlim, xlim, title, "Vert Z [mm]", "Counts");
    return hist;
  }

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

  TH2* init_xpos_vs_aang_hist(TString title) {
    auto th2 = &s13::ana::make_th2;
    auto hist = th2(300, -250, 250, 300, -350, 350, title,
                    "X [mm]", "aang [mrad]");
    return hist;
  }

  TH1* init_e_cut_dist_hist(double e_cut_width) {
    auto th1 = &s13::ana::make_th1;
    auto tit = s13::ana::tstrfmt("E cut dist. ( (E_{exp}-E_{sim})/E_{sim}= +/-%.1f)", e_cut_width/2);
    auto hist = th1(50, -1.5, 1.5, tit, "(E_{exp} - E_{sim}) / E_{sim}", "Counts");
    return hist;
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

  FitParams fit_kin_corr_peak(TH1* hist, double fit_region_width = 4) {
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

    double fit_rgn_hw = fit_region_width / 2;
    hist->Fit("gaus", "", "", -fit_rgn_hw, fit_rgn_hw);
    fit_params_raw = hist->GetFunction("gaus")->GetParameters();
    fit_errors_raw = hist->GetFunction("gaus")->GetParErrors();

    fit_params.SetParams(fit_params_raw);
    fit_params.SetErrors(fit_errors_raw);
    return fit_params;
  }

  void compute_bg_scaling_factors() {
    for (size_t i = 0; i < bin_num_; i++) {
      FitParams fit_bg = fit_kin_corr_peak(vertz_out_cut.at(i));
      FitParams fit_signal = fit_kin_corr_peak(vertz_in_cut.at(i));
      double scaling_fact = fit_bg.params[0] / fit_signal.params[0];
      bg_signal_factors.push_back(scaling_fact);
    }

    for (size_t i = 0; i < bin_num_; i++) {
      logger_.info("Scaling factor for theta [%.2f, %.2f]: %.2f",
                   bin_min_theta(i), bin_max_theta(i), bg_signal_factors[i]);
    }
  }

  void init_hists() {
    for (size_t i = 0; i < bin_num_; i++) {
      vertz_cut_dists.push_back(init_vertz_cut_dist_hist(bin_min_theta(i),
                                                         bin_max_theta(i), 300));
      vertz_aang_dists.push_back(init_aang_dist_hist());
      zpos_dists.push_back(init_zpos_dist_hist("Total"));
      zpos_in_cut_dists.push_back(init_zpos_dist_hist("In cut"));
      vertz_in_cut.push_back(init_kin_corr_hist("In cut"));
      vertz_in_cut_aang_corr.push_back(init_kin_corr_hist("In cut (corr. aang)"));
      vertz_out_cut.push_back(init_kin_corr_hist("Out cut"));
      vertz_total.push_back(init_kin_corr_hist("Total"));
      kin_corrs_total_2d.push_back(init_2d_kin_corr_hist("Total"));
      kin_corrs_in_cut_2d.push_back(init_2d_kin_corr_hist("In cut"));
    }

    for (size_t i = 0; i < aang_cut_widths.size(); ++i) {
      TString tit = s13::ana::tstrfmt("W: %.2f", aang_cut_widths[i]);
      // aang_dists_varcut.push_back(init_aang_dist_hist(tit, 300));
      aang_dists_varcut.push_back(init_zpos_dist_hist(tit, 400));
      aang_dists_varcut_sn_in.push_back(init_kin_corr_hist("IN"));
      aang_dists_varcut_sn_out.push_back(init_kin_corr_hist("OUT"));
      aang_dists_varcut_sn_tot.push_back(init_kin_corr_hist("TOT"));
    }

    aang_dist = init_aang_dist_hist("Total.", 800);
    corr_aang_dist = init_aang_dist_hist("Total. aang corrected.", 800);
    rdc_zpos_dist = init_zpos_dist_hist("Total");
    corr_zpos_dist = init_zpos_dist_hist("Total (aang corr.)");
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

  void compute_bg_scaling_coeff(Evt& evt) {
    if (!is_well_defined(evt.p_theta)) {
      return;
    }
    if (!(is_well_defined(evt.he_theta))) {
      return;
    }

    Espri espri = cuts_.IsLeftEspriEvent(evt) ? Espri::left : Espri::right;
    double aang = cuts_.IsLeftEspriEvent(evt) ? evt.esl_aang : evt.esr_aang;
    double aang_rdc = aang;
    // use corrected aang - for visual repr only, if not comment this out
    aang = correct_aang(aang, evt.p_theta_eff, espri);

    double xpos = cuts_.IsLeftEspriEvent(evt) ? evt.esl_xpos : evt.esr_xpos;
    double theta_diff = evt.he_theta - evt.he_theta_theor;
    int evt_bin_idx =  bin_idx(evt);

    if (std::abs(xpos) > 9000) {
      return;
    }
    if (( cuts_.IsLeftEspriEvent(evt) && evt.esl_de_raw < 240 ) ||
        ( cuts_.IsRightEspriEvent(evt) && evt.esr_de_raw < 240 )) {
      return;
    }
    if (evt_bin_idx == -1) {
      return;
    }

    double tgt_aang = std::atan(xpos/s13::gk_dist_rdc_target) * 1e3;
    double vertz = s13::ana::compute_vertz(xpos, aang, espri);
    aang_dist->Fill(aang_rdc - tgt_aang);
    corr_aang_dist->Fill(aang - tgt_aang);
    vertz_cut_dists[evt_bin_idx]->Fill(aang - tgt_aang);
    zpos_dists[evt_bin_idx]->Fill(vertz);
    vertz_aang_dists[evt_bin_idx]->Fill(aang);
    vertz_total[evt_bin_idx]->Fill(theta_diff);
    kin_corrs_total_2d[evt_bin_idx]->Fill(evt.he_theta, evt.p_theta_eff);

    auto loc_cuts = cuts_;
    loc_cuts.use_espri_aang_corr = true;
    if (loc_cuts.IsEspriAangCut(evt)) {
      vertz_in_cut_aang_corr[evt_bin_idx]->Fill(theta_diff);
    }
    loc_cuts.use_espri_aang_corr = false;
    if (loc_cuts.IsEspriAangCut(evt)) {
    // if (vertz > -250 && vertz < 250) {
      vertz_in_cut[evt_bin_idx]->Fill(theta_diff);
      zpos_in_cut_dists[evt_bin_idx]->Fill(vertz);
      kin_corrs_in_cut_2d[evt_bin_idx]->Fill(evt.he_theta, evt.p_theta_eff);
    } else {
      vertz_out_cut[evt_bin_idx]->Fill(theta_diff);
    }
  }

  void compute_varcut_sn(Evt& evt) {
    if (!is_well_defined(evt.p_theta)) {
      return;
    }
    if (!(is_well_defined(evt.he_theta))) {
      return;
    }

    Espri espri = cuts_.IsLeftEspriEvent(evt) ? Espri::left : Espri::right;
    double rdc_aang = cuts_.IsLeftEspriEvent(evt) ? evt.esl_aang : evt.esr_aang;
    double corr_aang = correct_aang(rdc_aang, evt.p_theta_eff, espri);
    // double aang = cuts_.use_espri_aang_corr ? corr_aang : rdc_aang;
    double xpos = cuts_.IsLeftEspriEvent(evt) ? evt.esl_xpos : evt.esr_xpos;
    double theta_diff = evt.he_theta - evt.he_theta_theor;

    if (std::abs(xpos) > 9000) {
      return;
    }
    if (( cuts_.IsLeftEspriEvent(evt) && evt.esl_de_raw < 240 ) ||
        ( cuts_.IsRightEspriEvent(evt) && evt.esr_de_raw < 240 )) {
      return;
    }

    // auto loc_cuts = cuts_;
    // if (loc_cuts.IsThetaCut(evt)) {
    //   return;
    // }

    // double tgt_aang = std::atan(xpos/s13::gk_dist_rdc_target) * 1e3;
    // double aang_diff = aang - tgt_aang;
    double corr_vertz = s13::ana::compute_vertz(xpos, corr_aang, espri);
    double rdc_vertz = s13::ana::compute_vertz(xpos, rdc_aang, espri);
    double vertz = cuts_.use_espri_aang_corr ? corr_vertz : rdc_vertz;

    for (size_t i = 0; i < aang_cut_widths.size(); ++i) {
      aang_dists_varcut[i]->Fill(vertz);
      aang_dists_varcut_sn_tot[i]->Fill(theta_diff);
      cuts_.espri_aang_width = aang_cut_widths[i];
      // if (cuts_.IsEspriAangCut(evt)) {
      if (std::abs(vertz) < aang_cut_widths[i]) {
        aang_dists_varcut_sn_in[i]->Fill(theta_diff);
      } else {
        aang_dists_varcut_sn_out[i]->Fill(theta_diff);
      }
    }

    corr_zpos_dist->Fill(corr_vertz);
    rdc_zpos_dist->Fill(rdc_vertz);

  }

  bool is_hodf_z2(Evt& evt) {
    double hodf_q = cuts_.get_hodf_q(evt);
    return hodf_q > 7;
  }

  bool is_hodf_he6(Evt& evt) {
    return cuts_.IsHodfHe6Cut(evt);
  }

  bool is_hodf_he4(Evt& evt) {
    return cuts_.get_hodf_hit_id(evt) >= 16 && cuts_.get_hodf_q(evt);
  }


  s13::misc::MessageLogger logger_;

  ElaCuts cuts_;

  size_t bin_num_;

};

void draw_aang_cut_marker_lines(double cut_width, double offset, double ymax) {
  s13::ana::draw_marker_line(cut_width / 2 + offset, 0,
                             cut_width / 2 + offset,
                             ymax);
  s13::ana::draw_marker_line(-cut_width / 2 + offset, 0,
                             -cut_width / 2 + offset,
                             ymax);
}

void draw_vertz_dists(s13::ana::RootScript& script, std::shared_ptr<VertzBgCoeffs> alg) {
  double aang_offset = alg->cuts().espri_selector == s13::ana::Espri::right ?
    ElaCuts::k_esr_aang_center_offset : ElaCuts::k_esl_aang_center_offset;
  double aang_hw = alg->cuts().espri_aang_width / 2;

  for (size_t i = 0; i < alg->bin_num(); i++) {
    if (i % 4 == 0) {
      script.NewPage(4, 3).cd(s13::ana::PadSeq::column);
    }

    script.cd();
    alg->vertz_cut_dists[i]->Draw("HIST");
    alg->vertz_cut_dists[i]->SetMaximum(alg->vertz_cut_dists[i]->GetMaximum() * 1.5);
    gPad->SetLogy();
    s13::ana::draw_marker_line(aang_hw + aang_offset, 0,
                               aang_hw + aang_offset,
                               alg->vertz_cut_dists[i]->GetMaximum());
    s13::ana::draw_marker_line(-aang_hw + aang_offset, 0,
                               -aang_hw + aang_offset,
                               alg->vertz_cut_dists[i]->GetMaximum());

    script.cd();
    auto s2 = s13::ana::draw_thstack({alg->zpos_dists[i],
          alg->zpos_in_cut_dists[i]}, "Vert Z dists.", "NOSTACK HIST");
    gPad->SetLogy();
    s2->SetMaximum(s2->GetMaximum() * 1.0);
    gPad->BuildLegend(0.15, 0.65, 0.40, 0.85);
    // alg->zpos_dists[i]->Draw();

    script.cd();
    s2 = s13::ana::draw_thstack({alg->vertz_in_cut[i],
          alg->vertz_out_cut[i]}, "1-d kin. corrs.", "NOSTACK HIST");
    s2->SetMaximum(s2->GetMaximum() * 1.0);
    gPad->BuildLegend(0.15, 0.65, 0.40, 0.85);

  }
}

void draw_bins(s13::ana::RootScript& script, std::shared_ptr<VertzBgCoeffs> alg) {
  double aang_offset = 0;
  double aang_hw = alg->cuts().espri_aang_width / 2;
  std::vector<double> aang_means;
  std::vector<double> aang_sigmas;

  for (size_t i = 0; i < alg->bin_num(); i++) {
    if (i % 4 == 0) {
      script.NewPage(4, 3).cd(s13::ana::PadSeq::column);
    }

    script.cd();
    alg->vertz_cut_dists[i]->Fit("gaus", "", "", -100, 100);
    auto fit_params = alg->vertz_cut_dists[i]->GetFunction("gaus")->GetParameters();
    aang_means.push_back(fit_params[1]);
    aang_sigmas.push_back(fit_params[2]);
    alg->vertz_cut_dists[i]->Draw();
    alg->vertz_cut_dists[i]->SetMaximum(alg->vertz_cut_dists[i]->GetMaximum() * 1.5);
    s13::ana::draw_marker_line(aang_hw + aang_offset, 0,
                               aang_hw + aang_offset,
                               alg->vertz_cut_dists[i]->GetMaximum());
    s13::ana::draw_marker_line(-aang_hw + aang_offset, 0,
                               -aang_hw + aang_offset,
                               alg->vertz_cut_dists[i]->GetMaximum());
    // s13::ana::draw_text(alg->cuts().GetTotalCutTitle(), 0.15, 0.75, 0.7, 0.9);
    gPad->SetLogy();

    // script.cd();
    // alg->vertz_aang_dists[i]->Draw();
    // gPad->SetLogy();

    script.cd();
    alg->vertz_in_cut[i]->Draw();
    double in_cut_fraction = alg->vertz_in_cut[i]->Integral()
      / alg->vertz_total[i]->Integral() * 100;
    s13::ana::draw_text(s13::ana::tstrfmt("IN/TOT: %.2f%%", in_cut_fraction));

    script.cd();
    alg->vertz_in_cut_aang_corr[i]->Draw();
    in_cut_fraction = alg->vertz_in_cut_aang_corr[i]->Integral()
      / alg->vertz_total[i]->Integral() * 100;
    s13::ana::draw_text(s13::ana::tstrfmt("IN/TOT: %.2f%%", in_cut_fraction));

    // auto s2 = s13::ana::draw_thstack({alg->vertz_in_cut[i],
    //       alg->vertz_out_cut[i]}, "1-d kin. corrs.", "NOSTACK");
    // s2->SetMaximum(s2->GetMaximum() * 0.8);
    // gPad->BuildLegend(0.15, 0.65, 0.40, 0.85);

    // script.cd();
    // auto in_cut_scaled = (TH1*)alg->vertz_in_cut[i]->Clone();
    // in_cut_scaled->SetTitle("In (scaled)");
    // in_cut_scaled->Scale(alg->bg_signal_factors[i]);
    // auto s3 = s13::ana::draw_thstack({in_cut_scaled,
    //       alg->vertz_out_cut[i]}, "1-d kin. corrs.", "NOSTACK HIST");
    // s3->SetMaximum(s3->GetMaximum() * 0.8);
    // gPad->BuildLegend(0.15, 0.65, 0.40, 0.85);

    // script.cd();
    // auto bg_dist = (TH1*)alg->vertz_out_cut[i]->Clone();
    // bg_dist->SetTitle("Clean BG");
    // bg_dist->Add(in_cut_scaled, -1);
    // bg_dist->Draw("HIST");
    // gPad->BuildLegend(0.15, 0.65, 0.40, 0.85);

    // auto kin_graph = s13::ana::BuildHe6ProtonKinematicsGraph();
    // kin_graph->SetLineWidth(2);
    // script.cd();
    // alg->kin_corrs_total_2d[i]->Draw("col");
    // kin_graph->Draw("SAME");
    // script.cd();
    // alg->kin_corrs_in_cut_2d[i]->Draw("col");
    // kin_graph->Draw("SAME");
  }

  std::cout << "Aang means: " << std::endl;
  for (auto el : aang_means) {
    std::cout << el << std::endl;
  }

  std::cout << "Aang sigmas: " << std::endl;
  for (auto el : aang_sigmas) {
    std::cout << el << std::endl;
  }

  script.NewPage(2,1).cd();
  alg->aang_dist->Fit("gaus", "", "", -100, 100);
  alg->aang_dist->Draw();
  gPad->SetLogy();

  script.cd();
  alg->corr_aang_dist->Fit("gaus", "", "", -100, 100);
  alg->corr_aang_dist->Draw();
  gPad->SetLogy();
}

void draw_aang_eff(s13::ana::RootScript& script, std::shared_ptr<VertzBgCoeffs> alg) {
  auto hist = alg->corr_zpos_dist;
  script.NewPage(1,1).cd();
  hist->Draw();
  hist->Fit("gaus", "", "", -40, 40);
  // double *fit_params = hist->GetFunction("gaus")->GetParameters();
  // double sigma = fit_params[2];
  // double vertz_hw = sigma * 3;
  double vertz_hw = 150;
  s13::ana::draw_marker_line(vertz_hw, 0, vertz_hw,
                             hist->GetMaximum());
  s13::ana::draw_marker_line(-vertz_hw, 0, -vertz_hw,
                             hist->GetMaximum());

  auto bin_num = [&hist](double x) -> int { return hist->GetXaxis()->FindBin(x); };
  std::cout << "min vertz: " << -vertz_hw << "; max vertz: "
            << vertz_hw << std::endl;
  std::cout << "min bin: " << bin_num(-vertz_hw) << "; max bin: "
            << bin_num(vertz_hw) << std::endl;

  double in_cut_int = hist->Integral(bin_num(-vertz_hw), bin_num(vertz_hw), "");
  std::cout << "IN cut integrand: " << in_cut_int << std::endl;
  double out_cut_int = hist->GetEntries() - in_cut_int;
  double aang_eff = out_cut_int / hist->GetEntries();
  s13::ana::draw_text(s13::ana::tstrfmt("aang_{eff}: %.2f%%", (1 - aang_eff) * 100));

  // gPad->SetLogy();
}

void draw_varcut(s13::ana::RootScript& script, std::shared_ptr<VertzBgCoeffs> alg) {
  for (size_t i = 0; i < alg->aang_cut_widths.size(); ++i) {
    if (i % 4 == 0) {
      script.NewPage(4, 3).cd(s13::ana::PadSeq::column);
    }

    double aang_hw = alg->aang_cut_widths[i] / 2;

    script.cd();
    alg->aang_dists_varcut[i]->Draw();
    gPad->SetLogy();
    alg->aang_dists_varcut[i]->SetMaximum(alg->aang_dists_varcut[i]->GetMaximum() * 1.5);
    s13::ana::draw_marker_line(aang_hw, 0, aang_hw,
                               alg->aang_dists_varcut[i]->GetMaximum());
    s13::ana::draw_marker_line(-aang_hw, 0, -aang_hw,
                               alg->aang_dists_varcut[i]->GetMaximum());

    script.cd();
    alg->aang_dists_varcut_sn_in[i]->Draw();
    double in_cut_fraction = alg->aang_dists_varcut_sn_in[i]->Integral()
      / alg->aang_dists_varcut_sn_tot[i]->Integral() * 100;
    s13::ana::draw_text(s13::ana::tstrfmt("IN/TOT: %.2f%%", in_cut_fraction));

    script.cd();
    auto s1 = s13::ana::draw_thstack({alg->aang_dists_varcut_sn_in[i],
          alg->aang_dists_varcut_sn_out[i], alg->aang_dists_varcut_sn_tot[i]},
      "1-d kin. corrs.", "NOSTACK");
    s1->SetMaximum(s1->GetMaximum() * 0.8);
    gPad->BuildLegend(0.15, 0.65, 0.40, 0.85);
  }

  // script.NewPage(2,1).cd();
  // // alg->rdc_zpos_dist->Fit("gaus", "", "", -150, 150);
  // alg->rdc_zpos_dist->Draw();
  // gPad->SetLogy();

  // script.cd();
  // alg->corr_zpos_dist->Draw();
  // // alg->corr_zpos_dist->Fit("gaus", "", "", -150, 150);
  // gPad->SetLogy();

  draw_aang_eff(script, alg);
}

void vertz_bg_coeffs() {
  auto cuts = s13::io::parse_cuts_config(std::fstream(cli_opts.cuts_config()));
  cuts.use_espri_aang_corr = true;
  auto alg =
    s13::ana::walk_alg<VertzBgCoeffs>(cuts,
                                      init_dataset_from_cli(cli_opts),
                                      cli_opts.use_mt());

  s13::ana::RootScript script(cli_opts.output_file_basename());
  gStyle->SetOptStat(1111111);
  gStyle->SetOptFit();
  // draw_bins(script, alg);
  draw_varcut(script, alg);
  // draw_vertz_dists(script, alg);
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

  vertz_bg_coeffs();
  std::cout << "exiting main\n" << std::flush;
}
